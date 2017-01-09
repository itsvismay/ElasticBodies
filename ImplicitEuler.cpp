#include "ImplicitEuler.h"

using namespace Eigen;
using namespace std;

typedef Eigen::Triplet<double> Trip;
typedef Matrix<double, 12, 1> Vector12d;

void ImplicitEuler::findgBlock(VectorXd& g_block, VectorXd& x, VectorXd& x_old, int ignorePast){
	//VectorXd g = (RegMass*x - RegMass*x_old)/h - RegMass*v_old - h*f;
	VectorXd g = (RegMass*x - RegMass*x_old) - (RegMass*v_old*h) - h*h*f;
	g_block = g.head(ignorePast*3);

}

void ImplicitEuler::renderNewtonsMethod(VectorXd& ext_force){
	//Implicit Code
	v_k.setZero();
	x_k.setZero();
	x_k = x_old;
	x_k(0) += 0.01;//2.1898e-09;
	v_k = v_old;
	clock_t t0 = clock();

	int ignorePastIndex = TV.rows() - fixedVerts.size();
	SparseMatrix<double> forceGradientStaticBlock;
	forceGradientStaticBlock.resize(3*ignorePastIndex, 3*ignorePastIndex);
	clock_t t1 = clock();


	SparseMatrix<double> RegMassBlock;
	RegMassBlock.resize(3*ignorePastIndex, 3*ignorePastIndex);
	RegMassBlock = RegMass.block(0, 0, 3*ignorePastIndex, 3*ignorePastIndex);
	clock_t t2 = clock();
	cout<<"t2 t1"<<endl;
	cout<<"SECONDS"<<double(t2 - t1)/CLOCKS_PER_SEC<<endl;

	forceGradient.setZero();
	bool Nan=false;
	int NEWTON_MAX = 100, i =0;
	// cout<<"--------"<<simTime<<"-------"<<endl;
	// cout<<"x_k"<<endl;
	// cout<<x_k<<endl<<endl;
	// cout<<"v_k"<<endl;
	// cout<<v_k<<endl<<endl;
	// cout<<"--------------------"<<endl;
	for( i=0; i<NEWTON_MAX; i++){
		clock_t t3 = clock();
		grad_g.setZero();
		ImplicitXtoTV(x_k, TVk);//TVk value changed in function
		ImplicitCalculateElasticForceGradient(TVk, forceGradient);
		ImplicitCalculateForces(TVk, forceGradient, x_k, f);

		// VectorXd g_block = x_k - x_old -h*v_old -h*h*InvMass*f;
		// grad_g = Ident - h*h*InvMass*forceGradient - h*rayleighCoeff*InvMass*forceGradient;


		//Block forceGrad and f to exclude the fixed verts
		forceGradientStaticBlock = forceGradient.block(0,0, 3*(ignorePastIndex), 3*ignorePastIndex);
		VectorXd g_block;
		findgBlock(g_block, x_k, x_old, ignorePastIndex);
		// VectorXd g = RegMass*x_k - RegMass*x_old - h*RegMass*v_old - h*h*f;
		// VectorXd g_block = g.head(ignorePastIndex*3);
		grad_g = RegMassBlock - h*h*forceGradientStaticBlock - h*rayleighCoeff*forceGradientStaticBlock;

		// Sparse Cholesky LL^T
		if(llt_solver.info() == Eigen::NumericalIssue){
			cout<<"Possibly using a non- pos def matrix in the LLT method"<<endl;
			exit(0);
		}


		llt_solver.factorize(grad_g);
		VectorXd deltaX = -1* llt_solver.solve(g_block);
		x_k.segment(0, 3*(ignorePastIndex)) += deltaX;

		//Sparse QR
		// SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> sqr;
		// sqr.compute(grad_g);
		// VectorXd deltaX = -1*sqr.solve(g_block);
		// x_k.segment(0, 3*(ignorePastIndex)) += deltaX;

		if(x_k != x_k){
			Nan = true;
			break;
		}

		cout<<"g norm"<<endl;
		cout<<g_block.squaredNorm()/convergence_scaling_paramter<<endl;
		if(g_block.squaredNorm()/convergence_scaling_paramter < sqrt(1e-11)/sqrt(convergence_scaling_paramter)){
			break;
		}

	}
	if(Nan){
		cout<<"ERROR: Newton's method doesn't converge"<<endl;
		cout<<i<<endl;
		exit(0);
	}
	if(i== NEWTON_MAX){
		cout<<"ERROR: Newton max reached"<<endl;
		cout<<i<<endl;
		exit(0);
	}
	v_old = (x_k - x_old)/h;
	x_old = x_k;
}

void ImplicitEuler::ImplicitCalculateForces( MatrixXd& TVk, SparseMatrix<double>& forceGradient, VectorXd& x_k, VectorXd& f){
	// //gravity
	f.setZero();

	for(unsigned int i=0; i<f.size()/3; i++){
		f(3*i+1) += massVector(3*i+1)*gravity;
	}

	//elastic
	if(solver.compare("newton")==0){
		for(unsigned int i=0; i<M.tets.size(); i++){
			Vector4i indices = M.tets[i].verticesIndex;
			MatrixXd F_tet = M.tets[i].oldComputeElasticForces(TVk, simTime%2);
			f.segment<3>(3*indices(0)) += F_tet.col(0);
			f.segment<3>(3*indices(1)) += F_tet.col(1);
			f.segment<3>(3*indices(2)) += F_tet.col(2);
			f.segment<3>(3*indices(3)) += F_tet.col(3);
		}
	}else{
		for(unsigned int i=0; i<M.tets.size(); i++){
			M.tets[i].computeElasticForces(TVk, f);
		}
	}


	// for(int k=0; k<f.rows(); k++){
	// 	if(fabs(external_f(k))>0.0001){
	// 		f(k) += external_f(k);
	// 	}
	// }
	// f += rayleighCoeff*forceGradient*(x_k - x_old)/h;
	return;
}

void ImplicitEuler::ImplicitCalculateElasticForceGradient(MatrixXd& TVk, SparseMatrix<double>& forceGradient){
	forceGradient.setZero();

	vector<Trip> triplets1;
	triplets1.reserve(12*12*M.tets.size());
	for(unsigned int i=0; i<M.tets.size(); i++){
		//Get P(dxn), dx = [1,0, 0...], then [0,1,0,....], and so on... for all 4 vert's x, y, z
		//P is the compute Force Differentials blackbox fxn

		Vector12d dx(12);
		dx.setZero();
		Vector4i indices = M.tets[i].verticesIndex;
		int kj;
		for(unsigned int j=0; j<12; j++){
			dx(j) = 1;
			MatrixXd dForces = M.tets[i].computeForceDifferentials(TVk, dx);
			kj = j%3;
			//row in order for dfxi/dxi ..dfxi/dzl
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[0], dForces(0,0)));
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[0]+1, dForces(1,0)));
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[0]+2, dForces(2,0)));

			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[1], dForces(0,1)));
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[1]+1, dForces(1,1)));
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[1]+2, dForces(2,1)));

			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[2], dForces(0,2)));
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[2]+1, dForces(1,2)));
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[2]+2, dForces(2,2)));

			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[3], dForces(0,3)));
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[3]+1, dForces(1,3)));
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[3]+2, dForces(2,3)));
			dx(j) = 0; //ASK check is this efficient?
		}
	}
	forceGradient.setFromTriplets(triplets1.begin(), triplets1.end());
	return;
}

void ImplicitEuler::render(VectorXd& ext_force){
	simTime+=1;
	cout<<"i"<<simTime<<endl;

	if(solver.compare("newton")==0){
		renderNewtonsMethod(ext_force);
	}else{
		cout<<"Solver not specified properly"<<endl;
		exit(0);
	}

	// cout<<"*******************"<<endl;
	// cout<< "New Pos"<<simTime<<endl;
	// cout<<x_old<<endl<<endl;
	// cout<< "New Vels"<<simTime<<endl;
	// cout<<v_old<<endl;
	// cout<<"*****************"<<endl<<endl;
	IntegratorAbstract::printInfo();

	ImplicitXtoTV(x_old, TV);

}


void ImplicitEuler::initializeIntegrator(double ph, SolidMesh& pM, MatrixXd& pTV, MatrixXi& pTT){
	IntegratorAbstract::initializeIntegrator(ph, pM, pTV, pTT);
	ZeroMatrix.resize(3*vertsNum, 3*vertsNum);
	ZeroMatrix.setZero();
	Ident.resize(3*vertsNum, 3*vertsNum);
	Ident.setIdentity();
	forceGradient.resize(3*vertsNum, 3*vertsNum);
	grad_g.resize(3*vertsNum, 3*vertsNum);
	x_k.resize(3*vertsNum);
	v_k.resize(3*vertsNum);
	external_f.resize(3*vertsNum);
	external_f.setZero();
	x_k.setZero();
	v_k.setZero();
	TVk = TV;
}

void ImplicitEuler::ImplicitXtoTV(VectorXd& x_tv, MatrixXd& TVk){
	TVk.setZero();
	for(unsigned int i=0; i < M.tets.size(); i++){
		Vector4i indices = M.tets[i].verticesIndex;
		TVk.row(indices(0)) = Vector3d(x_tv(3*indices(0)), x_tv(3*indices(0)+1), x_tv(3*indices(0) +2));
		TVk.row(indices(1)) = Vector3d(x_tv(3*indices(1)), x_tv(3*indices(1)+1), x_tv(3*indices(1) +2));
		TVk.row(indices(2)) = Vector3d(x_tv(3*indices(2)), x_tv(3*indices(2)+1), x_tv(3*indices(2) +2));
		TVk.row(indices(3)) = Vector3d(x_tv(3*indices(3)), x_tv(3*indices(3)+1), x_tv(3*indices(3) +2));
	}
	return;
}

void ImplicitEuler::ImplicitTVtoX(VectorXd& x_tv, MatrixXd& TVk){
	x_tv.setZero();
	for(unsigned int i = 0; i < M.tets.size(); i++){
		Vector4i indices = M.tets[i].verticesIndex;

		x_tv(3*indices(0)) = TVk.row(indices(0))[0];
		x_tv(3*indices(0)+1) = TVk.row(indices(0))[1];
		x_tv(3*indices(0)+2) = TVk.row(indices(0))[2];

		x_tv(3*indices(1)) = TVk.row(indices(1))[0];
		x_tv(3*indices(1)+1) = TVk.row(indices(1))[1];
		x_tv(3*indices(1)+2) = TVk.row(indices(1))[2];

		x_tv(3*indices(2)) = TVk.row(indices(2))[0];
		x_tv(3*indices(2)+1) = TVk.row(indices(2))[1];
		x_tv(3*indices(2)+2) = TVk.row(indices(2))[2];

		x_tv(3*indices(3)) = TVk.row(indices(3))[0];
		x_tv(3*indices(3)+1) = TVk.row(indices(3))[1];
		x_tv(3*indices(3)+2) = TVk.row(indices(3))[2];
	}
}
