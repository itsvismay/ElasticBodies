#include "ImplicitEuler.h"

using namespace Eigen;
using namespace std;

typedef Eigen::Triplet<double> Trip;
typedef Matrix<double, 12, 1> Vector12d;

void function1_grad(const real_1d_array &x, double &func, real_1d_array &grad, void *impe) 
{	
	ImplicitEuler* in = (ImplicitEuler*) impe;
	cout<<"-----------"<<endl;
    // this callback calculates f(x0,x1) = 100*(x0+3)^4 + (x1-3)^4
    // and its derivatives df/d0 and df/dx1
    //
    int n = 3*(in->vertsNum - in->fixedVerts.size());
    for(int i=0; i<n; i++){
    	in->x_k(i) = x[i];
    }
    // in->x_k(0)+=2.1898e-09;

    in->ImplicitXtoTV(in->x_k, in->TVk);//TVk value changed in function
	in->ImplicitCalculateElasticForceGradient(in->TVk, in->forceGradient); 
	in->ImplicitCalculateForces(in->TVk, in->forceGradient, in->x_k, in->f);
	VectorXd zero;
	zero.resize(n);
	zero.setZero();
	VectorXd g_block; 
	g_block.resize(n);
	g_block.setZero();
	in->findgBlock(g_block, in->x_k, in->x_old, n/3);

	cout<<"printing the middle val"<<endl;
	double fx0 = 0.5*in->x_k.transpose()*in->RegMass*in->x_k;
	double fx1 = in->x_old.transpose()*in->RegMass*in->x_k;
	double fx2 = in->h*in->v_old.transpose()*in->RegMass*in->x_k; 
	// cout<<fx0<<endl;
	// cout<<fx1<<endl;
	// cout<<fx2<<endl;
	func= fx0 - fx1 - fx2;  //big G function, anti-deriv of g
	// func = 0.5*in->h*((in->x_k - in->x_old)/in->h - in->v_old).transpose()*in->RegMass*((in->x_k - in->x_old)/in->h - in->v_old);
	double lastTerm = 0.0;

	for(int i=0; i<n; i++){
		// func += 1*(
		// 	(0.5*in->x_k(i)*in->massVector(i)*in->x_k(i) 
		// 	- in->massVector(i)*in->x_old(i)*in->x_k(i) 
		// 	- in->massVector(i)*in->h*in->v_old(i)*in->x_k(i));
		grad[i] = (1/in->convergence_scaling_paramter)*(g_block(i));
	}
	// cout<<lastTerm<<endl;
	// cout<<(lastTerm)-in->prevfx<<endl;
	// in->prevfx = lastTerm;
	cout<<func/in->convergence_scaling_paramter<<endl;
	double funcForce = 0;
	//force anti-derivative
	double strainE=0;
	for(unsigned int i=0; i<in->M.tets.size(); i++){
		strainE += in->M.tets[i].undeformedVol*in->M.tets[i].energyDensity;		
	}

	func+= in->h*(strainE);
	//damping anti-derivative
	// //func += in->h*rayleighCoeff*((in->x_k.dot(in->f) - strainE) - in->f.dot(in->x_old));
	//TODO Add gravity potential
	double grav =0;
	for(unsigned int i=0; i<in->x_k.size()/3; i++){
		 grav -= in->h*in->massVector(3*i+1)*in->x_k(3*i+1)*gravity;
		func -= in->h*in->massVector(3*i+1)*in->x_k(3*i+1)*gravity;
	}

	//This is the return value **
	func = (func/in->convergence_scaling_paramter);
	//****

	//******Lots of Logging***
    cout<<"n size"<<endl;
    cout<<n<<endl;
    cout<<"lbfgs x_k"<<endl;
	cout<<in->x_k<<endl;
	cout<<"x vals"<<endl;
	cout<<x.tostring(3)<<endl;
	cout<<"xold"<<endl;
	cout<<in->x_old<<endl;
    cout<<"lbfgs x_k- x_old"<<endl;
	cout<<in->x_k - in->x_old<<endl;
	cout<<"lbfgs force"<<endl;
	cout<<in->f<<endl;
	cout<<"lbfgs g"<<endl;
	double normgrad = 0;
	for(int i=0; i<n; i++){
		cout<<grad[i]<<endl;
		normgrad+=grad[i]*grad[i];
	}
	cout<<sqrt(normgrad)<<endl;
	cout<<"g norm/ scaling factor"<<endl;
	cout<<g_block.squaredNorm()/in->convergence_scaling_paramter<<endl;
	cout<<"Manually Computer G"<<endl;
	for(int i=0; i<1; i++){
		//RegMass*x - RegMass*x_old - h*RegMass*v_old - h*h*f;
		cout<<"Term 1: "<<(in->massVector(i)*in->x_k(i))/in->convergence_scaling_paramter<<endl;
		cout<<"Term 2: -"<< (in->massVector(i)*in->x_old(i))/in->convergence_scaling_paramter<<endl;
		cout<<"Term 3: -"<<(in->h*in->massVector(i)*in->v_old(i))/in->convergence_scaling_paramter<<endl;
		cout<<"term 4: -"<<(in->h*in->h*in->f(i))/in->convergence_scaling_paramter<<endl;
	}
	
	cout<<"func"<<endl;
	cout<<func<<endl;
	cout<<"grav"<<endl;
	cout<<grav<<endl;
	cout<<"Strain term"<<endl;
	cout<<in->h*in->h*strainE<<endl;

	//Using this to keep track of number of line search iterations
	in->simTime+=1;// not really "simulation time". 
}



void ImplicitEuler::findgBlock(VectorXd& g_block, VectorXd& x, VectorXd& x_old, int ignorePast){
	//VectorXd g = (RegMass*x - RegMass*x_old)/h - RegMass*v_old - h*f;
	VectorXd g = (RegMass*x - RegMass*x_old) - (RegMass*v_old*h) - h*h*f;
	g_block = g.head(ignorePast*3);
	cout<<"		First Term"<<endl;
	cout<<	RegMass*x<<endl;
	cout<<"--"<<endl;
	cout<<x <<endl;
	cout<<RegMass<<endl;
	cout<<"		Second Term"<<endl;
	cout<<	RegMass*x_old<<endl;
	cout<<"		Third Term"<<endl;
	cout<<	h*RegMass*v_old<<endl;
	cout<<"		Fourth Term"<<endl;
	cout<<	h*h*f<<endl;
}



int ImplicitEuler::alglibLBFGS(VectorXd& ext_force){
  	external_f = ext_force;
    int N = 3*(vertsNum);
    prevfx =0;
    // This example demonstrates minimization of f(x,y) = 100*(x+3)^4+(y-3)^4
    // using LBFGS method.
    //
    // Several advanced techniques are demonstrated:
    // * upper limit on step size
    // * restart from new point
    //
    real_1d_array x;
    double *positions= new double[N];
   	for(int i=0; i<N; i++){
   		positions[i] = x_old(i);
   	}

    x.setcontent(N, positions);
    double epsg = sqrt(1e-11)/sqrt(convergence_scaling_paramter);
    double epsf = 0;
    double epsx = 0;
    double stpmax = 0;
    ae_int_t maxits = 0;
    minlbfgsstate state;
    minlbfgsreport rep;

    // first run
    minlbfgscreate(1, x, state);
    minlbfgssetcond(state, epsg, epsf, epsx, maxits);
    // minlbfgssetstpmax(state, stpmax);
    alglib::minlbfgsoptimize(state, function1_grad,NULL, this);
    minlbfgsresults(state, x, rep);

    printf("TERMINATION TYPE: %d\n", int(rep.terminationtype)); // EXPECTED: 4
    cout<<epsg<<endl;
    for(int i=0; i<N; i++){
    	x_k(i) = x[i];
    }
    v_old = (x_k - x_old)/h;
    x_old = x_k;
    ImplicitXtoTV(x_old, TV);
    return 0;
}

void ImplicitEuler::renderNewtonsMethod(VectorXd& ext_force){
	//Implicit Code
	v_k.setZero();
	x_k.setZero();
	x_k = x_old;
	// x_k(0) += 0.0001;//2.1898e-09;
	v_k = v_old;

	int ignorePastIndex = TV.rows() - fixedVerts.size();
	SparseMatrix<double> forceGradientStaticBlock;
	forceGradientStaticBlock.resize(3*ignorePastIndex, 3*ignorePastIndex);	
	
	SparseMatrix<double> RegMassBlock;
	RegMassBlock.resize(3*ignorePastIndex, 3*ignorePastIndex);
	RegMassBlock = RegMass.block(0, 0, 3*ignorePastIndex, 3*ignorePastIndex);

	forceGradient.setZero();
	bool Nan=false;
	int NEWTON_MAX = 10, i =0;
	// cout<<"--------"<<simTime<<"-------"<<endl;
	// cout<<"x_k"<<endl;
	// cout<<x_k<<endl<<endl;
	// cout<<"v_k"<<endl;
	// cout<<v_k<<endl<<endl;
	// cout<<"--------------------"<<endl;
	for( i=0; i<NEWTON_MAX; i++){
		grad_g.setZero();
		ImplicitXtoTV(x_k, TVk);//TVk value changed in function
		ImplicitCalculateElasticForceGradient(TVk, forceGradient); 
		ImplicitCalculateForces(TVk, forceGradient, x_k, f);
		for(int k=0; k<f.rows(); k++){
			if(abs(ext_force(k))>0.0001){
				f(k) = 0.01*ext_force(k);
			}
		}
		// VectorXd g_block = x_k - x_old -h*v_old -h*h*InvMass*f;
		// grad_g = Ident - h*h*InvMass*forceGradient - h*rayleighCoeff*InvMass*forceGradient;
		

		//Block forceGrad and f to exclude the fixed verts
		forceGradientStaticBlock = forceGradient.block(0,0, 3*(ignorePastIndex), 3*ignorePastIndex);
		VectorXd g_block;
		findgBlock(g_block, x_k, x_old, ignorePastIndex);
		// VectorXd g = RegMass*x_k - RegMass*x_old - h*RegMass*v_old - h*h*f;
		// VectorXd g_block = g.head(ignorePastIndex*3);
		grad_g = RegMassBlock - h*h*forceGradientStaticBlock - h*rayleighCoeff*forceGradientStaticBlock;

		cout<<"newton f"<<endl;
		cout<<f<<endl;
		cout<<"newton g"<<endl;
		cout<<g_block/convergence_scaling_paramter<<endl;
		cout<<"newton x_k"<<endl;
		cout<<x_k<<endl;
	
		//solve for delta x
		// Conj Grad
		// ConjugateGradient<SparseMatrix<double>> cg;
		// cg.compute(grad_g);
		// VectorXd deltaX = -1*cg.solve(g);

		// Sparse Cholesky LL^T
		SimplicialLLT<SparseMatrix<double>> llt;
		llt.compute(grad_g);
		if(llt.info() == Eigen::NumericalIssue){
			cout<<"Possibly using a non- pos def matrix in the LLT method"<<endl;
			exit(0);
		}
		VectorXd deltaX = -1* llt.solve(g_block);
		x_k.segment(0, 3*(ignorePastIndex)) += deltaX;

		//Sparse QR 
		// SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> sqr;
		// sqr.compute(grad_g);
		// VectorXd deltaX = -1*sqr.solve(g_block);
		// x_k.segment(0, 3*(ignorePastIndex)) += deltaX;

		// CholmodSimplicialLLT<SparseMatrix<double>> cholmodllt;
		// cholmodllt.compute(grad_g);
		// VectorXd deltaX = -cholmodllt.solve(g_block);

		if(x_k != x_k){
			Nan = true;
			break;
		}

		if(g_block.squaredNorm()/convergence_scaling_paramter < 1e-11){
			cout<<"g norm"<<endl;
			cout<<g_block.squaredNorm()/convergence_scaling_paramter<<endl;
			break;
		}
		// exit(0);



		

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
	v_old.setZero();
	v_old = (x_k - x_old)/h;
	x_old = x_k;
}

void ImplicitEuler::ImplicitCalculateForces( MatrixXd& TVk, SparseMatrix<double>& forceGradient, VectorXd& x_k, VectorXd& f){
	// //gravity
	f.setZero();
	cout<<">>>>>>>>>>>>>>>>"<<endl;
	cout<<"TVK"<<endl;
	cout<<TVk-TV<<endl;
	for(unsigned int i=0; i<f.size()/3; i++){
		double vertex_mass = massVector(3*i+1);
		f(3*i+1) += vertex_mass*gravity;
	}
	cout<<"force after grav"<<endl;
	cout<<f<<endl;

	//elastic
	for(unsigned int i=0; i<M.tets.size(); i++){
		Vector4i indices = M.tets[i].verticesIndex;
		MatrixXd F_tet = M.tets[i].computeElasticForces(TVk, simTime%2);
		f.segment<3>(3*indices(0)) += F_tet.col(0);
		f.segment<3>(3*indices(1)) += F_tet.col(1);
		f.segment<3>(3*indices(2)) += F_tet.col(2);
		f.segment<3>(3*indices(3)) += F_tet.col(3);
	}
	// cout<<f<<endl<<endl;
	//damping
	f += rayleighCoeff*forceGradient*(x_k - x_old)/h;
	// cout<<f<<endl<<endl;
	cout<<">>>>>>>>>>>>>>>>"<<endl;
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

		
	}else if(solver.compare("lbfgs")==0){
		//renderLBFGS(ext_force);
		alglibLBFGS(ext_force);
		
	}else{
		cout<<"Solver not specified properly"<<endl;
		exit(0);
	}

	cout<<"*******************"<<endl;
	cout<< "New Pos"<<simTime<<endl;
	cout<<x_old<<endl<<endl;
	cout<< "New Vels"<<simTime<<endl;
	cout<<v_old<<endl;
	cout<<"*****************"<<endl<<endl;
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

