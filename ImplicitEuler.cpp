#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <pthread.h>
#include <fstream>
#include <math.h>
#include "Eigen/SPQRSupport"
#include <Eigen/CholmodSupport>

#include "ImplicitEuler.h"
#include "globals.h"

using namespace Eigen;
using namespace std;

typedef Eigen::Triplet<double> Trip;
typedef Matrix<double, 12, 1> Vector12d;

void ImplicitEuler::initializeIntegrator(int ph, SolidMesh& pM, MatrixXd& pTV){
	cout<<"here again"<<endl;
	IntegratorAbstract::initializeIntegrator(ph, pM, pTV);

	ZeroMatrix.resize(3*vertsNum, 3*vertsNum);
	ZeroMatrix.setZero();
	Ident = MatrixXd::Identity(3*vertsNum, 3*vertsNum).sparseView();
	forceGradient.resize(3*vertsNum, 3*vertsNum);
	grad_g.resize(3*vertsNum, 3*vertsNum);
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

void ImplicitEuler::ImplicitCalculateForces( MatrixXd& TVk, SparseMatrix<double>& forceGradient, VectorXd& x_k, VectorXd& f){
	// //gravity
	f.setZero();
	for(unsigned int i=0; i<M.tets.size(); i++){
		double vertex_mass = M.tets[i].undeformedVol/4;//assume const density 1
		Vector4i indices = M.tets[i].verticesIndex;
		f(3*indices(0)+1) += vertex_mass*gravity;
		f(3*indices(1)+1) += vertex_mass*gravity; 
		f(3*indices(2)+1) += vertex_mass*gravity;
		f(3*indices(3)+1) += vertex_mass*gravity;
	}

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
	return;
}

void ImplicitEuler::ImplicitCalculateElasticForceGradient(MatrixXd& TVk, SparseMatrix<double>& forceGradient){
	forceGradient.setZero();
	
	vector<Trip> triplets1;
	triplets1.reserve(3*vertsNum*3*vertsNum);	
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

void ImplicitEuler::render(){
	simTime+=1;
	cout<<"i"<<simTime<<endl;

	//Implicit Code
	// cout<<"diff in x"<<endl;
	// cout<<x_k - x_old<<endl<<endl;
	v_k.setZero();
	x_k.setZero();
	x_k = x_old;
	v_k = v_old;

	forceGradient.setZero();
	bool Nan=false;
	int NEWTON_MAX = 100, i =0;
	cout<<"--------"<<"-------"<<endl;
	cout<<"x_k"<<endl;
	cout<<x_k<<endl<<endl;
	cout<<"v_k"<<endl;
	cout<<v_k<<endl<<endl;
	cout<<"--------------------"<<endl;
	for( i=0; i<NEWTON_MAX; i++){
		grad_g.setZero();
	
		ImplicitXtoTV(x_k, TVk);//TVk value changed in function
		ImplicitCalculateElasticForceGradient(TVk, forceGradient); 
		ImplicitCalculateForces(TVk, forceGradient, x_k, f);

		// VectorXd g = x_k - x_old -h*v_old -h*h*InvMass*f;
		// grad_g = Ident - h*h*InvMass*forceGradient - h*rayleighCoeff*InvMass*forceGradient;
		
		VectorXd g = RegMass*x_k - RegMass*x_old - h*RegMass*v_old - h*h*f;
		grad_g = RegMass - h*h*forceGradient - h*rayleighCoeff*forceGradient;
	
		// cout<<"G"<<t<<endl;
		// cout<<g<<endl<<endl;
		// cout<<"G Gradient"<<t<<endl;
		// cout<<grad_g<<endl;

		//solve for delta x
		// Conj Grad
		// ConjugateGradient<SparseMatrix<double>> cg;
		// cg.compute(grad_g);
		// VectorXd deltaX = -1*cg.solve(g);

		// Sparse Cholesky LL^T
		// SimplicialLLT<SparseMatrix<double>> llt;
		// llt.compute(grad_g);
		// VectorXd deltaX = -1* llt.solve(g);

		//Sparse QR 
		SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> sqr;
		sqr.compute(grad_g);
		VectorXd deltaX = -1*sqr.solve(g);

		// CholmodSimplicialLLT<SparseMatrix<double>> cholmodllt;
		// cholmodllt.compute(grad_g);
		// VectorXd deltaX = -cholmodllt.solve(g);
		

		x_k+=deltaX;
		if(x_k != x_k){
			Nan = true;
			break;
		}
		if(g.squaredNorm()<.00000001){
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
	v_old.setZero();
	v_old = (x_k - x_old)/h;
	x_old = x_k;
	cout<<"*******************"<<endl;
	cout<< "New Pos"<<simTime<<endl;
	cout<<x_old<<endl<<endl;
	cout<< "New Vels"<<simTime<<endl;
	cout<<v_old<<endl;
	cout<<"*****************"<<endl<<endl;
	ImplicitXtoTV(x_old, TV);

	IntegratorAbstract::printInfo();
}
// void Simulation::renderImplicit(){


// 	//LBFGS
// 	t+=1;
// 	int N=3*vertices;
// 	int i, ret = 0;
//     lbfgsfloatval_t fx;
//     lbfgsfloatval_t *x = lbfgs_malloc(N);
//     lbfgs_parameter_t param;
//     if (x == NULL) {
//         printf("ERROR: Failed to allocate a memory block for variables.\n");
//     }

//     /* Initialize the variables. */
//     x_k.setZero();
//     for (i = 0;i < N; i++) {
       
//        x[i] = x_old(i);
//     }
//     x_k = x_old;
//     v_k = v_old;
//     cout<<"--------"<<t<<"-------"<<endl;
// 	cout<<"x_old"<<endl;
// 	cout<<x_old<<endl<<endl;
// 	cout<<"v_old"<<endl;
// 	cout<<v_old<<endl<<endl;
// 	cout<<"--------------------"<<endl;

//     /* Initialize the parameters for the L-BFGS optimization. */
//     lbfgs_parameter_init(&param);
//     //param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;
//     // param.gtol = 0.0001;
//     // param.ftol = 0.000001;
//     param.epsilon = 0.00001;
//     /*
//         Start the L-BFGS optimization; this will invoke the callback functions
//         evaluate() and progress() when necessary.
//      */
//     ret = lbfgs(N, x, &fx, evaluate, progress, this, &param);
    
//     /* Report the result. */
//     cout<<"Ret X val"<<endl;
//     for(i=0; i<N; i++){
//     	printf("%f\n", x[i]);;
//     	x_k(i) = x[i];
//     }
//     v_old = (x_k - x_old)/h;
//     x_old = x_k;
//     ImplicitXtoTV(x_old, TV);
    

//     cout<<"*******************"<<endl;
// 	cout<< "New Pos"<<t<<endl;
// 	cout<<x_old<<endl<<endl;
// 	cout<< "New Vels"<<t<<endl;
// 	cout<<v_old<<endl;
// 	cout<<"*****************"<<endl<<endl;

//     lbfgs_free(x);
// }

