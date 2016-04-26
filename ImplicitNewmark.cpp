#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <pthread.h>
#include <fstream>
#include <math.h>
#include <lbfgs.h>
#include "Eigen/SPQRSupport"
#include <Eigen/CholmodSupport>

#include "ImplicitNewmark.h"
#include "globals.h"

using namespace Eigen;
using namespace std;

typedef Eigen::Triplet<double> Trip;
typedef Matrix<double, 12, 1> Vector12d;

static lbfgsfloatval_t evaluate(void *impe, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step){
	ImplicitEuler* in = (ImplicitEuler*) impe;

	//from x to x_k
	for(int i=0; i< n; i++){
		in->x_k(i) = x[i];
	}
	

	in->ImplicitXtoTV(in->x_k, in->TVk);//TVk value changed in function
	in->ImplicitCalculateElasticForceGradient(in->TVk, in->forceGradient); 
	in->ImplicitCalculateForces(in->TVk, in->forceGradient, in->x_k, in->f);

	lbfgsfloatval_t fx = 0.0;
	for(int i=0; i<n; i++){
		fx+= 0.5*x[i]*in->massVector(i)*x[i] - in->massVector(i)*in->x_old(i)*x[i] - in->massVector(i)*in->h*in->v_old(i)*x[i]; //big G function, anti-deriv of g
		g[i] = in->massVector(i)*x[i] - in->massVector(i)*in->x_old(i) - in->massVector(i)*in->h*in->v_old(i) - in->h*in->h*in->f(i);
	}
	//force anti-derivative
	double strainE=0;
	for(unsigned int i=0; i<in->M.tets.size(); i++){
		strainE += in->M.tets[i].undeformedVol*in->M.tets[i].energyDensity;		
	}

	// ////////////////
	// double gnorm =0;
	// for(int i=0; i<n; i++){
	// 	gnorm+=g[i]*g[i];
	// }
	// // cout<<"Gnorm "<<sqrt(gnorm)<<endl;
	// // double xnorm =0;
	// // for(int i=0; i<n; i++){
	// // 	xnorm+=x[i]*x[i];
	// // }
	// // cout<<"xnorm "<<sqrt(xnorm)<<endl;

	// // for(int i=0; i<n; i++){
	// // 	xnorm+=x[i]*x[i];
	// // }
	// ////////////////

	fx+= in->h*in->h*(strainE);
	//damping anti-derivative
	fx += in->h*rayleighCoeff*((in->x_k.dot(in->f) - strainE) - in->f.dot(in->x_old));
	// cout<<"energy "<<fx<<endl;
	return fx;
}

static int progress(void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls){	

	printf("Iteration %d:\n", k);
    printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
    printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    printf("\n");
    return 0;	
}

void ImplicitNewmark::initializeIntegrator(double ph, SolidMesh& pM, MatrixXd& pTV){
	IntegratorAbstract::initializeIntegrator(ph, pM, pTV);
	ZeroMatrix.resize(3*vertsNum, 3*vertsNum);
	ZeroMatrix.setZero();
	Ident = MatrixXd::Identity(3*vertsNum, 3*vertsNum).sparseView();
	forceGradient.resize(3*vertsNum, 3*vertsNum);
	grad_g.resize(3*vertsNum, 3*vertsNum);
	x_k.resize(3*vertsNum);
	v_k.resize(3*vertsNum);
	x_k.setZero();
	v_k.setZero();
	TVk = TV;
}

void ImplicitNewmark::NewmarkTVtoX(VectorXd& x_tv, MatrixXd& TVk){
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

void ImplicitNewmark::NewmarkCalculateForces( MatrixXd& TVk, SparseMatrix<double>& forceGradient, VectorXd& x_k, VectorXd& f){
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