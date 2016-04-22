#ifndef IMPLICITEULER__H
#define IMPLICITEULER__H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <pthread.h>
#include <fstream>
#include <math.h>

#include "IntegratorsAbstract.h"

using namespace Eigen;
using namespace std;

class ImplicitEuler: public IntegratorAbstract{

public:
	//constants
	SparseMatrix<double> ZeroMatrix;
	SparseMatrix<double> Ident;

	SparseMatrix<double> forceGradient;
	SparseMatrix<double> grad_g;

	VectorXd x_k, v_k;
	MatrixXd TVk;

	void render();
	void renderNewtonsMethod();
	void renderLBFGS();
	
	void initializeIntegrator(double ph, SolidMesh& pM, MatrixXd& pTV);
	void ImplicitCalculateElasticForceGradient(MatrixXd& TVk, SparseMatrix<double>& forceGradient);
	void ImplicitCalculateForces( MatrixXd& TVk, SparseMatrix<double>& forceGradient, VectorXd& x_k, VectorXd& f);
	void ImplicitTVtoX(VectorXd& x_tv, MatrixXd& TVk);
	void ImplicitXtoTV(VectorXd& x_tv, MatrixXd& TVk);


};

#endif