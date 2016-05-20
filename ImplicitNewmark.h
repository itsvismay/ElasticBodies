#ifndef IMPLICIT_NEWMARK__H
#define IMPLICIT_NEWMARK__H

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

class ImplicitNewmark: public IntegratorAbstract{

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
	
	
	void initializeIntegrator(double ph, SolidMesh& pM, MatrixXd& pTV, MatrixXi& pTT);
	void NewmarkCalculateElasticForceGradient(MatrixXd& TVk, SparseMatrix<double>& forceGradient);
	void NewmarkCalculateForces( MatrixXd& TVk, SparseMatrix<double>& forceGradient, VectorXd& x_k, VectorXd& f);
	void NewmarkTVtoX(VectorXd& x_tv, MatrixXd& TVk);
	void NewmarkXtoTV(VectorXd& x_tv, MatrixXd& TVk);


};

#endif