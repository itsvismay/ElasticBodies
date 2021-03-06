#ifndef implicit_newmark_h
#define implicit_newmark_h

#include "IntegratorsAbstract.h"

class ImplicitNewmark: public IntegratorAbstract{

public:
	//constants
	SparseMatrix<double> ZeroMatrix;
	SparseMatrix<double> Ident;

	SparseMatrix<double> forceGradient;
	SparseMatrix<double> grad_g;

	VectorXd x_k, v_k, f_old;
	MatrixXd TVk;
	double gamma = 0.5;
	double beta =0.25;
	
	void render(VectorXd& ext_force);
	void renderNewtonsMethod();
	void renderLBFGS();
	
	
	void initializeIntegrator(double ph, SolidMesh& pM, MatrixXd& pTV, MatrixXi& pTT);
	void NewmarkCalculateElasticForceGradient(MatrixXd& TVk, SparseMatrix<double>& forceGradient);
	void NewmarkCalculateForces( MatrixXd& TVk, SparseMatrix<double>& forceGradient, VectorXd& x_k, VectorXd& f);
	void NewmarkTVtoX(VectorXd& x_tv, MatrixXd& TVk);
	void NewmarkXtoTV(VectorXd& x_tv, MatrixXd& TVk);


};

#endif