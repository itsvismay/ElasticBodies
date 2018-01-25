#ifndef implicit_newmark_h
#define implicit_newmark_h

#include "IntegratorsAbstract.h"
#include "../../alglib-cpp/src/optimization.h"

using namespace alglib;

class ImplicitNewmark: public IntegratorAbstract{

public:
	//constants
	SparseMatrix<double> ZeroMatrix;
	SparseMatrix<double> Ident;

	SparseMatrix<double> grad_g;

	double gamma = 0.5;
	double beta =0.25;

	VectorXd x_k, v_k, f_old;
	MatrixXd TVk;
	int bfgsIterations =0;
	int formulation = 1; // 0 for old, 1 for new

	void render(VectorXd& ext_force);
	void renderNewtonsMethod(VectorXd& ext_force);
	void renderLBFGS();
	int alglibBFGS(VectorXd &ext_force);

	void findgBlock(VectorXd& g_block, VectorXd& x, VectorXd& x_old, int ignorePast, double gamma, double beta);
	void find_dEnergyBlock(VectorXd& g_block, VectorXd& y_k, int ignorePastIndex);
	void find_d_dEnergyBlock(SparseMatrix<double>& grad_g_block, SparseMatrix<double>& forceGradientStaticBlock, SparseMatrix<double>& RegMassBlock);
	double find_Energy();

	void initializeIntegrator(double ph, SolidMesh& pM, MatrixXd& pTV, MatrixXi& pTT);
	void NewmarkCalculateElasticForceGradient(MatrixXd& TVk, SparseMatrix<double>& forceGradient);
	void NewmarkCalculateForces( MatrixXd& TVk, SparseMatrix<double>& forceGradient, VectorXd& x_k, VectorXd& f);
	void NewmarkTVtoX(VectorXd& x_tv, MatrixXd& TVk);
	void NewmarkXtoTV(VectorXd& x_tv, MatrixXd& TVk);


};

#endif
