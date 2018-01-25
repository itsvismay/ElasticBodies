#ifndef implicit_euler_h
#define implicit_euler_h

#include "IntegratorsAbstract.h"
#include "../../alglib-cpp/src/optimization.h"


using namespace Eigen;
using namespace std;
using namespace alglib;

class ImplicitEuler: public IntegratorAbstract{

public:
	//constants
	SparseMatrix<double> ZeroMatrix;
	SparseMatrix<double> Ident;


	SparseMatrix<double> grad_g;

	double gamma = 0.5;
	double beta = 0.25;

	VectorXd x_k, v_k;
	MatrixXd TVk;
	int bfgsIterations = 0;
	int formulation = 1; //0 for old , 1 for new

	void render(VectorXd& ext_force);
	void renderNewtonsMethod(VectorXd& ext_force);
	void renderLBFGS(VectorXd& ext_force);
	int alglibLBFGSVismay(VectorXd& ext_force);
	void findgBlock(VectorXd& g_block, VectorXd& x, VectorXd& x_old, int ignorePast);
	void find_dEnergyBlock(VectorXd& g_block, VectorXd& y_k, int ignorePastIndex);
 	void find_d_dEnergyBlock(SparseMatrix<double>& grad_g_block, SparseMatrix<double>& forceGradientStaticBlock, SparseMatrix<double>& RegMassBlock);

	double ImplicitCalculateEnergyFunction();
	void initializeIntegrator(double ph, SolidMesh& pM, MatrixXd& pTV, MatrixXi& pTT);
	void ImplicitCalculateElasticForceGradient(MatrixXd& TVk, SparseMatrix<double>& forceGradient);
	void ImplicitCalculateForces( MatrixXd& TVk, SparseMatrix<double>& forceGradient, VectorXd& x_k, VectorXd& f);
	void ImplicitTVtoX(VectorXd& x_tv, MatrixXd& TVk);
	void ImplicitXtoTV(VectorXd& x_tv, MatrixXd& TVk);


};

#endif
