#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <pthread.h>
#include <fstream>
#include <math.h>
#include <lbfgs.h>

#include "Verlet.h"


using namespace Eigen;
using namespace std;

class Simulation{

public:
	// int t=0;
	// double timestep;
	// SparseMatrix<double> InvMass;
	// SparseMatrix<double> RegMass;
	// vector<int> fixedVertices;
	SolidMesh M;
	// int vertices;
	// VectorXd x_old, v_old, vertex_masses, x_k, v_k, f;
	// vector<int> mapV2TV;
	// MatrixXd TV, TVk;
	IntegratorAbstract* integrator;

	// //used in render implicit
	// SparseMatrix<double> ZeroMatrix;
	// SparseMatrix<double> grad_g;
	// SparseMatrix<double> Ident;
	// SparseMatrix<double> forceGradient;


	Simulation(void);
	void render();
	void renderExplicit();
	void renderImplicit();
	void renderNewmark();
	void createMassMatrix();
	void createXFromTet();
	void createForceVector();
	void createTVFromTet();
	void initializeSimulation(double deltaT, 
		char method, 
		MatrixXi& TT, 
		MatrixXd& TV, 
		vector<int>& map);
	
	void setXIntoTV(VectorXd& x_new);

	void calculateGravity();
	void calculateElasticForce();

	void fixVertices(int fixed);
	bool isFixed(int vert);

	void ImplicitCalculateElasticForceGradient(MatrixXd& TVk, SparseMatrix<double>& forceGradient);
	void ImplicitCalculateForces(MatrixXd& TVk, SparseMatrix<double>& forceGradient, VectorXd& v_k, VectorXd& f);
	void ImplicitXtoTV(VectorXd& x_tv, MatrixXd& TVk);
	void ImplicitXfromTV(VectorXd& x_n1, MatrixXd& TVk);
	void ImplicitTVtoX(VectorXd& x_tv, MatrixXd& TVk);

	// static lbfgsfloatval_t evaluate(void *instance, 
	// 	const lbfgsfloatval_t *x, 
	// 	lbfgsfloatval_t *g, 
	// 	const int n, 
	// 	const lbfgsfloatval_t step);
	// static int progress(void *instance,
	//     const lbfgsfloatval_t *x,
	//     const lbfgsfloatval_t *g,
	//     const lbfgsfloatval_t fx,
	//     const lbfgsfloatval_t xnorm,
	//     const lbfgsfloatval_t gnorm,
	//     const lbfgsfloatval_t step,
	//     int n,
	//     int k,
	//     int ls);	
};