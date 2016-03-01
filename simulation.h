#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <pthread.h>
#include <fstream>
#include <math.h>


#include "solidmesh.h"

using namespace Eigen;
using namespace std;
class Simulation{

public:
	int t=0;
	double timestep;
	SparseMatrix<double> InvMass;
	vector<int> fixedVertices;
	SolidMesh M;
	int vertices;
	VectorXd x_old, v_old, vertex_masses, x_k, v_k;
	VectorXd f;
	vector<int> mapV2TV;
	MatrixXd TV, TVk;
	char integrator;

	//used in render implicit
	SparseMatrix<double> ZeroMatrix;
	SparseMatrix<double> global_gradForces;
	SparseMatrix<double> grad_g;
	SparseMatrix<double> Ident;
	SparseMatrix<double> forceGradient;


	Simulation(void);
	void render();
	void renderExplicit();
	void renderImplicit();
	void createInvMassMatrix();
	void createXFromTet();
	void createForceVector();
	void createTVFromTet();
	void initializeSimulation(double deltaT, char method, MatrixXi& TT_One, MatrixXd& TV_One, vector<int>& map);
	
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

	// void insertToSpringSet(int i1, int i2);

	
};