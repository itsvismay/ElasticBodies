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
	VectorXd x_old, v_old, vertex_masses;
	VectorXd f;
	vector<int> mapV2TV;
	MatrixXd TV;
	SparseMatrix<double> ZeroMatrix;
	char integrator;

	//used in render implicit
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

	SparseMatrix<double> ImplicitCalculateElasticForceGradient(MatrixXd& TVk);
	VectorXd ImplicitCalculateForces(MatrixXd& TVk, SparseMatrix<double>& forceGradient, VectorXd& v_k);
	MatrixXd ImplicitXtoTV(VectorXd& x_tv);
	void ImplicitXfromTV(VectorXd& x_n1, MatrixXd& TVk);
};