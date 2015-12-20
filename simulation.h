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
	SolidMesh M;
	int vertices;
	VectorXd x_old, v_old, vertex_masses;
	VectorXd f;
	vector<int> mapV2TV;
	MatrixXd TV;
	SparseMatrix<double> ZeroMatrix;


	Simulation(void);
	void renderExplicit();
	void renderImplicit();
	void createInvMassMatrix();
	void createXFromTet();
	void createForceVector();
	void createTVFromTet();
	void initializeSimulation(double deltaT, MatrixXi& TT_One, MatrixXd& TV_One, vector<int>& map);
	
	void setXIntoTV(VectorXd& x_new);

	void calculateGravity();
	void calculateElasticForce();

	void fixVertices(int fixed);

	SparseMatrix<double> ImplicitCalculateElasticForceGradient(MatrixXd& TVk);
	VectorXd ImplicitCalculateForces(MatrixXd& TVk, SparseMatrix<double>& forceGradient, VectorXd& v_k);
	MatrixXd ImplicitXtoTV(VectorXd& x_tv);
	void ImplicitXfromTV(VectorXd& x_n1, MatrixXd& TVk);
};