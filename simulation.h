#ifndef SIMULATION__H
#define SIMULATION__H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <pthread.h>
#include <math.h>



#include "Verlet.h"
#include "ImplicitEuler.h"
#include "ImplicitNewmark.h"


using namespace Eigen;
using namespace std;

class Simulation{

public:
	SolidMesh M;
	IntegratorAbstract* integrator;
	vector<int> mapV2TV;
	int iters;
	MatrixXd B;

	Simulation(void);
	int initializeSimulation(double deltaT, int iterations, char method, MatrixXi& TT, MatrixXd& TV, MatrixXd& B, vector<int>& moveVertices);
	
	int reIndexClampedVertices(vector<int>& moveVertices, MatrixXd& TV, MatrixXi& TT);
	void setInitPosition(vector<int> moveVertices, MatrixXd& TV, MatrixXi& TT);
	void setTVtoX(VectorXd &x, MatrixXd &TV);
	void calculateElasticForces(VectorXd &f, MatrixXd &TV);
	void calculateForceGradient(MatrixXd &TVk, SparseMatrix<double>& forceGradient);
	void headless();
	void render();

};

#endif