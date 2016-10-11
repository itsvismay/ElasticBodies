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

#define HOME_SAVED_PATH "/u/vismay/ElasticBodies/"
#define OUTPUT_SAVED_PATH "/scratch/cluster/vismay/"

// #define HOME_SAVED_PATH "/home/vismay/ElasticBodies/"
// #define OUTPUT_SAVED_PATH "/home/vismay/ElasticBodies/"

class Simulation{

public:
	SolidMesh M;
	IntegratorAbstract* integrator;
	vector<int> mapV2TV;
	int iters;
	MatrixXd sB;

	//TODO: cleanup this code (its so terrible)
	MatrixXd thisTV;
	MatrixXd thisx;

	Simulation(void);
	int initializeSimulation(double deltaT, int iterations, char method, MatrixXi& TT, MatrixXd& TV, MatrixXd& B, vector<int>& moveVertices, vector<int> fixVertices, double youngs, double poissons);
	
	void binarySearchYoungs(vector<int> moveVertices, MatrixXd& TV, MatrixXi& TT, int fv, MatrixXd& B);
	void staticSolveStepNewtonsMethod(double move_step, int ignorePastIndex, vector<int>& moveVertices, MatrixXd& TV,  MatrixXi& TT);
	void syntheticTests(vector<int> moveVertices, MatrixXd& TV, MatrixXi& TT, int fv, MatrixXd& B);
	void reIndexTVandTT(vector<int> newVertsIndices, int sizeFixed, int sizeMove,MatrixXd& TV, MatrixXi& TT, MatrixXd& newTV, MatrixXi& newTT);
	
	void staticSolveStepLBFGS(double move_step, int ignorePastIndex, vector<int>& moveVertices, MatrixXd& TV,  MatrixXi& TT);

	void setInitPosition(vector<int> moveVertices, MatrixXd& TV, MatrixXi& TT, int fv, MatrixXd& B);
	void printObj(int numberOfPrints, MatrixXd& TV, MatrixXi& TT, MatrixXd& B);
	void setTVtoX(VectorXd &x, MatrixXd &TV);
	void calculateElasticForces(VectorXd &f, MatrixXd &TV);
	void calculateForceGradient(MatrixXd &TVk, SparseMatrix<double>& forceGradient);
	void headless();
	void render();

};

#endif
