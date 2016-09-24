#ifndef CONSISTENCY__H
#define CONSISTENCY__H

#include <igl/writeOBJ.h>

#include "simulation.h"
#include "globals.h"


using namespace Eigen;
using namespace std;

#define TUTORIAL_SHARED_PATH "/u/vismay/ElasticBodies"
#define CONSISTENCY_TEST_SAVE_PATH "/u/vismay/ElasticBodies" 

class ConsistencyTest{

public:


	double printThisOften, printForThisManySeconds;
	int spaceIterations, timeIterations;

	ConsistencyTest(void);

	void runVerletTestRow(int spaceStep, MatrixXd& TV, MatrixXd& B, MatrixXi& TT, string dt);
	void runImpEulerTestRow(int spaceStep, MatrixXd& TV, MatrixXd& B, MatrixXi& TT, string dt);
	void runNewmarkTestRow(int spaceStep, MatrixXd& TV, MatrixXd& B, MatrixXi& TT, string dt);
	void runTestRow(int spaceStep, MatrixXd& TV, MatrixXd& B, MatrixXi& TT, string dt);
	void runAllTests();
	void test(double timestep, char method, string printToHere, MatrixXi cTT, MatrixXd cTV, MatrixXd cB );
	void printOBJ(int numberOfPrints, string printToHere, MatrixXd& cB, Simulation& cSim, MatrixXi& cTT);

	bool checkAllAccuracy();
	bool checkVerletAccuracy(double timestep, int iterations, char method, MatrixXi TT, MatrixXd TV, MatrixXd B );
	bool checkEulerAccuracy(double timestep, int iterations, char method, MatrixXi TT, MatrixXd TV, MatrixXd B );
	bool checkNewmarkAccuracy(double timestep, int iterations, char method, MatrixXi TT, MatrixXd TV, MatrixXd B );
};

#endif
