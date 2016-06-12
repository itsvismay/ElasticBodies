#ifndef CONSISTENCY__H
#define CONSISTENCY__H

#include <igl/writeOBJ.h>

#include "simulation.h"
#include "globals.h"


using namespace Eigen;
using namespace std;

#define TUTORIAL_SHARED_PATH "../"
#define CONSISTENCY_TEST_SAVE_PATH  "../"

class ConsistencyTest{

public:
	MatrixXi cTT;
	MatrixXi cTF;
	MatrixXd cTV;
	MatrixXd cB;
	
	double h;
	int cmethod;
	Simulation cSim;
	vector<int> cMapV2TV;

	double printThisOften, printForThisManySeconds;

	ConsistencyTest(void);
	void runTimeTests(Simulation& sim);
	void runSpaceTests(Simulation& sim);
	
	void test(double timestep, char method, string printToHere);
	void printOBJ(int numberOfPrints, string printToHere);

};

#endif
