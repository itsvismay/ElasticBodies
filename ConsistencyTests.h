#ifndef CONSISTENCY__H
#define CONSISTENCY__H

#include <igl/writeOBJ.h>

#include "simulation.h"
#include "globals.h"


using namespace Eigen;
using namespace std;

class ConsistencyTest{

public:
	double h;
	int cmethod;
	MatrixXi cTT;
	MatrixXd cTV;
	Simulation cSim;
	vector<int> cMapV2TV;

	ConsistencyTest(void);
	void initializeTest(double deltaT, char method, MatrixXi& TT, MatrixXd& TV, Simulation& sim, vector<int> mapV2TV);
	void runTests();


};

#endif