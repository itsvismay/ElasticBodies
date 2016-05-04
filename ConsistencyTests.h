#ifndef CONSISTENCY__H
#define CONSISTENCY__H

#include <igl/writeOBJ.h>

#include "simulation.h"
#include "globals.h"


using namespace Eigen;
using namespace std;

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

	ConsistencyTest(void);
	void runTimeTests(Simulation& sim);
	void runSpaceTests(Simulation& sim);
	
	void timeTest(double timestep, double printThisOften, char method, string printToHere);
	void printOBJ(double number, string printToHere);

};

#endif