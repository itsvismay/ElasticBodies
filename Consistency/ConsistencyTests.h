#ifndef consistency_tests_h
#define consistency_tests_h

#include "simulation.h"
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/writeOBJ.h>
#include <igl/barycenter.h>
#include <igl/readOFF.h>
#include <igl/readOBJ.h>

class ConsistencyTest{

public:

	double printThisOften, printForThisManySeconds;
	int spaceIterations, timeIterations;

	ConsistencyTest(void);

	void replaceWithMain();
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
