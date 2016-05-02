#include "ConsistencyTests.h"

using namespace Eigen;
using namespace std;


ConsistencyTest::ConsistencyTest(void){}

void ConsistencyTest::initializeTest(double deltaT, char method, MatrixXi& TT, MatrixXd& TV, Simulation& sim, vector<int> mapV2TV){
	h = deltaT;
	cmethod = method;
	cTT = TT;
	cTV = TV;
	cSim = sim;
	cMapV2TV = mapV2TV;
	runTests();
}

void ConsistencyTest::runTests(){
	cout<<"Begin Tests"<<endl;
	int iters = 2;
	MatrixXd V;
	V.resize(cMapV2TV.size(), 3);

	cSim.initializeSimulation(h, iters, cmethod, cTT, cTV);

	// //SET V from TV
	// for(unsigned int i=0; i<cSim.mapV2TV.size(); i++){
	// 	V.row(i) = cSim.integrator->TV.row(cMapV2TV[i]);
	// 	cout<<cMapV2TV[i]<<endl;
	// }

	return;
}