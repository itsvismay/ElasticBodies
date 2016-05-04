#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/barycenter.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>

#include <ctime>
#include "ConsistencyTests.h"
#include <fstream>
using namespace Eigen;
using namespace std;

#ifndef TUTORIAL_SHARED_PATH
#define TUTORIAL_SHARED_PATH "../shared"
#endif
ConsistencyTest::ConsistencyTest(void){}

void ConsistencyTest::runSpaceTests(Simulation& sim){
	double explicitTimestep = 1e-5;
	double implicitTimestep = 1e-3;
	cSim = sim;

	//Time stuff------
	time_t now = time(0);
	string dt = ctime(&now);//local time, replace all spaces and new lines
	dt.erase('\n');
	replace(dt.begin(), dt.end(), ' ', '-');
	//-----------------

	MatrixXd V;
	MatrixXi F;
	MatrixXd B;
	igl::readOBJ(TUTORIAL_SHARED_PATH "/beam.obj", V, F);

	int retT = -1;
	
	//SPACE TESTS Coarse
	retT= igl::copyleft::tetgen::tetrahedralize(V,F,"-pqa1500", cTV, cTT, cTF);
	igl::barycenter(cTV, cTT, cB);
	if(retT ==0){
		//Explicit
		test(explicitTimestep, 0.5, 'e', "../TestsResults/ConsistencyTests/"+dt+"/explicit/space/coarse"+to_string(explicitTimestep)+"/");
		//Implicit
		test(implicitTimestep, 0.5, 'i', "../TestsResults/ConsistencyTests/"+dt+"/implicit/space/coarse"+to_string(implicitTimestep)+"/");
	}

	// SPACE TESTS Middle
	retT = igl::copyleft::tetgen::tetrahedralize(V,F,"-pqa300", cTV, cTT, cTF);
	igl::barycenter(cTV, cTT, cB);
	if(retT ==0){
		//Explicit
		test(explicitTimestep, 0.5, 'e', "../TestsResults/ConsistencyTests/"+dt+"/explicit/space/middle:"+to_string(explicitTimestep)+"/");
		//Implicit
		test(implicitTimestep, 0.5, 'i', "../TestsResults/ConsistencyTests/"+dt+"/implicit/space/middle"+to_string(implicitTimestep)+"/");
	}
	//SPACE TESTS Fine
	retT = igl::copyleft::tetgen::tetrahedralize(V,F,"-pqa200", cTV, cTT, cTF);
	igl::barycenter(cTV, cTT, cB);
	if(retT ==0){
		//Explicit
		test(explicitTimestep, 0.5, 'e', "../TestsResults/ConsistencyTests/"+dt+"/explicit/space/fine:"+to_string(explicitTimestep)+"/");
		//Implicit
		test(implicitTimestep, 0.5, 'i', "../TestsResults/ConsistencyTests/"+dt+"/implicit/space/fine"+to_string(implicitTimestep)+"/");
	}

	return;
}

void ConsistencyTest::runTimeTests(Simulation& sim){
	cSim = sim;
	//Time stuff------
	time_t now = time(0);
	string dt = ctime(&now);//local time, replace all spaces and new lines
	dt.erase('\n');
	replace(dt.begin(), dt.end(), ' ', '-');
	//-----------------

	MatrixXd V;
	MatrixXi F;
	MatrixXd B;
	igl::readOBJ(TUTORIAL_SHARED_PATH "/beam.obj", V, F);
	igl::copyleft::tetgen::tetrahedralize(V,F,"-pq", cTV, cTT, cTF);
	igl::barycenter(cTV, cTT, cB);


	//EXPLICIT TIME TESTS
	double explicitTimestep = 1e-3;
	for(int i=0; i<3; i++){
		test(explicitTimestep, 0.5, 'e', "../TestsResults/ConsistencyTests/"+dt+"/explicit/time/timestep:"+to_string(explicitTimestep)+"/");
		explicitTimestep*=.1;
	}

	//IMPLICIT TIME TESTS
	double implicitTimestep = 1e-1;
	for(int i=0; i<3; i++){
		test(implicitTimestep, 0.5, 'i', "../TestsResults/ConsistencyTests/"+dt+"/implicit/time/timestep:"+to_string(implicitTimestep)+"/");
		implicitTimestep*=.1;
	}

	

	return;
}

void ConsistencyTest::test(double timestep, double printThisOften, char method, string printToHere){
	double seconds =0;
	int iters =0;

	cSim.initializeSimulation(timestep, 1, method, cTT, cTV, cB);
	cSim.integrator->fixVertices(0);
	cSim.integrator->fixVertices(1);
	double printForThisManySeconds = 2;
	while(seconds<printForThisManySeconds){
		iters+=1;
		cSim.render();
		if(iters*timestep>=printThisOften){
			seconds+=printThisOften;
			iters =0;
			printOBJ(seconds, printToHere);
		}
	}
	return;

}

void ConsistencyTest::printOBJ(double seconds, string printToHere){

	double refinement = 9;
	double t = ((refinement - 1)+1) / 9.0;


	VectorXd v = cB.col(2).array() - cB.col(2).minCoeff();
	v /= v.col(0).maxCoeff();

	vector<int> s;
	for (unsigned i=0; i<v.size();++i){
		if (v(i) < t){
			s.push_back(i);
		}
	}

	MatrixXd V_temp(s.size()*4,3);
	MatrixXi F_temp(s.size()*4,3);

	for (unsigned i=0; i<s.size();++i)
	{
		V_temp.row(i*4+0) = cSim.integrator->TV.row(cTT(s[i],0));
		V_temp.row(i*4+1) = cSim.integrator->TV.row(cTT(s[i],1));
		V_temp.row(i*4+2) = cSim.integrator->TV.row(cTT(s[i],2));
		V_temp.row(i*4+3) = cSim.integrator->TV.row(cTT(s[i],3));
		F_temp.row(i*4+0) << (i*4)+0, (i*4)+1, (i*4)+3;
		F_temp.row(i*4+1) << (i*4)+0, (i*4)+2, (i*4)+1;
		F_temp.row(i*4+2) << (i*4)+3, (i*4)+2, (i*4)+0;
		F_temp.row(i*4+3) << (i*4)+1, (i*4)+2, (i*4)+3;
	}

	cout<<printToHere + to_string(seconds)<<endl;
	system(("mkdir -p "+printToHere).c_str());
	igl::writeOBJ(printToHere + to_string(seconds)+".obj", V_temp, F_temp);

	return;
}