#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/barycenter.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>

#include <sys/wait.h>
#include <ctime>
#include "ConsistencyTests.h"
#include <fstream>
using namespace Eigen;
using namespace std; 

ConsistencyTest::ConsistencyTest(void){
	printThisOften = 0.1;
	printForThisManySeconds = 5;
	// spaceIterations = 1;
	timeIterations = 1;
}

bool ConsistencyTest::checkAllAccuracy(){
	MatrixXd TV;
	MatrixXi TT;
	MatrixXd B;

	TT.resize(2, 4);
	TT<<  0, 1, 2, 3,
			4, 0, 2, 3;

	TV.resize(5, 3);
	TV << 10, 0, 0, //affect
			0, 10, 0,
			0, 0, 10,
			0, 0, 0,
			0, -10, 0;

	int timestep =0.01;
	// return 	checkVerletAccuracy(timestep, 1, 'e', TT, TV, B)
	// 		&&
			// checkNewmarkAccuracy(timestep, 1, 'n', TT, TV, B)
			// &&
		return	checkEulerAccuracy(timestep, 1, 'i', TT, TV, B);
}

bool ConsistencyTest::checkVerletAccuracy(double timestep, int iterations, char method, MatrixXi TT, MatrixXd TV, MatrixXd B ){
// 	vector<int> moveVertices;
// 	vector<int> fixedVertices;
// 	Simulation vSim;
// 	vSim.initializeSimulation(timestep, iterations, method, TT, TV, B, moveVertices, fixedVertices);
// 	vSim.render();
// 	cout<<"Verlet"<<endl<<vSim.integrator->TV<<endl;
// 	return 1;
}
bool ConsistencyTest::checkEulerAccuracy(double timestep, int iterations, char method, MatrixXi TT, MatrixXd TV, MatrixXd B ){
	vector<int> moveVertices;
	vector<int> fixedVertices;
	gravity = 10;
	Simulation vSim;
	cout<<TV<<endl;
	cout<<TT<<endl;
	vSim.initializeSimulation(timestep, iterations, method, TT, TV, B, moveVertices, fixedVertices, 2e6, 0.35);
	vSim.render();
	cout<<"Euler"<<endl<<vSim.integrator->TV<<endl;
	return 1;
}
bool ConsistencyTest::checkNewmarkAccuracy(double timestep, int iterations, char method, MatrixXi TT, MatrixXd TV, MatrixXd B ){
// 	vector<int> moveVertices;
// 	vector<int> fixedVertices;
// 	Simulation vSim;
// 	vSim.initializeSimulation(timestep, iterations, method, TT, TV, B, moveVertices, fixedVertices);
// 	vSim.render();
// 	cout<<"Newmark"<<endl<<vSim.integrator->TV<<endl;
// 	return 1;
}

void ConsistencyTest::runAllTests(){
	checkAllAccuracy();
	exit(0);

 	//Time stuff------
 	time_t now = time(0);
 	string dt = ctime(&now);//local time, replace all spaces and new lines
 	dt.erase('\n');
 	replace(dt.begin(), dt.end(), ' ', '-');
 	//-----------------
 	int spaceStep = 84355;

 	MatrixXd V;
 	MatrixXi F;
 	MatrixXd B;
 	MatrixXd TV;
 	MatrixXi TT;
 	MatrixXi TF;
 	igl::readOBJ(TUTORIAL_SHARED_PATH "shared/springTruncd.obj", V, F);
	
// 	pid_t pids[spaceIterations];

	igl::copyleft::tetgen::tetrahedralize(V,F,"-pRq10", TV, TT, TF);
	igl::barycenter(TV, TT, B);

	runTestRow(spaceStep, TV, B, TT, dt);

// 	// for(int i=0; i<spaceIterations; i++){
// 	// 	if((pids[i] = fork())<0){
// 	// 		cout<< "Fork error"<<endl;
// 	// 		abort();
// 	// 	}else if(pids[i] == 0){
// 	// 		igl::copyleft::tetgen::tetrahedralize(V,F,("-pqa"+to_string(spaceStep)).c_str(), TV, TT, TF);
// 	// 		igl::barycenter(TV, TT, B);
// 	// 		runTestRow(spaceStep, TV, B, TT, dt);
// 	// 		exit(0);
// 	// 	}
// 	// 		spaceStep*=.5;
// 	// }

// 	// //wait for children to exit
// 	// int status;
// 	// pid_t pid;
// 	// while(spaceIterations>0){
// 	// 	pid = wait(&status);
// 	// 	cout<<"Exited --"<<pid<<" "<<status<<endl;
// 	// 	spaceIterations--;
// 	// }
// 	cout<<"------------------------"<<endl;
// 	cout<<"ALL EXITED"<<endl;
}

void ConsistencyTest::runTestRow(int spaceStep, MatrixXd& TV, MatrixXd& B, MatrixXi& TT, string dt){
	// runVerletTestRow(spaceStep, TV, B, TT, dt);
	runImpEulerTestRow(spaceStep, TV, B , TT, dt);
	// runNewmarkTestRow(spaceStep, TV, B, TT, dt);
}

// void ConsistencyTest::runVerletTestRow(int spaceStep, MatrixXd& TV, MatrixXd& B, MatrixXi& TT, string dt){
// 	double verletTimestep = 1e-5;
// 	for(int i=0; i<timeIterations; i++){
// 		cout<<"verlet" +to_string(spaceStep)+"time"+to_string(verletTimestep)<<endl;
// 		// test(verletTimestep, 'e', CONSISTENCY_TEST_SAVE_PATH"TestsResults/ConsistencyTests/"+dt+"/explicit/space"+to_string(spaceStep)+"/timestep:"+to_string(verletTimestep)+"/", TT, TV, B);
// 		verletTimestep*=.1;
// 	}
// }

void ConsistencyTest::runImpEulerTestRow(int spaceStep, MatrixXd& TV, MatrixXd& B, MatrixXi& TT, string dt){
	double implicitTimestep = 1e-1;
	for(int i=0; i<timeIterations; i++){
		cout<<"euler" +to_string(spaceStep)+"time"+to_string(implicitTimestep)<<endl;
		test(implicitTimestep, 'i', CONSISTENCY_TEST_SAVE_PATH"TestsResults/ConsistencyTests/"+dt+"/implicit/svk/space"+to_string(spaceStep)+"/timestep:"+to_string(implicitTimestep)+"/", TT, TV, B);
		implicitTimestep*=.1;
	}
}

// void ConsistencyTest::runNewmarkTestRow(int spaceStep, MatrixXd& TV, MatrixXd& B, MatrixXi& TT, string dt){
// 	double newmarkTimestep = 1e-1;
// 	for(int i=0; i<timeIterations; i++){
// 		cout<<"newmark" +to_string(spaceStep)+"time"+to_string(newmarkTimestep)<<endl;
// 		// test(newmarkTimestep, 'n', CONSISTENCY_TEST_SAVE_PATH"TestsResults/ConsistencyTests/"+dt+"/newmark/space"+to_string(spaceStep)+"/timestep:"+to_string(newmarkTimestep)+"/", TT, TV, B);
// 		newmarkTimestep*=.1;
// 	}
// }

void ConsistencyTest::test(double timestep, char method, string printToHere, MatrixXi cTT, MatrixXd cTV, MatrixXd cB){
	double seconds =0;
	int iters =0;
	int numberOfPrints =0;
	vector<int> moveVertices;
	vector<int> fixedVertices;

	Simulation cSim;

	//fix vertices
	for(int i=0; i<cTV.rows(); i++){
		if(cTV.row(i)[1]<=41 && cTV.row(i)[1]>39){
			fixedVertices.push_back(i);
		}
	}
	cSim.initializeSimulation(timestep, 1, method, cTT, cTV, cB, moveVertices, fixedVertices, 1.2e6, 0.35);

	while(seconds<printForThisManySeconds){
		iters+=1;
		cSim.render();
		if(iters*timestep>=printThisOften){
			seconds+=printThisOften;
			iters =0;
			printOBJ(numberOfPrints, printToHere, cB, cSim, cTT);
			numberOfPrints+=1;
		}
	}
	return;
}

void ConsistencyTest::printOBJ(int numberOfPrints, string printToHere, MatrixXd& cB, Simulation& cSim, MatrixXi& cTT){

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
		V_temp.row(i*4+0) = cSim.integrator->TV.row(cSim.integrator->TT(s[i],0));
		V_temp.row(i*4+1) = cSim.integrator->TV.row(cSim.integrator->TT(s[i],1));
		V_temp.row(i*4+2) = cSim.integrator->TV.row(cSim.integrator->TT(s[i],2));
		V_temp.row(i*4+3) = cSim.integrator->TV.row(cSim.integrator->TT(s[i],3));
		F_temp.row(i*4+0) << (i*4)+0, (i*4)+1, (i*4)+3;
		F_temp.row(i*4+1) << (i*4)+0, (i*4)+2, (i*4)+1;
		F_temp.row(i*4+2) << (i*4)+3, (i*4)+2, (i*4)+0;
		F_temp.row(i*4+3) << (i*4)+1, (i*4)+2, (i*4)+3;
	}

	cout<<printToHere + to_string(numberOfPrints)<<endl;
	system(("mkdir -p "+printToHere).c_str());
	igl::writeOBJ(printToHere + to_string(numberOfPrints)+".obj", V_temp, F_temp);

	return;
}
