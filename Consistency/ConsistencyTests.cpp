#include "ConsistencyTests.h"

ConsistencyTest::ConsistencyTest(void){
	
}

// bool ConsistencyTest::checkAllAccuracy(){
// 	MatrixXd TV;
// 	MatrixXi TT;
// 	MatrixXd B;

// 	TT.resize(2, 4);
// 	TT<<  0, 1, 2, 3,
// 			4, 0, 2, 3;

// 	TV.resize(5, 3);
// 	TV << 10, 10, 0, //affect
// 			0, 20, 0,
// 			0, 10, 10,
// 			0, 10, 0,
// 			0, 0, 0;

// 	// return 	checkVerletAccuracy(timestep, 1, 'e', TT, TV, B)
// 	// 		&&
// 	return checkNewmarkAccuracy(timestep, 1, 'n', TT, TV, B);
// 			// &&
// 		//return	checkEulerAccuracy(timestep, 1, 'i', TT, TV, B);
// }

bool ConsistencyTest::checkVerletAccuracy(double timestep, int iterations, char method, MatrixXi TT, MatrixXd TV, MatrixXd B ){
// 	vector<int> moveVertices;
// 	vector<int> fixedVertices;
// 	Simulation vSim;
// 	vSim.initializeSimulation(timestep, iterations, method, TT, TV, B, moveVertices, fixedVertices);
// 	vSim.render();
// 	cout<<"Verlet"<<endl<<vSim.integrator->TV<<endl;
// 	return 1;
}

// bool ConsistencyTest::checkEulerAccuracy(double timestep, int iterations, char method, MatrixXi TT, MatrixXd TV, MatrixXd B ){
// 	vector<int> moveVertices;
// 	vector<int> fixedVertices;
// 	gravity = -0;
// 	rayleighCoeff =0;
// 	Simulation vSim;

// 	cout<<"timestep \n"<<timestep<<endl;
// 	vSim.initializeSimulation(timestep, iterations, method, TT, TV, B, moveVertices, fixedVertices, 2e6, 0.35);
// 	vSim.integrator->v_old(0)+=1;
// 	for(int i=0; i<10; i++)
// 		vSim.render();

// 	VectorXd realXfromMathematica(15);
// 	realXfromMathematica<< 10.0292,
// 					          10,
// 					    0.012497,
// 					   0.0291688,
// 					          20,
// 					 -0.00416129,
// 					   0.0125097,
// 					          10,
// 					     9.99583,
// 					   0.0291678,
// 					          10,
// 					 -0.00416102,
// 					   0.0291688,
// 					-7.95354e-08,
// 					 -0.00416129;

// 	bool q = (vSim.integrator->x_old - realXfromMathematica).norm()<0.0001;
// 	return q;
// }

// bool ConsistencyTest::checkNewmarkAccuracy(double timestep, int iterations, char method, MatrixXi TT, MatrixXd TV, MatrixXd B ){
// 	vector<int> moveVertices;
// 	vector<int> fixedVertices;
// 	gravity = -0;
// 	rayleighCoeff =0;
// 	Simulation vSim;

// 	vSim.integrator->v_old(0)+=1;
// 	vSim.initializeSimulation(timestep, iterations, method, TT, TV, B, moveVertices, fixedVertices, 2e6, 0.35);
// 	vSim.render();

// 	cout<<"Newmark"<<endl<<vSim.integrator->x_old<<endl;
// 	return 1;
// }

void ConsistencyTest::runAllTests(){
	// if(checkAllAccuracy()){
	// 	cout<<"Accuracy Tests Passed"<<endl;
	// }
	// exit(0);
}

void ConsistencyTest::test(double timestep, char method, string printToHere, MatrixXi cTT, MatrixXd cTV, MatrixXd cB){
	system(("mkdir -p "+printToHere).c_str());
	system(("(git log -1; echo ''; echo 'Ran Test On:'; date;) >>"+printToHere+"log.txt").c_str());	
	
	double seconds =0;
	int iters =0;
	int numberOfPrints =0;
	vector<int> moveVertices;
	vector<int> fixedVertices;

	Simulation cSim;

	//fix vertices
	for(int i=0; i<cTV.rows(); i++){
		if(cTV.row(i)[1]<=41 && cTV.row(i)[1]>30){
			fixedVertices.push_back(i);
		}
	}
	cSim.initializeSimulation(timestep, 1, method, cTT, cTV, cB, moveVertices, fixedVertices, 2e6, 0.35);
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
	system(("(echo 'Finished Test On:'; date;)>>"+printToHere+"log.txt").c_str());
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

void ConsistencyTest::replaceWithMain()
{
	cout<<"########START CONSISTENCY TESTS######"<<endl;

	cout<<"----Constants for the test"<<endl;
	//Const across all tests
	printThisOften = 0.1;
	printForThisManySeconds = 10;
	timeIterations = 1;
	gravity = -10;
	rayleighCoeff = 0;
	material_model = "neo";
	solver = "newton";
	objectName = "spring";


	int spaceStep = 0;
	char method = 'i';
	string tetmesh_code = "-pRq1.7";
	string spaceDescription = tetmesh_code;
	double implicitTimestep = 1e-1;
	cout<<"----Constants for the test"<<endl;

	//Time stuff------
 	time_t now = time(0);
 	string dt = ctime(&now);//local time, replace all spaces and new lines
 	dt.erase('\n');
 	replace(dt.begin(), dt.end(), ' ', '-');
 	//-----------------

 	MatrixXd V;
 	MatrixXi F;
 	MatrixXd B;
 	MatrixXd TV;
 	MatrixXi TT;
 	MatrixXi TF;

 	igl::readOBJ(TUTORIAL_SHARED_PATH "shared/"+objectName+".obj", V, F);
	igl::copyleft::tetgen::tetrahedralize(V,F,tetmesh_code, TV, TT, TF);
	igl::barycenter(TV, TT, B);
	spaceStep =TT.rows();
	
	if(method == 'i'){
		string printHere = CONSISTENCY_TEST_SAVE_PATH"TestsResults/ConsistencyTests/objectName:"+objectName+"/implicit_euler@"+material_model+"@"+solver+"/"+to_string(spaceStep)+"tets@"+tetmesh_code+"/timestep:"+to_string(implicitTimestep)+"/";
		test(implicitTimestep, method, printHere, TT, TV, B);	
	}else if(method == 'n'){
		cout<<"Not Yet, fill in the same thing for newmark"<<endl;
	}else if(method = 'e'){
		cout<<"Use an implicit method instead"<<endl;
	}else{
		cout<<"Lolz"<<endl;
		exit(0);
	}
	
	return;
}
