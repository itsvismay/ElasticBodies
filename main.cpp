#include <igl/viewer/Viewer.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/barycenter.h>
#include <fstream>

#include "simulation.h"
#include "globals.h"
#include "ConsistencyTests.h"

ofstream momentumFile;
ofstream energyFile;
ofstream strainEnergyFile;
ofstream kineticEnergyFile;
ofstream gravityEnergyFile;

double rayleighCoeff;
double gravity;
bool headless;


// Input polygon
Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd B;

// Tetrahedralized interior
Eigen::MatrixXd TV;
Eigen::MatrixXi TT;
Eigen::MatrixXi TF;

#ifndef TUTORIAL_SHARED_PATH
#define TUTORIAL_SHARED_PATH "../shared"
#endif

MatrixXi TT_One_G;
MatrixXd TV_One_G;

Simulation Sim;

bool drawLoopTest(igl::viewer::Viewer& viewer){
	viewer.data.clear();
	Sim.render();
	

	viewer.data.add_points(Sim.integrator->TV, RowVector3d(1,0,0));
	viewer.data.add_edges(Sim.integrator->TV.row(0), Sim.integrator->TV.row(1), RowVector3d(1,0,0));
	viewer.data.add_edges(Sim.integrator->TV.row(0), Sim.integrator->TV.row(2), RowVector3d(1,0,0));
	viewer.data.add_edges(Sim.integrator->TV.row(0), Sim.integrator->TV.row(3), RowVector3d(1,0,0));
	viewer.data.add_edges(Sim.integrator->TV.row(1), Sim.integrator->TV.row(2), RowVector3d(1,0,0));
	viewer.data.add_edges(Sim.integrator->TV.row(1), Sim.integrator->TV.row(3), RowVector3d(1,0,0));
	viewer.data.add_edges(Sim.integrator->TV.row(2), Sim.integrator->TV.row(3), RowVector3d(1,0,0));

	viewer.data.add_edges(Sim.integrator->TV.row(4), Sim.integrator->TV.row(3), RowVector3d(0,1,0));
	viewer.data.add_edges(Sim.integrator->TV.row(4), Sim.integrator->TV.row(0), RowVector3d(0,0,1));
	viewer.data.add_edges(Sim.integrator->TV.row(4), Sim.integrator->TV.row(2), RowVector3d(0,0,0));

	// viewer.data.add_edges(Sim.integrator->TV.row(5), Sim.integrator->TV.row(3), RowVector3d(0,1,0));
	// viewer.data.add_edges(Sim.integrator->TV.row(5), Sim.integrator->TV.row(0), RowVector3d(0,0,1));
	// viewer.data.add_edges(Sim.integrator->TV.row(5), Sim.integrator->TV.row(1), RowVector3d(0,0,0));

	// viewer.data.add_edges(Sim.integrator->TV.row(6), Sim.integrator->TV.row(3), RowVector3d(0,1,0));
	// viewer.data.add_edges(Sim.integrator->TV.row(6), Sim.integrator->TV.row(0), RowVector3d(0,0,1));
	// viewer.data.add_edges(Sim.integrator->TV.row(6), Sim.integrator->TV.row(5), RowVector3d(0,0,0));


	viewer.data.add_edges(RowVector3d(0,0,0), RowVector3d(100,0,0), RowVector3d(1,1,1));
	viewer.data.add_edges(RowVector3d(0,0,0), RowVector3d(0,100,0), RowVector3d(1,1,1));
	//viewer.data.add_edges(RowVector3d(0,0,0), RowVector3d(0,0,100), RowVector3d(1,1,1));
	
	return false;
}

bool drawLoop(igl::viewer::Viewer& viewer){
	Sim.render();

	// // for(unsigned int i=0; i<Sim.mapV2TV.size(); i++){
	// // 	V.row(i) = Sim.integrator->TV.row(Sim.mapV2TV[i]);
	// // }

	// viewer.data.add_edges(RowVector3d(0,0,0), RowVector3d(100,0,0), RowVector3d(1,1,1));
	// viewer.data.add_edges(RowVector3d(0,0,0), RowVector3d(0,100,0), RowVector3d(1,1,1));
	// //viewer.data.add_edges(RowVector3d(0,0,0), RowVector3d(0,0,100), RowVector3d(1,1,1));
	
	// viewer.data.clear();
	// viewer.data.set_mesh(V, F);
	// viewer.data.set_face_based(false);

	double refinement = 9;
	double t = ((refinement - 1)+1) / 9.0;

	VectorXd v = B.col(2).array() - B.col(2).minCoeff();
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
		V_temp.row(i*4+0) = Sim.integrator->TV.row(TT(s[i],0));
		V_temp.row(i*4+1) = Sim.integrator->TV.row(TT(s[i],1));
		V_temp.row(i*4+2) = Sim.integrator->TV.row(TT(s[i],2));
		V_temp.row(i*4+3) = Sim.integrator->TV.row(TT(s[i],3));
		F_temp.row(i*4+0) << (i*4)+0, (i*4)+1, (i*4)+3;
		F_temp.row(i*4+1) << (i*4)+0, (i*4)+2, (i*4)+1;
		F_temp.row(i*4+2) << (i*4)+3, (i*4)+2, (i*4)+0;
		F_temp.row(i*4+3) << (i*4)+1, (i*4)+2, (i*4)+3;
	}
	viewer.data.clear();
	
	viewer.data.add_edges(RowVector3d(0,0,0), RowVector3d(100,0,0), RowVector3d(1,1,1));
	viewer.data.add_edges(RowVector3d(0,0,0), RowVector3d(0,100,0), RowVector3d(1,1,1));
	//viewer.data.add_edges(RowVector3d(0,0,0), RowVector3d(0,0,100), RowVector3d(1,1,1));
	viewer.data.set_mesh(V_temp,F_temp);
	viewer.data.set_face_based(true);
	return false;
}

void useFullObject(bool headless, double timestep, int iterations, char method){
	vector<int> mapV2TV;
	// Load a surface mesh
	// igl::readOFF(TUTORIAL_SHARED_PATH "/3holes.off",V,F);
	igl::readOBJ(TUTORIAL_SHARED_PATH "/beam.obj", V, F);

	// Tetrahedralize the interior
	// igl::copyleft::tetgen::tetrahedralize(V,F,"-pq2/0", TV,TT,TF);
	igl::copyleft::tetgen::tetrahedralize(V,F,"pq1.414Y", TV,TT,TF);
	
	// // Compute barycenters
	igl::barycenter(TV, TT,B);

	Sim.initializeSimulation(timestep,iterations, method, TT, TV, B);
	Sim.mapV2TV = mapV2TV;
	//fix vertices
	Sim.integrator->fixVertices(0);
	Sim.integrator->fixVertices(1);
	// Sim.integrator->fixVertices(2);
	if(headless){
		Sim.headless();
	}else{
		igl::viewer::Viewer viewer;
		viewer.callback_pre_draw = &drawLoop;
		viewer.launch();
	}
	// if(!headless){
	// 	igl::viewer::Viewer viewer;
	// 	viewer.callback_pre_draw = &drawLoop;
	// 	viewer.launch();
	// }else{
		// while(integrator->simTime<iterations){
		// 	Sim.render();
		// 	for(unsigned int i=0; i<Sim.mapV2TV.size(); i++){
		// 		V.row(i) = Sim.integrator->TV.row(Sim.mapV2TV[i]);
		// 	}
		// 	if(Sim.integrator->simTime%1== 0){
		// 		cout<<"--"<<endl;
		// 		cout<<int(Sim.integrator->simTime*timestep*100)%10<<endl;
		// 		// igl::writeOBJ("../output1/object"+to_string(Sim.t)+".obj", V, F);
		// 	}
		// }	
	// }
}

void useMyObject(bool headless, double timestep, int iterations, char method){
	vector<int> mapV2TV;


	// TT_One_G.resize(1, 4);
	// TT_One_G<< 1, 2, 3, 0;

	// TV_One_G.resize(4, 3);
	// TV_One_G << 0, 0, 10, //affect
	// 			10, 0, 0,
	// 			0, 10, 0,
	// 			0, 0, 0;

	TT_One_G.resize(2, 4);
	TT_One_G<<  0, 1, 2, 3,
				4, 0, 2, 3;

	TV_One_G.resize(5, 3);
	TV_One_G << 10, 0, 0, //affect
				0, 10, 0,
				0, 0, 10,
				0, 0, 0,
				0, -10, 0;


				
	Sim.initializeSimulation(timestep, iterations, method, TT_One_G, TV_One_G, B);
	Sim.integrator->fixVertices(1);
	if(headless){
		Sim.headless();
	}else{
		igl::viewer::Viewer viewer;
		viewer.callback_pre_draw = &drawLoopTest;
		viewer.launch();
	}
		
}

void consistencyTests( double timestep, int iterations, char method){	

	ConsistencyTest c;
	c.runTimeTests(Sim);
	
}

int main(int argc, char *argv[])
{
	double timestep = 0;
	char method;
	int object = 0;
	char runHeadless;
	int iterations;
	string line;
	ifstream configFile ("../config.txt");
	if(configFile.is_open()){
		getline(configFile, line);
		timestep = stod(line.c_str());
		cout<<timestep<<endl;
		
		getline(configFile, line);
		runHeadless = line.c_str()[0];
		cout<<line<<endl;
		
		getline(configFile, line);
		method = line.c_str()[0];
		cout<<method<<endl;

		getline(configFile, line);
		iterations = atoi(line.c_str());
		cout<<iterations<<endl;

		getline(configFile, line);
		rayleighCoeff = stod(line.c_str());
		cout<<rayleighCoeff<<endl;

		getline(configFile, line);
		gravity = stod(line.c_str());
		cout<<gravity<<endl;

		getline(configFile, line);
		object = atoi(line.c_str());
		cout<<object<<endl;
	}else{
		cout<<"Config file not found"<<endl;
		return 0;
	}
	
    ///////////////////
	cout<<"###########################My Code ###################"<<endl;
	headless = false;
	if(runHeadless=='t'){
		headless = true;
	}
	momentumFile.open("../PythonScripts/momentum.txt");
	energyFile.open("../PythonScripts/energy.txt");
	strainEnergyFile.open("../PythonScripts/senergy.txt");
	kineticEnergyFile.open("../PythonScripts/kenergy.txt");
	gravityEnergyFile.open("../PythonScripts/genergy.txt");
	
	if(object ==0){
		useMyObject(headless, timestep, iterations, method);	
	}else if(object ==1){
		useFullObject(headless, timestep, iterations, method);
	}else if(object == 2){
		consistencyTests(timestep, iterations, method);
	}else{
		cout<<"What do you want to run?"<<endl;
	}
	energyFile.close();
	strainEnergyFile.close();
	kineticEnergyFile.close();
	gravityEnergyFile.close();
	momentumFile.close();
	cout<<"###########################My Code ###################"<<endl;
	

	return 0;
}

// - do consistency test with/without damping
// - - find size of timestep, so as I decrease timestep, trajectory
// - see if I need to use implicit midpoint newmark

// test damping
// - find real damping parameters for spring

// fix mapping  from tetmesh to rendered surface faces
// redo implicit methods fixVertices method (fix velocites every time step)