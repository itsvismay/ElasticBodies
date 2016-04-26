#include <igl/viewer/Viewer.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/barycenter.h>

#include "simulation.h"
#include "globals.h"

ofstream momentumFile;
ofstream energyFile;
ofstream strainEnergyFile;
ofstream kineticEnergyFile;
ofstream gravityEnergyFile;

double rayleighCoeff;
double gravity;


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
	viewer.data.clear();
	Sim.render();

	for(unsigned int i=0; i<Sim.mapV2TV.size(); i++){
		V.row(i) = Sim.integrator->TV.row(Sim.mapV2TV[i]);
	}

	viewer.data.add_edges(RowVector3d(0,0,0), RowVector3d(100,0,0), RowVector3d(1,1,1));
	viewer.data.add_edges(RowVector3d(0,0,0), RowVector3d(0,100,0), RowVector3d(1,1,1));
	//viewer.data.add_edges(RowVector3d(0,0,0), RowVector3d(0,0,100), RowVector3d(1,1,1));
	
	viewer.data.set_mesh(V, F);
	viewer.data.set_face_based(false);
	return false;
}

void useFullObject(bool headless, double timestep, int iterations, char method){
	vector<int> mapV2TV;
	// Load a surface mesh
	// igl::readOFF(TUTORIAL_SHARED_PATH "/3holes.off",V,F);
	igl::readOBJ(TUTORIAL_SHARED_PATH "/wheel.obj", V, F);
	V = V;
	// Tetrahedralize the interior
	igl::copyleft::tetgen::tetrahedralize(V,F,"-pq2/0", TV,TT,TF);
	// // Compute barycenters
	igl::barycenter(TV, TT,B);
	// //constructs the map from V to TV
	for(unsigned int i=0; i< V.rows(); i++){
		for(unsigned int j=0; j<TV.rows(); j++){
			if( (V.row(i)-TV.row(j)).squaredNorm() <= 0.00001){
				mapV2TV.push_back(j);
				break;
			}
		}
	}
	Sim.initializeSimulation(timestep, method, TT, TV, mapV2TV);
	if(!headless){
		igl::viewer::Viewer viewer;
		viewer.callback_pre_draw = &drawLoop;
		viewer.launch();
	}else{
		while(Sim.integrator->simTime<iterations){
			Sim.render();
			for(unsigned int i=0; i<Sim.mapV2TV.size(); i++){
				V.row(i) = Sim.integrator->TV.row(Sim.mapV2TV[i]);
			}
			if(Sim.integrator->simTime%1== 0){
				cout<<"--"<<endl;
				cout<<int(Sim.integrator->simTime*timestep*100)%10<<endl;
				// igl::writeOBJ("../output1/object"+to_string(Sim.t)+".obj", V, F);
			}
		}
	}
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


				
	Sim.initializeSimulation(timestep, method, TT_One_G, TV_One_G, mapV2TV);
	igl::viewer::Viewer viewer;
	bool boolVariable = true;
	double timeVariable = 0.001;
	if(!headless){
		viewer.callback_pre_draw = &drawLoopTest;
		viewer.launch();
	}else{
		while(Sim.integrator->simTime<iterations){
			Sim.render();
		}
	}
}

void consistencyTests(bool headless, double timestep, int iterations, char method){
	vector<int> mapV2TV;
	// Load a surface mesh
	// igl::readOFF(TUTORIAL_SHARED_PATH "/3holes.off",V,F);
	igl::readOBJ(TUTORIAL_SHARED_PATH "/wheel.obj", V, F);
	V = V;
	// Tetrahedralize the interior
	igl::copyleft::tetgen::tetrahedralize(V,F,"-pq2/0", TV,TT,TF);
	// // Compute barycenters
	igl::barycenter(TV, TT,B);
	// //constructs the map from V to TV
	for(unsigned int i=0; i< V.rows(); i++){
		for(unsigned int j=0; j<TV.rows(); j++){
			if( (V.row(i)-TV.row(j)).squaredNorm() <= 0.00001){
				mapV2TV.push_back(j);
				break;
			}
		}
	}
	Sim.initializeSimulation(timestep, method, TT, TV, mapV2TV);
	if(!headless){
		igl::viewer::Viewer viewer;
		viewer.callback_pre_draw = &drawLoop;
		viewer.launch();
	}else{
		while(Sim.integrator->simTime<iterations){
			Sim.render();
			for(unsigned int i=0; i<Sim.mapV2TV.size(); i++){
				V.row(i) = Sim.integrator->TV.row(Sim.mapV2TV[i]);
			}
			if(Sim.integrator->simTime%1== 0){
				cout<<"--"<<endl;
				cout<<int(Sim.integrator->simTime*timestep*100)%10<<endl;
				// igl::writeOBJ("../output1/object"+to_string(Sim.t)+".obj", V, F);
			}
		}
	}
}

int main(int argc, char *argv[])
{
	double timestep = 0;
	char headless;
	char method;
	int iterations = 0;
	int object = 0;
	string line;
	ifstream configFile ("../config.txt");
	if(configFile.is_open()){
		getline(configFile, line);
		timestep = stod(line.c_str());
		cout<<timestep<<endl;
		
		getline(configFile, line);
		headless = line.c_str()[0];
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
	momentumFile.open("../PythonScripts/momentum.txt");
	energyFile.open("../PythonScripts/energy.txt");
	strainEnergyFile.open("../PythonScripts/senergy.txt");
	kineticEnergyFile.open("../PythonScripts/kenergy.txt");
	gravityEnergyFile.open("../PythonScripts/genergy.txt");
    ///////////////////
	cout<<"###########################My Code ###################"<<endl;
	bool runHeadless = false;
	if(headless=='t'){
		runHeadless = true;
	}
	if(object ==0){
		useMyObject(runHeadless, timestep, iterations, method);	
	}else{
		useFullObject(runHeadless, timestep, iterations, method);
	}
	cout<<"###########################My Code ###################"<<endl;
	momentumFile.close();
	energyFile.close();
	strainEnergyFile.close();
	kineticEnergyFile.close();
	gravityEnergyFile.close();

	return 0;
}

//Get Full SVK and Full Neohookean working


// read lbfgs and api
// houdini

// keep CG
// get eigen sparse solvers (or suitesparse) - sparse cholesky, sparse QR
// try profiler with real lambda and mu values

// test implicit euler with real lambda and mu
// - test with direct solver
// - do consistency test with/without damping
// - - find size of timestep, so as I decrease timestep, trajectory
// - see if I need to use implicit midpoint newmark

// test damping
// - find real damping parameters for spring
// - 