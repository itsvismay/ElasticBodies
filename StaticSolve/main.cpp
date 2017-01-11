#include "simulation.h"

#include <igl/viewer/Viewer.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/barycenter.h>
#include <igl/slice.h>

ofstream momentumFile;
ofstream energyFile;
ofstream strainEnergyFile;
ofstream kineticEnergyFile;
ofstream gravityEnergyFile;
ofstream optimizationFile;

double rayleighCoeff;
double gravity;
string material_model;
bool headless;
double youngs;
double poissons;
string tetgen_code;
string solver;
string objectName;


// Input polygon
Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd B;

// Tetrahedralized interior
Eigen::MatrixXd ForcesTV; 
Eigen::MatrixXd FixedTV;
Eigen::MatrixXd TV;
Eigen::MatrixXi TT;
Eigen::MatrixXi TF;


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
	viewer.data.add_edges(RowVector3d(0,0,0), RowVector3d(0,100,0), RowVector3d(0,0,0));
	//viewer.data.add_edges(RowVector3d(0,0,0), RowVector3d(0,0,100), RowVector3d(1,1,1));
	
	return false;
}

bool drawLoop(igl::viewer::Viewer& viewer){
	// Sim.render();

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
		V_temp.row(i*4+0) = Sim.integrator->TV.row(Sim.integrator->TT(s[i],0));
		V_temp.row(i*4+1) = Sim.integrator->TV.row(Sim.integrator->TT(s[i],1));
		V_temp.row(i*4+2) = Sim.integrator->TV.row(Sim.integrator->TT(s[i],2));
		V_temp.row(i*4+3) = Sim.integrator->TV.row(Sim.integrator->TT(s[i],3));
		F_temp.row(i*4+0) << (i*4)+0, (i*4)+1, (i*4)+3;
		F_temp.row(i*4+1) << (i*4)+0, (i*4)+2, (i*4)+1;
		F_temp.row(i*4+2) << (i*4)+3, (i*4)+2, (i*4)+0;
		F_temp.row(i*4+3) << (i*4)+1, (i*4)+2, (i*4)+3;
	}
	viewer.data.clear();
	
	viewer.data.add_edges(RowVector3d(0,0,0), RowVector3d(20,0,0), RowVector3d(1,0,0));
	viewer.data.add_edges(RowVector3d(0,0,0), RowVector3d(0,40,0), RowVector3d(0,1,0));
	viewer.data.add_edges(RowVector3d(0,0,0), RowVector3d(0,0,60), RowVector3d(0,0,1));
	
	// int totalRows = Sim.integrator->TV.rows();
	// for (int i=0; i<Sim.putForceOnTheseVerts.rows(); i++) {
	// 	ForcesTV.row(i) = Sim.integrator->TV.row(Sim.putForceOnTheseVerts(i));
	// }

	viewer.data.add_points(ForcesTV, RowVector3d(1,0,0));
	viewer.data.add_points(FixedTV, RowVector3d(0,1,0));
	viewer.data.set_mesh(V_temp,F_temp);
	viewer.data.set_face_based(true);
	return false;
}

void useFullObject(bool headless, double timestep, int iterations, char method){
	// Load a surface mesh
	// igl::readOFF(TUTORIAL_SHARED_PATH "shared/"+objectName+".off", V, F);
    igl::readOBJ(TUTORIAL_SHARED_PATH "shared/"+objectName+".obj", V, F);
	// Tetrahedralize the interior
	igl::copyleft::tetgen::tetrahedralize(V,F, tetgen_code, TV,TT,TF);

	
	vector<int> moveVertices;
	vector<int> fixedVertices;

	//***************************
	Sim.initializeSimulation(timestep,iterations, method, TT, TV, B, moveVertices, fixedVertices, youngs, poissons);
	ForcesTV.resize(Sim.putForceOnTheseVerts.rows(), 3);
	cout<<Sim.putForceOnTheseVerts.rows()<<endl;
	
	FixedTV.resize(Sim.integrator->fixedVerts.size(), 3);
	FixedTV.setZero();
	for(int i=0; i<Sim.integrator->fixedVerts.size(); i++){
		FixedTV.row(i) = Sim.integrator->TV.row(Sim.integrator->fixedVerts[i]);
	}
	// cout<<"Fixed"<<endl;
	// cout<<FixedTV<<endl;
	if(headless){
		Sim.headless();
	}else{
		igl::viewer::Viewer viewer;
		viewer.callback_pre_draw = &drawLoop;
		viewer.launch();
	}

}

void useMyObject(bool headless, double timestep, int iterations, char method){
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
	TV_One_G << 10, 10, 0, //affect
				0, 20, 0,
				0, 10, 10,
				0, 10, 0,
				0, 0, 0;
				
	// TV_One_G << 0, 0, 0, //affect
	// 			1, 10, 0,
	// 			2, 0, 10,
	// 			3, 0, 0,
	// 			4, -10, 0;

	vector<int> moveVertices;
	vector<int> fixedVertices;

	// move vertices
	//moveVertices.push_back(0);

	// fix vertices
	// fixedVertices.push_back(1);

				
	Sim.initializeSimulation(timestep, iterations, method, TT_One_G, TV_One_G, B, moveVertices, fixedVertices, youngs, poissons);
	
	if(headless){
		Sim.headless();
	}else{
		igl::viewer::Viewer viewer;
		viewer.core.is_animating = true;
		viewer.callback_pre_draw = &drawLoopTest;
		viewer.launch();
	}
		
}



int main(int argc, char *argv[])
{
	double timestep = 0;
	char method;
	int object = 0;
	char runHeadless;
	int iterations;
	string line;

	ifstream configFile (HOME_SAVED_PATH "StaticSolve/config.txt");
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

		getline(configFile, line);
		youngs = stod(line.c_str());
		cout<<youngs<<endl;

		getline(configFile, line);
		poissons = stod(line.c_str());
		cout<<poissons<<endl;

		getline(configFile, line);
		tetgen_code = line.c_str();
		cout<<tetgen_code<<endl;

		getline(configFile, line);
		material_model = line.c_str();
		cout<<material_model<<endl;

		getline(configFile, line);
		solver = line.c_str();
		cout<<solver<<endl;

		getline(configFile, line);
		objectName = line.c_str();
		cout<<objectName<<endl;
	}else{
		cout<<"Check yo self: Config file not found"<<endl;
		return 0;
	}
    
	///////////////////
	cout<<"###########################My Code ###################"<<endl;
	headless = false;
	if(runHeadless=='t'){
		headless = true;
	}
	momentumFile.open("../Scripts/momentum.txt");
	energyFile.open("../Scripts/energy.txt");
	strainEnergyFile.open("../Scripts/senergy.txt");
	kineticEnergyFile.open("../Scripts/kenergy.txt");
	gravityEnergyFile.open("../Scripts/genergy.txt");
	optimizationFile.open("../TestsResults/opt.txt");
	
	if(object ==0){
		useMyObject(headless, timestep, iterations, method);	
	}else if(object ==1){
		useFullObject(headless, timestep, iterations, method);
	}else if(object == 2){
		//consistencyTests(timestep, iterations, method);
	}else{
		cout<<"What do you want to run?"<<endl;
	}
	energyFile.close();
	strainEnergyFile.close();
	kineticEnergyFile.close();
	gravityEnergyFile.close();
	momentumFile.close();
	optimizationFile.close();
	cout<<"###########################My Code ###################"<<endl;
	

	return 0;
}
