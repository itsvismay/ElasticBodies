#include <iostream>
#include <vector>
#include <pthread.h>
#include <fstream>
#include <string>
#include <math.h>
#include <lbfgs.h>
#include <set>
#include "Eigen/SPQRSupport"
#include <Eigen/CholmodSupport>

#include "simulation.h"
#include "globals.h"



using namespace Eigen;
using namespace std;





Simulation::Simulation(void){}

void Simulation::initializeSimulation(double deltaT, char method, MatrixXi& TT, MatrixXd& TV, vector<int> map){
	if (method =='e'){
		integrator = new Verlet();		
	}else if(method == 'i'){
		// integrator = new ImplicitEuler();
	}
	else if(method == 'n'){
		integrator = new ImplicitNewmark();
	}
	else{
		cout<<"Not supported yet"<<endl;
		exit(0);
	}
	mapV2TV = map;
	M.initializeMesh(TT, TV);

	integrator->initializeIntegrator(deltaT, M, TV);	
}

void Simulation::render(){

	integrator->render();

}