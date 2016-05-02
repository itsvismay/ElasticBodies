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

int Simulation::initializeSimulation(double deltaT, int iterations, char method, MatrixXi& TT, MatrixXd& TV){
	iters = iterations;
	if (method =='e'){
		integrator = new Verlet();
		cout<<"Initialized Verlet"<<endl;	
	}else if(method == 'i'){
		integrator = new ImplicitEuler();
		cout<<"Initialized Implicit Euler"<<endl;
	}
	else if(method == 'n'){
		integrator = new ImplicitNewmark();
		cout<<"Initialized Implicit Newmark"<<endl;
	}
	else{
		cout<<"Method not supported yet"<<endl;
		exit(0);
	}
	// mapV2TV = map;
	M.initializeMesh(TT, TV);

	integrator->initializeIntegrator(deltaT, M, TV);

	return 1;
}

void Simulation::headless(){
	while(integrator->simTime<iters){
		integrator->render();
	}
}
void Simulation::render(){
	integrator->render();
}

