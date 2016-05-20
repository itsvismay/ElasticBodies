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

int Simulation::initializeSimulation(double deltaT, int iterations, char method, MatrixXi& TT, MatrixXd& TV, MatrixXd& B){
	iters = iterations;
	B = B;
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
	//Initialize Solid Mesh
	M.initializeMesh(TT, TV);

	integrator->initializeIntegrator(deltaT, M, TV, TT);

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

int Simulation::reIndexClampedVertices(vector<int>& moveVertices){
	//Re-index clamped vertices
	for(int i=0; i<moveVertices.size(); i++){
		//re-index vertices
		//take TV row at index moveVertices(i)
		//replace with TV row at index TV.rows - (moveVertices.size() -i)
		//Vector3d ro1 = TV.row(moveVertices(i));
		cout<<"move"<<moveVertices[i]<<endl;
		Vector3d ro2 = integrator->TV.row(integrator->TV.rows() - moveVertices.size()+ i);
		integrator->TV.row(integrator->TV.rows() - moveVertices.size()+ i) = integrator->TV.row(moveVertices[i]);
		integrator->TV.row(moveVertices[i]) = ro2;

		//re-index all the tet pointers
		int v1 = moveVertices[i];
		int v2 = integrator->TV.rows() - moveVertices.size() +i;
		for(int k=0; k< integrator->TT.rows(); k++){
			for(int j=0; j<4; j++){
				if (integrator->TT.row(k)[j]== v1){
					integrator->TT.row(k)[j] = v2;
				}else if(integrator->TT.row(k)[j]==v2){
					integrator->TT.row(k)[j] = v1;
				}
			}
		}
	}

	return 1;
}

void Simulation::setInitPosition(vector<int> moveVertices){
	reIndexClampedVertices(moveVertices);
}
