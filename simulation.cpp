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
		integrator = new ImplicitEuler();
	}else if(method == 'n'){
		integrator = new ImplicitNewmark();
	}else{
		cout<<"Not supported yet"<<endl;
		exit(0);
	}
	mapV2TV = map;
	M.initializeMesh(TT, TV);

	integrator->initializeIntegrator(deltaT, M, TV);	
}


// }


// void Simulation::renderNewmark(){
// 	t+=1;
// 	v_k.setZero();
// 	x_k.setZero();
// 	x_k = x_old;
// 	forceGradient.setZero();
// 	bool Nan = false;
// 	int NEWTON_MAX = 100, i=0;
// 	double gamma = 0.5;
// 	double beta =0.25;
// 	VectorXd f_old = f;
// 	for(i=0; i<NEWTON_MAX; i++){
// 		grad_g.setZero();
// 		ImplicitXtoTV(x_k, TVk);
// 		ImplicitCalculateElasticForceGradient(TVk, forceGradient);
// 		ImplicitCalculateForces(TVk, forceGradient, v_k, f);
		
// 		VectorXd g = x_k - x_old - timestep*v_old - (timestep*timestep/2)*(1-2*beta)*InvMass*f_old - (timestep*timestep*beta)*InvMass*f;
// 		grad_g = Ident - timestep*timestep*beta*InvMass*(forceGradient+(rayleighCoeff/timestep)*forceGradient);
		
// 		SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> sqr;
// 		sqr.compute(grad_g);
// 		VectorXd deltaX = -1*sqr.solve(g);
// 		// v_k+=deltaX/timestep;
// 		x_k+=deltaX;
// 		//v_k = (x_k- x_old)/timestep;
// 		if(x_k != x_k){
// 			Nan = true;
// 			break;
// 		}
// 		if(g.squaredNorm()<.00000001){
// 			break;
// 		}
// 	}
// 	if(Nan){
// 		cout<<"ERROR NEWMARK: Newton's method doesn't converge"<<endl;
// 		cout<<i<<endl;
// 		exit(0);
// 	}
// 	if(i== NEWTON_MAX){
// 		cout<<"ERROR NEWMARK: Newton max reached"<<endl;
// 		cout<<i<<endl;
// 		exit(0);
// 	}
// 	// v_old.setZero();
// 	v_old = v_old + timestep*(1-gamma)*InvMass*f_old + timestep*gamma*InvMass*f;
// 	x_old = x_k;
// 	cout<<"*******************"<<endl;
// 	cout<< "New Pos"<<t<<endl;
// 	cout<<x_old<<endl<<endl;
// 	cout<< "New Vels"<<t<<endl;
// 	cout<<v_old<<endl;
// 	cout<<"*****************"<<endl<<endl;
// 	ImplicitXtoTV(x_old, TV);
// }


void Simulation::render(){

	integrator->render();

}