#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <pthread.h>
#include <fstream>
#include <math.h>

#include "Verlet.h"
#include "globals.h"

using namespace Eigen;
using namespace std;

void Verlet::initializeIntegrator(double ph, SolidMesh& pM, MatrixXd& pTV){
	cout<<"Here"<<endl;
	IntegratorAbstract::initializeIntegrator(ph, pM, pTV);
}

void Verlet::setXIntoTV(VectorXd& x_new){
	for(unsigned int i=0; i < M.tets.size(); i++){
		Vector4i indices = M.tets[i].verticesIndex;
		TV.row(indices(0)) = Vector3d(x_new(3*indices(0)), x_new(3*indices(0)+1), x_new(3*indices(0) +2));
		TV.row(indices(1)) = Vector3d(x_new(3*indices(1)), x_new(3*indices(1)+1), x_new(3*indices(1) +2));
		TV.row(indices(2)) = Vector3d(x_new(3*indices(2)), x_new(3*indices(2)+1), x_new(3*indices(2) +2));
		TV.row(indices(3)) = Vector3d(x_new(3*indices(3)), x_new(3*indices(3)+1), x_new(3*indices(3) +2)); 
	}
}

void Verlet::createForceVector(){
	f.setZero();
	calculateGravity();
	calculateElasticForce();
}

void Verlet::calculateGravity(){
	for(unsigned int i=0; i<M.tets.size(); i++){
		double vertex_mass = M.tets[i].undeformedVol/4;//assume const density 1
		Vector4i indices = M.tets[i].verticesIndex;

		f(3*indices(0)+1) += vertex_mass*gravity;
		f(3*indices(1)+1) += vertex_mass*gravity; 
		f(3*indices(2)+1) += vertex_mass*gravity;
		f(3*indices(3)+1) += vertex_mass*gravity;
	}
}

void Verlet::calculateElasticForce(){
	for(unsigned int i=0; i<M.tets.size(); i++){
		Vector4i indices = M.tets[i].verticesIndex;

		MatrixXd F_tet = M.tets[i].computeElasticForces(TV, simTime%2);
		f.segment<3>(3*indices(0)) += F_tet.col(0);
		f.segment<3>(3*indices(1)) += F_tet.col(1);
		f.segment<3>(3*indices(2)) += F_tet.col(2);
		f.segment<3>(3*indices(3)) += F_tet.col(3);
	}
}

void Verlet::render(){
	simTime+=1;
	cout<<"e"<<simTime<<endl;
	// cout<<"XOld"<<endl;
	// cout<<x_old<<endl;
	// cout<<"VOld"<<endl;
	// cout<<v_old<<endl;
	// cout<<"fOld"<<endl;
	// cout<<f<<endl;

	x_old = x_old + h*v_old;
	if(x_old != x_old){
		cout<<"NAN Exiting"<<endl;
		exit(0);
	}
	setXIntoTV(x_old);
	createForceVector();
	v_old = v_old + h*InvMass*f;

	IntegratorAbstract::printInfo();
}
