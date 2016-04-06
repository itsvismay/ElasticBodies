#include <iostream>
#include <vector>
#include <pthread.h>
#include <fstream>
#include <string>
#include <math.h>
#include "simulation.h"
#include <lbfgs.h>
#include <set>
#include "Eigen/SPQRSupport"
#include <Eigen/CholmodSupport>
#include "globals.h"


using namespace Eigen;
using namespace std;
typedef Eigen::Triplet<double> Trip;
typedef Matrix<double, 12, 1> Vector12d;


Simulation::Simulation(void){}

void Simulation::initializeSimulation(double deltaT, char method, MatrixXi& TT_One, MatrixXd& TV_One, vector<int>& map){
	integrator = method;
	timestep = deltaT;
	mapV2TV = map;
	TV = TV_One;
	vertices = TV_One.rows();
	TVk.resize(vertices, 3);
	ZeroMatrix.resize(3*vertices, 3*vertices);
	ZeroMatrix.setZero();
	global_gradForces.resize(3*vertices, 3*vertices);
	forceGradient.resize(3*vertices, 3*vertices);
	Ident = MatrixXd::Identity(3*vertices, 3*vertices).sparseView();
	grad_g.resize(3*vertices, 3*vertices);

	x_k.resize(3*vertices);
	x_k.setZero();
	v_k.resize(3*vertices);
	v_k.setZero();

	x_old.resize(3*vertices);
	x_old.setZero();
	v_old.resize(3*vertices);
	v_old.setZero();
	v_old(0) =1;
	// v_old(1) =1;
	v_old(2) =1;
	v_old(3) =1;
	// v_old = VectorXd::Random(3*vertices);

	f.resize(3*vertices);
	f.setZero();
	
	M.initializeMesh(TT_One, TV);	
	createInvMassMatrix(); //creates InvMass
	createXFromTet(); //creates x_old
}

void Simulation::createInvMassMatrix(){
	vertex_masses.resize(3*vertices);
	vertex_masses.setZero();

	InvMass.resize(3*vertices, 3*vertices);
	InvMass.setZero();

	for(unsigned int i=0; i<M.tets.size(); i++){
		double vol = (M.tets[i].undeformedVol/4);// TODO: add density 1/Mass for inv matrix,  assume density is 1, so volume=mass
		Vector4i indices = M.tets[i].verticesIndex;
		//TODO use segment here

		vertex_masses(3*indices(0)) += vol;
		vertex_masses(3*indices(0)+1) += vol;
		vertex_masses(3*indices(0)+2) += vol;

		vertex_masses(3*indices(1)) += vol;
		vertex_masses(3*indices(1)+1) += vol;
		vertex_masses(3*indices(1)+2) += vol;

		vertex_masses(3*indices(2)) += vol;
		vertex_masses(3*indices(2)+1) += vol;
		vertex_masses(3*indices(2)+2) += vol;

		vertex_masses(3*indices(3)) += vol;
		vertex_masses(3*indices(3)+1) += vol;
		vertex_masses(3*indices(3)+2) += vol;
	}

	for(int i=0; i<vertices*3; i++){
		InvMass.coeffRef(i, i) = 1/vertex_masses(i);
	}
	cout<<vertex_masses<<endl;

}

void Simulation::fixVertices(int fixed){

	fixedVertices.push_back(fixed);
	//TODO use segments here too
	vertex_masses(3*fixed) = 1000000000000;
	vertex_masses(3*fixed+1) = 1000000000000;
	vertex_masses(3*fixed+2) = 1000000000000;

	InvMass.coeffRef(3*fixed, 3*fixed) = 0;
	InvMass.coeffRef(3*fixed+1, 3*fixed+1) = 0;
	InvMass.coeffRef(3*fixed+2, 3*fixed+2) = 0;
	v_old.segment<3>(3*fixed)*=0;

}

void Simulation::createXFromTet(){
	x_old.setZero();
	for(unsigned int i = 0; i < M.tets.size(); i++){
		Vector4i indices = M.tets[i].verticesIndex;

		x_old(3*indices(0)) = TV.row(indices(0))[0];
		x_old(3*indices(0)+1) = TV.row(indices(0))[1];
		x_old(3*indices(0)+2) = TV.row(indices(0))[2];

		x_old(3*indices(1)) = TV.row(indices(1))[0];
		x_old(3*indices(1)+1) = TV.row(indices(1))[1];
		x_old(3*indices(1)+2) = TV.row(indices(1))[2];

		x_old(3*indices(2)) = TV.row(indices(2))[0];
		x_old(3*indices(2)+1) = TV.row(indices(2))[1];
		x_old(3*indices(2)+2) = TV.row(indices(2))[2];

		x_old(3*indices(3)) = TV.row(indices(3))[0];
		x_old(3*indices(3)+1) = TV.row(indices(3))[1];
		x_old(3*indices(3)+2) = TV.row(indices(3))[2];
	}
}

void Simulation::setXIntoTV(VectorXd& x_new){
	for(unsigned int i=0; i < M.tets.size(); i++){
		Vector4i indices = M.tets[i].verticesIndex;
		TV.row(indices(0)) = Vector3d(x_new(3*indices(0)), x_new(3*indices(0)+1), x_new(3*indices(0) +2));
		TV.row(indices(1)) = Vector3d(x_new(3*indices(1)), x_new(3*indices(1)+1), x_new(3*indices(1) +2));
		TV.row(indices(2)) = Vector3d(x_new(3*indices(2)), x_new(3*indices(2)+1), x_new(3*indices(2) +2));
		TV.row(indices(3)) = Vector3d(x_new(3*indices(3)), x_new(3*indices(3)+1), x_new(3*indices(3) +2)); 
	}
}

void Simulation::createForceVector(){
	f.setZero();
	calculateGravity();
	calculateElasticForce();
}

void Simulation::calculateGravity(){
	for(unsigned int i=0; i<M.tets.size(); i++){
		double vertex_mass = M.tets[i].undeformedVol/4;//assume const density 1
		Vector4i indices = M.tets[i].verticesIndex;

		f(3*indices(0)+1) += vertex_mass*gravity;

		f(3*indices(1)+1) += vertex_mass*gravity; 

		f(3*indices(2)+1) += vertex_mass*gravity;

		f(3*indices(3)+1) += vertex_mass*gravity;
	}
}

void Simulation::calculateElasticForce(){
	for(unsigned int i=0; i<M.tets.size(); i++){
		Vector4i indices = M.tets[i].verticesIndex;

		MatrixXd F_tet = M.tets[i].computeElasticForces(TV, t%2);
		f.segment<3>(3*indices(0)) += F_tet.col(0);
		f.segment<3>(3*indices(1)) += F_tet.col(1);
		f.segment<3>(3*indices(2)) += F_tet.col(2);
		f.segment<3>(3*indices(3)) += F_tet.col(3);
	}
}

void Simulation::ImplicitXtoTV(VectorXd& x_tv, MatrixXd& TVk){
	TVk.setZero();
	for(unsigned int i=0; i < M.tets.size(); i++){
		Vector4i indices = M.tets[i].verticesIndex;
		TVk.row(indices(0)) = Vector3d(x_tv(3*indices(0)), x_tv(3*indices(0)+1), x_tv(3*indices(0) +2));
		TVk.row(indices(1)) = Vector3d(x_tv(3*indices(1)), x_tv(3*indices(1)+1), x_tv(3*indices(1) +2));
		TVk.row(indices(2)) = Vector3d(x_tv(3*indices(2)), x_tv(3*indices(2)+1), x_tv(3*indices(2) +2));
		TVk.row(indices(3)) = Vector3d(x_tv(3*indices(3)), x_tv(3*indices(3)+1), x_tv(3*indices(3) +2)); 
	}
	return;
}

void Simulation::ImplicitTVtoX(VectorXd& x_tv, MatrixXd& TVk){
	x_tv.setZero();
	for(unsigned int i = 0; i < M.tets.size(); i++){
		Vector4i indices = M.tets[i].verticesIndex;

		x_tv(3*indices(0)) = TVk.row(indices(0))[0];
		x_tv(3*indices(0)+1) = TVk.row(indices(0))[1];
		x_tv(3*indices(0)+2) = TVk.row(indices(0))[2];

		x_tv(3*indices(1)) = TVk.row(indices(1))[0];
		x_tv(3*indices(1)+1) = TVk.row(indices(1))[1];
		x_tv(3*indices(1)+2) = TVk.row(indices(1))[2];

		x_tv(3*indices(2)) = TVk.row(indices(2))[0];
		x_tv(3*indices(2)+1) = TVk.row(indices(2))[1];
		x_tv(3*indices(2)+2) = TVk.row(indices(2))[2];

		x_tv(3*indices(3)) = TVk.row(indices(3))[0];
		x_tv(3*indices(3)+1) = TVk.row(indices(3))[1];
		x_tv(3*indices(3)+2) = TVk.row(indices(3))[2];
	}
}

void Simulation::ImplicitCalculateForces( MatrixXd& TVk, SparseMatrix<double>& forceGradient, VectorXd& x_k, VectorXd& f){
	// //gravity
	f.setZero();
	for(unsigned int i=0; i<M.tets.size(); i++){
		double vertex_mass = M.tets[i].undeformedVol/4;//assume const density 1
		Vector4i indices = M.tets[i].verticesIndex;
		f(3*indices(0)+1) += vertex_mass*gravity;
		f(3*indices(1)+1) += vertex_mass*gravity; 
		f(3*indices(2)+1) += vertex_mass*gravity;
		f(3*indices(3)+1) += vertex_mass*gravity;
	}

	//elastic
	for(unsigned int i=0; i<M.tets.size(); i++){
		Vector4i indices = M.tets[i].verticesIndex;
		MatrixXd F_tet = M.tets[i].computeElasticForces(TVk, t%2);
		f.segment<3>(3*indices(0)) += F_tet.col(0);
		f.segment<3>(3*indices(1)) += F_tet.col(1);
		f.segment<3>(3*indices(2)) += F_tet.col(2);
		f.segment<3>(3*indices(3)) += F_tet.col(3);
	}
	// cout<<f<<endl<<endl;
	//damping
	f += rayleighCoeff*forceGradient*(x_k - x_old)/timestep;
	// cout<<f<<endl<<endl;
	return;
}

void Simulation::ImplicitCalculateElasticForceGradient(MatrixXd& TVk, SparseMatrix<double>& forceGradient){
	forceGradient.setZero();
	
	vector<Trip> triplets1;
	triplets1.reserve(3*vertices*3*vertices);	
	for(unsigned int i=0; i<M.tets.size(); i++){
		//Get P(dxn), dx = [1,0, 0...], then [0,1,0,....], and so on... for all 4 vert's x, y, z
		//P is the compute Force Differentials blackbox fxn

		Vector12d dx(12);
		dx.setZero();
		Vector4i indices = M.tets[i].verticesIndex;
		int kj;
		for(unsigned int j=0; j<12; j++){
			dx(j) = 1;
			MatrixXd dForces = M.tets[i].computeForceDifferentials(TVk, dx);
			kj = j%3;
			//row in order for dfxi/dxi ..dfxi/dzl
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[0], dForces(0,0)));
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[0]+1, dForces(1,0)));
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[0]+2, dForces(2,0)));

			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[1], dForces(0,1)));
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[1]+1, dForces(1,1)));
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[1]+2, dForces(2,1)));

			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[2], dForces(0,2)));
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[2]+1, dForces(1,2)));
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[2]+2, dForces(2,2)));

			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[3], dForces(0,3)));
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[3]+1, dForces(1,3)));
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[3]+2, dForces(2,3)));
			dx(j) = 0; //ASK check is this efficient?
		}
	}
	forceGradient.setFromTriplets(triplets1.begin(), triplets1.end());

	return;
}
void Simulation::ImplicitXfromTV(VectorXd& x_n1, MatrixXd& TVk){
	x_n1.setZero();
	for(unsigned int i = 0; i < M.tets.size(); i++){
		Vector4i indices = M.tets[i].verticesIndex;
		x_n1.segment<3>(indices(0)) = TVk.row(indices(0));
		x_n1.segment<3>(indices(1)) = TVk.row(indices(1));
		x_n1.segment<3>(indices(2)) = TVk.row(indices(2));
		x_n1.segment<3>(indices(3)) = TVk.row(indices(3));
	}
}
void Simulation::renderExplicit(){
	t+=1;
	//Explicit Code
	x_old = x_old + timestep*v_old;
	if(x_k != x_k){
		cout<<"NAN"<<endl;
		exit(0);
	}
	setXIntoTV(x_old);
	createForceVector();
	v_old = v_old + timestep*InvMass*f;

}
void Simulation::renderImplicit(){
	t+=1;
	//Implicit Code
	// cout<<"diff in x"<<endl;
	// cout<<x_k - x_old<<endl<<endl;
	v_k.setZero();
	x_k.setZero();
	x_k = x_old;
	v_k = v_old;

	forceGradient.setZero();
	bool Nan=false;
	int NEWTON_MAX = 100, i =0;
	cout<<"--------"<<t<<"-------"<<endl;
	cout<<"x_k"<<endl;
	cout<<x_k<<endl<<endl;
	cout<<"v_k"<<endl;
	cout<<v_k<<endl<<endl;
	cout<<"--------------------"<<endl;
	for( i=0; i<NEWTON_MAX; i++){
		grad_g.setZero();
	
		ImplicitXtoTV(x_k, TVk);//TVk value changed in function
		ImplicitCalculateElasticForceGradient(TVk, forceGradient); 
		ImplicitCalculateForces(TVk, forceGradient, x_k, f);

		VectorXd g = x_k - x_old -timestep*v_old -timestep*timestep*InvMass*f;
		grad_g = Ident - timestep*timestep*InvMass*forceGradient - timestep*rayleighCoeff*InvMass*forceGradient;
		
		// cout<<"G"<<t<<endl;
		// cout<<g<<endl<<endl;
		// cout<<"G Gradient"<<t<<endl;
		// cout<<grad_g<<endl;

		//solve for delta x
		// Conj Grad
		// ConjugateGradient<SparseMatrix<double>> cg;
		// cg.compute(grad_g);
		// VectorXd deltaX = -1*cg.solve(g);

		// Sparse Cholesky LL^T
		// SimplicialLLT<SparseMatrix<double>> llt;
		// llt.compute(grad_g);
		// VectorXd deltaX = -1* llt.solve(g);

		//Sparse QR 
		SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> sqr;
		sqr.compute(grad_g);
		VectorXd deltaX = -1*sqr.solve(g);

		// CholmodSimplicialLLT<SparseMatrix<double>> cholmodllt;
		// cholmodllt.compute(grad_g);
		// VectorXd deltaX = -cholmodllt.solve(g);
		

		x_k+=deltaX;
		if(x_k != x_k){
			Nan = true;
			break;
		}
		if(g.squaredNorm()<.00000001){
			break;
		}
	}
	if(Nan){
		cout<<"ERROR: Newton's method doesn't converge"<<endl;
		cout<<i<<endl;
		exit(0);
	}
	if(i== NEWTON_MAX){
		cout<<"ERROR: Newton max reached"<<endl;
		cout<<i<<endl;
		exit(0);
	}
	v_old.setZero();
	v_old = (x_k - x_old)/timestep;
	x_old = x_k;
	cout<<"*******************"<<endl;
	cout<< "New Pos"<<t<<endl;
	cout<<x_old<<endl<<endl;
	cout<< "New Vels"<<t<<endl;
	cout<<v_old<<endl;
	cout<<"*****************"<<endl<<endl;
	ImplicitXtoTV(x_old, TV);

	
}

//TODO: Optimize this using hashing
bool Simulation::isFixed(int vert){
	for(unsigned int j=0; j<fixedVertices.size(); j++){
		if(vert == fixedVertices[j]){
			return true;
		}
	}
	return false;
}

void Simulation::render(){
	if(integrator == 'e'){
		cout<<endl;
		cout<<'e'<<t<<endl;
		renderExplicit();
	}else{
		cout<<endl;
		cout<<'i'<<t<<endl;
		renderImplicit();
	}

	////////////////////////////////////
	double TotalEnergy = 0;
	double gravityE =0;
	double kineticE =0;
	double strainE = 0;
	for(int i=0; i<vertices; i++){
		if(!isFixed(i)){
			int k=3*i;
			gravityE +=  vertex_masses(k)*-1*gravity*(x_old(k));
			kineticE += 0.5*vertex_masses(k)*v_old(k)*v_old(k);
			k++;
			gravityE +=  vertex_masses(k)*-1*gravity*(x_old(k));
			kineticE += 0.5*vertex_masses(k)*v_old(k)*v_old(k);
			k++;
			gravityE +=  vertex_masses(k)*-1*gravity*(x_old(k));
			kineticE += 0.5*vertex_masses(k)*v_old(k)*v_old(k);
		}		
	}
	for(unsigned int i=0; i<M.tets.size(); i++){
		strainE += M.tets[i].undeformedVol*M.tets[i].energyDensity;		
	}

	TotalEnergy+= gravityE + kineticE + strainE;
	// cout<<endl<<"Grav E"<<endl;
	// cout<<strainE<<endl;
	cout<<"Tot E"<<endl;
	cout<<TotalEnergy<<endl;
	cout<<"Strain E"<<endl;
	cout<<strainE<<endl;
	energyFile<<t<<", "<<TotalEnergy<<"\n";
	strainEnergyFile<<t<<", "<<strainE<<"\n";
	kineticEnergyFile<<t<<", "<<kineticE<<"\n";
	gravityEnergyFile<<t<<", "<<gravityE<<"\n";

	////////////////////////////////////
	
	// //////Momentum Code////////////
	// if(momentumFile.is_open()){
	// 	double xp=0;
	// 	double yp=0;
	// 	double zp=0;
	// 	cout<<"-----------"<<endl;
	// 	for(int i=0; i<v_old.rows(); ){
	// 		xp += v_old(i)*vertex_masses(i);
	// 		cout<<"x"<<endl;
	// 		cout<<v_old(i)<<endl;
	// 		i++;
	// 		yp += v_old(i)*vertex_masses(i);
	// 		cout<<"y"<<endl;
	// 		cout<<v_old(i)<<endl;
	// 		i++;
	// 		zp += v_old(i)*vertex_masses(i);
	// 		cout<<"z"<<endl;
	// 		cout<<v_old(i)<<endl;
	// 		i++;
	// 	}
	// 	cout<<"- - - - -- "<<endl;
	// 	momentumFile<<t<<","<<xp<<","<<yp<<","<<zp<<"\n";
	// }else{
	// 	cout<<"no open file"<<endl;
	// }
	// /////////////////
}