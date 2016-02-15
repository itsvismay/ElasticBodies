#include <igl/viewer/Viewer.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/barycenter.h>
#include <iostream>
#include <vector>
#include <pthread.h>
#include <fstream>
#include <string>
#include <math.h>
#include "simulation.h"
#include <set>
#ifndef TUTORIAL_SHARED_PATH
#define TUTORIAL_SHARED_PATH "../shared"
#endif

using namespace Eigen;
using namespace std;
typedef Eigen::Triplet<double> Trip;
typedef Matrix<double, 12, 1> Vector12d;
double rayleighCoeff;
double gravity = 0;

// Input polygon
Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd B;
ofstream momentumFile;
ofstream energyFile;
ofstream strainEnergyFile;
ofstream kineticEnergyFile;
ofstream gravityEnergyFile;
// Tetrahedralized interior
Eigen::MatrixXd TV;
Eigen::MatrixXi TT;
Eigen::MatrixXi TF;

class Spring{
public:
	int v1, v2;
	float restLen;
	float stiffness;

	Spring(int vert1, int vert2, float rL);
};

Spring::Spring(int vert1, int vert2, float rL){
		v1 = vert1;
		v2 = vert2;
		stiffness = 500000;
		restLen = rL;
	
}
vector<Spring> springList;
set< pair<int, int> > springSet;

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
	v_old = VectorXd::Random(3*vertices)*10;

	f.resize(3*vertices);
	f.setZero();
	
	M.initializeMesh(TT_One, TV);	
	createInvMassMatrix(); //creates InvMass
	createXFromTet(); //creates x_old

	/////////////////////////////////
	for(int i=0; i<M.tets.size(); i++){
		Vector4i indices = M.tets[i].verticesIndex;
		insertToSpringSet(indices[0],  indices[1]);
		insertToSpringSet(indices[0], indices[2]);
		insertToSpringSet(indices[0] , indices[3]);
		insertToSpringSet(indices[1], indices[2]);
		insertToSpringSet(indices[1], indices[3]);
		insertToSpringSet(indices[2], indices[3]);
		// Spring s1(indices[0], indices[1]);
		// springList.push_back(s1);
		// Spring s2(indices[0], indices[2]);
		// springList.push_back(s2);
		// Spring s3(indices[0], indices[3]);
		// springList.push_back(s3);
		// Spring s4(indices[1], indices[2]);
		// springList.push_back(s4);
		// Spring s5(indices[1], indices[3]);
		// springList.push_back(s5);
		// Spring s6(indices[2], indices[3]);
		// springList.push_back(s6);


	}
	/////////////////////////////////
}

void Simulation::insertToSpringSet(int i1, int i2){
	pair<int, int> pair1(i1,  i2);
	pair<int, int> pair2(i2,  i1);

	pair<set< pair<int,int> >::iterator, bool> ret;
	ret = springSet.insert(pair1);

	if(ret.second==false){
		//item exists in set
		return;
	}
	springSet.insert(pair2);
	Spring s(i1, i2, (TV.row(i2)- TV.row(i1)).norm());
	springList.push_back(s);
	return;

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
	//UNCOMMENT
	// calculateGravity();
	// calculateElasticForce();

	////////////////////////////////////////
	for(int i=0; i<springList.size(); i++){
		int v1 = springList[i].v1;
		int v2 = springList[i].v2;
		Vector3d p1 = TV.row(v1);
		Vector3d p2 = TV.row(v2);
		float curlen = (p2-p1).norm();
		f.segment<3>(3*v1) += (springList[i].stiffness*(curlen - springList[i].restLen)/curlen)*(p2-p1);
		f.segment<3>(3*v2) -= (springList[i].stiffness*(curlen - springList[i].restLen)/curlen)*(p2-p1);
	}
	for(int i=0; i<f.rows(); i++){
		i++;
		f(i) = gravity*vertex_masses(i/3);
		i++;
	}
	////////////////////////////////////////
	// cout<<endl<<"Force"<<endl;

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

		f(3*indices(0)) += F_tet.col(0)[0];
		f(3*indices(0)+1) += F_tet.col(0)[1];
		f(3*indices(0)+2) += F_tet.col(0)[2];

		f(3*indices(1)) += F_tet.col(1)[0];
		f(3*indices(1)+1) += F_tet.col(1)[1]; 
		f(3*indices(1)+2) += F_tet.col(1)[2];

		f(3*indices(2)) += F_tet.col(2)[0];
		f(3*indices(2)+1) += F_tet.col(2)[1];
		f(3*indices(2)+2) += F_tet.col(2)[2];

		f(3*indices(3)) += F_tet.col(3)[0];
		f(3*indices(3)+1) += F_tet.col(3)[1];
		f(3*indices(3)+2) += F_tet.col(3)[2];

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

void Simulation::ImplicitCalculateForces( MatrixXd& TVk, SparseMatrix<double>& forceGradient, VectorXd& v_k, VectorXd& f){
	// //gravity
	f.setZero();
	//UNCOMMENT
	// for(unsigned int i=0; i<M.tets.size(); i++){
	// 	double vertex_mass = M.tets[i].undeformedVol/4;//assume const density 1
	// 	Vector4i indices = M.tets[i].verticesIndex;

	// 	f(3*indices(0)+1) += vertex_mass*gravity;

	// 	f(3*indices(1)+1) += vertex_mass*gravity; 

	// 	f(3*indices(2)+1) += vertex_mass*gravity;

	// 	f(3*indices(3)+1) += vertex_mass*gravity;
	// }
	// //elastic
	// for(unsigned int i=0; i<M.tets.size(); i++){
	// 	Vector4i indices = M.tets[i].verticesIndex;
	// 	MatrixXd F_tet = M.tets[i].computeElasticForces(TVk, t%2);
	// 	f.segment<3>(3*indices(0)) += F_tet.col(0);
	// 	f.segment<3>(3*indices(1)) += F_tet.col(1);
	// 	f.segment<3>(3*indices(2)) += F_tet.col(2);
	// 	f.segment<3>(3*indices(3)) += F_tet.col(3);
	// }

	// //damping
	// f+= rayleighCoeff*forceGradient*v_k;

	////////////////////////////////////////
	for(int i=0; i<springList.size(); i++){
		int v1 = springList[i].v1;
		int v2 = springList[i].v2;
		Vector3d p1 = TVk.row(v1);
		Vector3d p2 = TVk.row(v2);
		float curlen = (p2-p1).norm();
		f.segment<3>(3*v1) += (springList[i].stiffness*(curlen - springList[i].restLen)/curlen)*(p2-p1);
		f.segment<3>(3*v2) -= (springList[i].stiffness*(curlen - springList[i].restLen)/curlen)*(p2-p1);

	}
	for(int i=0; i<f.rows(); i++){
		i++;
		f(i) += gravity*vertex_masses(i/3);
		i++;
	}
	////////////////////////////////////////
	return;
}

void Simulation::ImplicitCalculateElasticForceGradient(MatrixXd& TVk, SparseMatrix<double>& forceGradient){
	// vector<Trip> triplets1;
	// triplets1.reserve(CHANGE);	
	// for(unsigned int i=0; i<M.tets.size(); i++){
	// 	//Get P(dxn), dx = [1,0, 0...], then [0,1,0,....], and so on... for all 4 vert's x, y, z
	// 	//P is the compute Force Differentials blackbox fxn
	// 	MatrixXd local_gradForces(12, 12);
	// 	local_gradForces.setZero();

	// 	Vector12d dx(12);
	// 	dx.setZero();
	// 	Vector4i indices = M.tets[i].verticesIndex;
	// 	int kj;
	// 	for(unsigned int j=0; j<12; j++){
	// 		dx(j) = 1;
	// 		MatrixXd dForces = M.tets[i].computeForceDifferentials(TVk, dx);
	// 		kj = j%3;
	// 		//row in order for dfxi/dxi ..dfxi/dzl
	// 		triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[0], dForces(0,0)));
	// 		triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[0]+1, dForces(1,0)));
	// 		triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[0]+2, dForces(2,0)));

	// 		triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[1], dForces(0,1)));
	// 		triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[1]+1, dForces(1,1)));
	// 		triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[1]+2, dForces(2,1)));

	// 		triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[2], dForces(0,2)));
	// 		triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[2]+1, dForces(1,2)));
	// 		triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[2]+2, dForces(2,2)));

	// 		triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[3], dForces(0,3)));
	// 		triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[3]+1, dForces(1,3)));
	// 		triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[3]+2, dForces(2,3)));
	// 		dx(j) = 0; //ASK check is this efficient?
	// 	}
	// }
	// forceGradient.setFromTriplets(triplets1.begin(), triplets1.end());

	//OLD METHOD OF SETTING GLOBAL GRAD F, USING LOCAL GRAD F
	//UNCOMMENT
	// for(unsigned int i=0; i<M.tets.size(); i++){
	// 	//Get P(dxn), dx = [1,0, 0...], then [0,1,0,....], and so on... for all 4 vert's x, y, z
	// 	//P is the compute Force Differentials blackbox fxn
	// 	MatrixXd local_gradForces(12, 12);
	// 	local_gradForces.setZero();

	// 	Vector12d dx(12);
	// 	dx.setZero();
	// 	for(unsigned int j=0; j<12; j++){
	// 		dx(j) = 1;
	// 		MatrixXd dForces = M.tets[i].computeForceDifferentials(TVk, dx);
	// 		local_gradForces.col(j) = Map<VectorXd>(dForces.data(), dForces.cols()*dForces.rows());;// TODO values are really big. Check if correct
	// 		dx(j) = 0; //ASK check is this efficient?
	// 	}


	// 	//TODO optimize this somehow. Create Triplet list from local grad forces
		vector<Trip> triplets;
		triplets.reserve(3*vertices*3*vertices);//estimation of entries
	// 	Vector4i indices = M.tets[i].verticesIndex;
	// 	for(unsigned int ki=0; ki<4; ki++){
	// 		for(unsigned int kj=0; kj<4; kj++){
	// 			triplets.push_back(Trip(3*indices[ki], 3*indices[kj], local_gradForces(3*ki, 3*kj))); //dfxi/dxi
	// 			triplets.push_back(Trip(3*indices[ki], 3*indices[kj]+1, local_gradForces(3*ki, 3*kj+1))); //dfxi/dyi
	// 			triplets.push_back(Trip(3*indices[ki], 3*indices[kj]+2, local_gradForces(3*ki, 3*kj+2))); //dfxi/dzi

	// 			triplets.push_back(Trip(3*indices[ki]+1,  3*indices[kj], local_gradForces(3*ki+1, 3*kj))); //dfyi/dxi
	// 			triplets.push_back(Trip(3*indices[ki]+1,  3*indices[kj]+1, local_gradForces(3*ki+1, 3*kj+1))); //dfyi/dyi
	// 			triplets.push_back(Trip(3*indices[ki]+1,  3*indices[kj]+2, local_gradForces(3*ki+1, 3*kj+2))); //dfyi/dzi

	// 			triplets.push_back(Trip(3*indices[ki]+2,  3*indices[kj], local_gradForces(3*ki+2, 3*kj))); //dfzi/dxi
	// 			triplets.push_back(Trip(3*indices[ki]+2,  3*indices[kj]+1, local_gradForces(3*ki+2, 3*kj+1))); //dfzi/dyi
	// 			triplets.push_back(Trip(3*indices[ki]+2,  3*indices[kj]+2, local_gradForces(3*ki+2, 3*kj+2))); //dfzi/dzi
	// 		}
	// 	}
	// 	//Now insert triplets into global matrix
	// 	forceGradient.setFromTriplets(triplets.begin(), triplets.end());
	// }
	for(int i=0; i<springList.size(); i++){
		int v1 = springList[i].v1;
		int v2 = springList[i].v2;
		Vector3d p1 = TVk.row(v1);
		Vector3d p2 = TVk.row(v2);
		double curlen = (p2-p1).norm();
		
		//figure out placement into global dF matrix
		float K = springList[i].stiffness;
		float L = springList[i].restLen;

		Matrix3d I = MatrixXd::Identity(3,3);
		Matrix3d localdF = -K*(1-L/curlen)*I - K*L*(p2-p1)*(p2-p1).transpose()/curlen/curlen/curlen;
		for(int j=0; j<3; j++){
			for(int q=0; q<3; q++){
				triplets.push_back(Trip(3*v1+j, 3*v1+q, localdF.coeff(j,q)));
                triplets.push_back(Trip(3*v2+j, 3*v2+q, localdF.coeff(j,q)));
                triplets.push_back(Trip(3*v1+j, 3*v2+q, -localdF.coeff(j,q)));
                triplets.push_back(Trip(3*v2+j, 3*v1+q, -localdF.coeff(j,q)));
			}
		}
		
	}
	forceGradient.setFromTriplets(triplets.begin(), triplets.end());
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
	setXIntoTV(x_old);
	createForceVector();
	v_old = v_old + timestep*InvMass*f;

}
void Simulation::renderImplicit(){
	t+=1;
	//Implicit Code
	v_k.setZero();
	x_k.setZero();
	v_k = v_old;
	x_k = x_old;
	forceGradient.setZero();
	bool Nan=false;
	int NEWTON_MAX = 100, i =0;

	for( i=0; i<NEWTON_MAX; i++){
		grad_g.setZero();
	
		ImplicitXtoTV(x_k, TVk);//TVk value changed in function
		ImplicitCalculateElasticForceGradient(TVk, forceGradient); 
		ImplicitCalculateForces(TVk, forceGradient, v_k, f);

		VectorXd g = x_k - (x_old + timestep*v_k + timestep*timestep*InvMass*f);
		grad_g = Ident - timestep*timestep*InvMass*forceGradient;

		//solve for delta v
		// Conj Grad
		// ConjugateGradient<SparseMatrix<double>> cg;
		// cg.compute(grad_g);
		// VectorXd deltaX = -1*cg.solve(g);

		//Sparse Cholesky LL^T
		SimplicialLLT<SparseMatrix<double>> llt;
		llt.compute(grad_g);
		VectorXd deltaX = -1* llt.solve(g);


		// x_k = x_old+timestep*v_k;
		// v_k = v_k + deltaV;
		x_k+=deltaX;

		if(v_k != v_k){
			Nan = true;
			break;
		}
		if(g.squaredNorm()<0.00001){
			break;
		}
	}
	if(Nan || i== NEWTON_MAX){
		cout<<"ERROR: Newton's method doesn't converge"<<endl;
		cout<<i<<endl;
		exit(0);
	}

	// x_old = x_old+timestep*v_k;
	// v_old = v_k;
	v_old = (x_k - x_old)/timestep;
	x_old = x_k;

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
	double strainE = 0;
	//UNCOMMENT
	// for(unsigned int i=0; i<M.tets.size(); i++){
	// 	strainE += M.tets[i].undeformedVol*M.tets[i].energyDensity;		
	// }
	for(unsigned int i=0; i<springList.size(); i++){
		int v1 = springList[i].v1;
		int v2 = springList[i].v2;
		Vector3d p1 = TV.row(v1);
		Vector3d p2 = TV.row(v2);
		float curlen = (p2-p1).norm() - springList[i].restLen;
		strainE += 0.5*springList[i].stiffness*(curlen*curlen);
	}

	TotalEnergy+= gravityE + kineticE + strainE;
	// cout<<endl<<"Grav E"<<endl;
	// cout<<strainE<<endl;
	cout<<v_old.squaredNorm()<<endl;
	cout<<"Tot E"<<endl;
	cout<<TotalEnergy<<endl;
	cout<<"Strain E"<<endl;
	cout<<strainE<<endl;
	energyFile<<t<<", "<<TotalEnergy<<"\n";
	strainEnergyFile<<t<<", "<<strainE<<"\n";
	kineticEnergyFile<<t<<", "<<kineticE<<"\n";
	gravityEnergyFile<<t<<", "<<gravityE<<"\n";

	////////////////////////////////////
	
	////////////////////
	// if(momentumFile.is_open()){
	// 	double xp=0;
	// 	double yp=0;
	// 	double zp=0;
	// 	for(int i=0; i<v_old.rows(); ){
	// 		xp += v_old(i)*vertex_masses(i);
	// 		i++;
	// 		yp += v_old(i)*vertex_masses(i);
	// 		i++;
	// 		zp += v_old(i)*vertex_masses(i);
	// 		i++;
	// 	}
	// 	momentumFile<<t<<","<<xp<<","<<yp<<","<<zp<<"\n";
	// }else{
	// 	cout<<"no open file"<<endl;
	// }
	///////////////////
}


MatrixXi TT_One_G;
MatrixXd TV_One_G;

Simulation Sim;

bool drawLoopTest(igl::viewer::Viewer& viewer){
	viewer.data.clear();
	Sim.render();
	

	viewer.data.add_points(Sim.TV, RowVector3d(1,0,0));
	viewer.data.add_edges(Sim.TV.row(0), Sim.TV.row(1), RowVector3d(1,0,0));
	viewer.data.add_edges(Sim.TV.row(0), Sim.TV.row(2), RowVector3d(1,0,0));
	viewer.data.add_edges(Sim.TV.row(0), Sim.TV.row(3), RowVector3d(1,0,0));
	viewer.data.add_edges(Sim.TV.row(1), Sim.TV.row(2), RowVector3d(1,0,0));
	viewer.data.add_edges(Sim.TV.row(1), Sim.TV.row(3), RowVector3d(1,0,0));
	viewer.data.add_edges(Sim.TV.row(2), Sim.TV.row(3), RowVector3d(1,0,0));

	// viewer.data.add_edges(Sim.TV.row(4), Sim.TV.row(3), RowVector3d(0,1,0));
	// viewer.data.add_edges(Sim.TV.row(4), Sim.TV.row(0), RowVector3d(0,0,1));
	// viewer.data.add_edges(Sim.TV.row(4), Sim.TV.row(2), RowVector3d(0,0,0));

	// viewer.data.add_edges(Sim.TV.row(5), Sim.TV.row(3), RowVector3d(0,1,0));
	// viewer.data.add_edges(Sim.TV.row(5), Sim.TV.row(0), RowVector3d(0,0,1));
	// viewer.data.add_edges(Sim.TV.row(5), Sim.TV.row(1), RowVector3d(0,0,0));

	// viewer.data.add_edges(Sim.TV.row(6), Sim.TV.row(3), RowVector3d(0,1,0));
	// viewer.data.add_edges(Sim.TV.row(6), Sim.TV.row(0), RowVector3d(0,0,1));
	// viewer.data.add_edges(Sim.TV.row(6), Sim.TV.row(5), RowVector3d(0,0,0));


	viewer.data.add_edges(RowVector3d(0,0,0), RowVector3d(100,0,0), RowVector3d(1,1,1));
	viewer.data.add_edges(RowVector3d(0,0,0), RowVector3d(0,100,0), RowVector3d(1,1,1));
	//viewer.data.add_edges(RowVector3d(0,0,0), RowVector3d(0,0,100), RowVector3d(1,1,1));
	
	return false;
}

bool drawLoop(igl::viewer::Viewer& viewer){
	viewer.data.clear();
	cout<<Sim.t<<endl;
	Sim.render();

	for(unsigned int i=0; i<Sim.mapV2TV.size(); i++){
		V.row(i) = Sim.TV.row(Sim.mapV2TV[i]);
	}

	viewer.data.add_edges(RowVector3d(0,0,0), RowVector3d(100,0,0), RowVector3d(1,1,1));
	viewer.data.add_edges(RowVector3d(0,0,0), RowVector3d(0,100,0), RowVector3d(1,1,1));
	//viewer.data.add_edges(RowVector3d(0,0,0), RowVector3d(0,0,100), RowVector3d(1,1,1));
	
	viewer.data.set_mesh(V, F);
	viewer.data.set_face_based(false);
	return false;
}

void useFullObject(bool headless, float timestep, int iterations, char method){
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
	// Sim.fixVertices(1);
	if(!headless){
		igl::viewer::Viewer viewer;
		viewer.callback_pre_draw = &drawLoop;
		viewer.launch();
	}else{
		while(Sim.t<iterations){
			Sim.render();
			for(unsigned int i=0; i<Sim.mapV2TV.size(); i++){
				V.row(i) = Sim.TV.row(Sim.mapV2TV[i]);
			}
			if(Sim.t%1000 == 0){
				//igl::writeOBJ("../output2/object"+to_string(Sim.t)+".obj", V, F);
			}
		}
	}
}

void useMyObject(bool headless, float timestep, int iterations, char method){
	vector<int> mapV2TV;


	TT_One_G.resize(1, 4);
	TT_One_G<< 0, 1, 2, 3;

	TV_One_G.resize(4, 3);
	TV_One_G << 0, 0, 10, //affect
				10, 0, 0,
				0, 10, 0,
				0, 0, 0;


	// TT_One_G.resize(4, 4);
	// TT_One_G<< 0, 2, 1, 3,
	// 			4, 2, 0, 3,
	// 			5, 3, 0, 1,
	// 			6, 3, 0, 5;

	// TV_One_G.resize(7, 3);
	// TV_One_G << 0, 0, 10, //affect
	// 			0, 10, 0,
	// 			10, 0, 0,
	// 			0, 0, 0,
	// 			10, -10, 0,
	// 			-10, 10, 0,
	// 			-10, 0, 0;
				
	Sim.initializeSimulation(timestep, method, TT_One_G, TV_One_G, mapV2TV);

	igl::viewer::Viewer viewer;
	bool boolVariable = true;
	double timeVariable = 0.001;
	if(!headless){
		viewer.callback_init = [&](igl::viewer::Viewer& viewer)
		{
			viewer.ngui-> addGroup("My settings");

			viewer.ngui-> addVariable("bool", boolVariable);
			viewer.ngui-> addVariable("double", timeVariable);

			// viewer.ngui-> addButton("Print Hello", [](){cout<<"Hello"<<endl;});

			viewer.screen-> performLayout();
			return false;
		};
		viewer.callback_pre_draw = &drawLoopTest;
		viewer.launch();
	}else{
		while(Sim.t<iterations){
			Sim.render();
		}
	}
}

int main(int argc, char *argv[])
{
	float timestep = 0;
	char headless;
	char method;
	int iterations = 0;
	int object = 0;
	string line;
	ifstream configFile ("../config.txt");
	if(configFile.is_open()){
		getline(configFile, line);
		timestep = atof(line.c_str());
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
		rayleighCoeff = atof(line.c_str());
		cout<<rayleighCoeff<<endl;

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
