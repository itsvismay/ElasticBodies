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
double gravity = -9.8;

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
	for(unsigned int i=0; i<M.tets.size(); i++){
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

	////////////////////////////////////////
	for(unsigned int i=0; i<springList.size(); i++){
		int v1 = springList[i].v1;
		int v2 = springList[i].v2;
		Vector3d p1 = TV.row(v1);
		Vector3d p2 = TV.row(v2);
		float curlen = (p2-p1).norm();
		f.segment<3>(3*v1) += (springList[i].stiffness*(curlen - springList[i].restLen)/curlen)*(p2-p1);
		f.segment<3>(3*v2) -= (springList[i].stiffness*(curlen - springList[i].restLen)/curlen)*(p2-p1);
	}
	for(int i=0; i<vertices; i++){
		f(3*i+1) += gravity*vertex_masses(3*i);
	}
	////////////////////////////////////////
	// cout<<endl<<"Force"<<endl;

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

	for(unsigned int i=0; i<springList.size(); i++){
		int v1 = springList[i].v1;
		int v2 = springList[i].v2;
		Vector3d p1 = TVk.row(v1);
		Vector3d p2 = TVk.row(v2);
		float curlen = (p2-p1).norm();
		f.segment<3>(3*v1) += (springList[i].stiffness*(curlen - springList[i].restLen)/curlen)*(p2-p1);
		f.segment<3>(3*v2) -= (springList[i].stiffness*(curlen - springList[i].restLen)/curlen)*(p2-p1);

	}
	for(int i=0; i<vertices; i++){
		f(3*i+1) += gravity*vertex_masses(3*i);
	}
	////////////////////////////////////////
	return;
}

void Simulation::ImplicitCalculateElasticForceGradient(MatrixXd& TVk, SparseMatrix<double>& forceGradient){


	// 	//TODO optimize this somehow. Create Triplet list from local grad forces
	vector<Trip> triplets;
	triplets.reserve(3*vertices*3*vertices);//estimation of entries

	for(unsigned int i=0; i<springList.size(); i++){
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

		// solve for delta v
		// Conj Grad
		// ConjugateGradient<SparseMatrix<double>> cg;
		// cg.compute(grad_g);
		// VectorXd deltaX = -1*cg.solve(g);

		//Sparse Cholesky LL^T
		// SimplicialLLT<SparseMatrix<double>> llt;
		// llt.compute(grad_g);
		// VectorXd deltaX = -1* llt.solve(g);

		//Sparse QR 
		SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> sqr;
		sqr.compute(grad_g);
		VectorXd deltaX = -1*sqr.solve(g);

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
	for(int i=0; i<3*vertices; i++){
		if(!isFixed(i)){
			int k=3*i;
			// gravityE +=  massVector(k)*-1*gravity*(x_old(k));
			kineticE += 0.5*vertex_masses(k)*v_old(k)*v_old(k);
			
			k++;
			gravityE +=  vertex_masses(k)*-1*gravity*(x_old(k));
			kineticE += 0.5*vertex_masses(k)*v_old(k)*v_old(k);
			
			k++;
			// gravityE +=  massVector(k)*-1*gravity*(x_old(k));
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
	if(momentumFile.is_open()){
		double xp=0;
		double yp=0;
		double zp=0;
		for(int i=0; i<v_old.rows(); ){
			xp += v_old(i)*vertex_masses(i);
			i++;
			yp += v_old(i)*vertex_masses(i);
			i++;
			zp += v_old(i)*vertex_masses(i);
			i++;
		}
		momentumFile<<t<<","<<xp<<","<<yp<<","<<zp<<"\n";
	}else{
		cout<<"no open file"<<endl;
	}
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
	

	viewer.data.add_edges(RowVector3d(0,0,0), RowVector3d(100,0,0), RowVector3d(1,1,1));
	viewer.data.add_edges(RowVector3d(0,0,0), RowVector3d(0,100,0), RowVector3d(1,1,1));
	
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
	igl::copyleft::tetgen::tetrahedralize(V,F,"-pq1.414Y", TV,TT,TF);
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
	Sim.fixVertices(1);
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


	// TT_One_G.resize(2, 4);
	// TT_One_G<<  0, 1, 2, 3,
	// 			4, 0, 2, 3;

	// TV_One_G.resize(5, 3);
	// TV_One_G << 10, 0, 0, //affect
	// 			0, 10, 0,
	// 			0, 0, 10,
	// 			0, 0, 0,
	// 			0, -10, 0;
				
	Sim.initializeSimulation(timestep, method, TT_One_G, TV_One_G, mapV2TV);
	Sim.fixVertices(1);
	// Sim.fixVertices(2);
	// Sim.fixVertices(3);
	igl::viewer::Viewer viewer;
	if(!headless){
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