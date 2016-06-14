#include <iostream>
#include <vector>
#include <pthread.h>
#include <fstream>
#include <string>
#include <math.h>
#include <lbfgs.h>
#include <set>
#include <ctime>
#include "Eigen/SPQRSupport"
#include <Eigen/CholmodSupport>

#include "simulation.h"
#include "globals.h"



using namespace Eigen;
using namespace std;
typedef Eigen::Triplet<double> Trip;
typedef Matrix<double, 12, 1> Vector12d;

Simulation::Simulation(void){}

int Simulation::initializeSimulation(double deltaT, int iterations, char method, MatrixXi& TT, MatrixXd& TV, MatrixXd& B, vector<int>& moveVertices, vector<int> fixVertices){
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

	//TODO: Make this shit more efficient
	//Hash maps or something
	vector<int> vertexNewIndices;
	for(int i=0; i<TV.rows(); i++){
		bool flag =false;
		for(int j=0; j<fixVertices.size(); j++){
			if(i==fixVertices[j]){
				flag = true;
			}
		}
		for(int j=0; j<moveVertices.size(); j++){
			if(i==moveVertices[j]){
				flag = true;
			}
		}
		if(!flag){
			vertexNewIndices.push_back(i);
		}
	}
	for(int j=0; j<fixVertices.size(); j++){
		vertexNewIndices.push_back(fixVertices[j]);
	}
	for(int j=0; j<moveVertices.size(); j++){
		vertexNewIndices.push_back(moveVertices[j]);
	}

	vector<int> newfixIndices;
	for(int i= vertexNewIndices.size() - (moveVertices.size() + fixVertices.size()); i<(vertexNewIndices.size()-moveVertices.size()); i++){
		newfixIndices.push_back(i);
	}
	vector<int> newMoveIndices;
	for(int i= vertexNewIndices.size() - moveVertices.size(); i<vertexNewIndices.size(); i++){
		newMoveIndices.push_back(i);
	}
	for(int i=0; i<vertexNewIndices.size(); i++){
		// cout<<vertexNewIndices[i]<<endl;
	}

	// cout<<"Before"<<endl;
	// cout<<TV<<endl;
	// cout<<TT<<endl;

	// reIndexTV(vertexNewIndices, TV, TT);
	// cout<<"After"<<endl;
	// cout<<TV<<endl;
	// cout<<TT<<endl;
	//Initialize Solid Mesh
	M.initializeMesh(TT, TV);
	
	// setInitPosition(moveVertices, TV, TT, fixVertices.size());

	integrator->initializeIntegrator(deltaT, M, TV, TT);
	integrator->fixVertices(fixVertices);
	return 1;
}

void Simulation::headless(){
	clock_t begin = clock();

	while(integrator->simTime<iters){
		integrator->render();
	}

	clock_t end = clock();
	cout<<"Seconds Elapsed: "<<double(end-begin)/CLOCKS_PER_SEC<<endl;
}
void Simulation::render(){
	integrator->render();
}

void Simulation::reIndexTV(vector<int> newVertsIndices, MatrixXd& TV, MatrixXi& TT){
	MatrixXd newTV = TV;
	for(int i=0; i<newVertsIndices.size(); i++){
		newTV.row(i) = TV.row(newVertsIndices[i]);
	}
	TV = newTV;
	for(int k=0; k<TT.rows(); k++){
		for(int j=0; j<4; j++){
			int ind = TT.row(k)[j];
			TT.row(k)[j] = newVertsIndices[ind];
		}
	}
}

// int Simulation::reIndexClampedVertices(vector<int>& moveVertices, MatrixXd& TV, MatrixXi& TT){
// 	//Re-index clamped vertices
// 	for(int i=0; i<moveVertices.size(); i++){
// 		//re-index vertices
// 		//take TV row at index moveVertices(i)
// 		//replace with TV row at index TV.rows - (moveVertices.size() -i)
// 		//Vector3d ro1 = TV.row(moveVertices(i));
// 		Vector3d ro2 = TV.row(TV.rows() - moveVertices.size()+ i);
// 		TV.row(TV.rows() - moveVertices.size()+ i) = TV.row(moveVertices[i]);
// 		TV.row(moveVertices[i]) = ro2;

// 		//re-index all the tet pointers
// 		int v1 = moveVertices[i];
// 		int v2 = TV.rows() - moveVertices.size() +i;
// 		for(int k=0; k< TT.rows(); k++){
// 			for(int j=0; j<4; j++){
// 				if (TT.row(k)[j]== v1){
// 					TT.row(k)[j] = v2;
// 				}else if(TT.row(k)[j]==v2){
// 					TT.row(k)[j] = v1;
// 				}
// 			}
// 		}
// 	}

// 	return 1;
// }

void Simulation::setTVtoX(VectorXd &x, MatrixXd &TV){
	//TV to X
	for(unsigned int i = 0; i < M.tets.size(); i++){
		Vector4i indices = M.tets[i].verticesIndex;

		x(3*indices(0)) = TV.row(indices(0))[0];
		x(3*indices(0)+1) = TV.row(indices(0))[1];
		x(3*indices(0)+2) = TV.row(indices(0))[2];

		x(3*indices(1)) = TV.row(indices(1))[0];
		x(3*indices(1)+1) = TV.row(indices(1))[1];
		x(3*indices(1)+2) = TV.row(indices(1))[2];

		x(3*indices(2)) = TV.row(indices(2))[0];
		x(3*indices(2)+1) = TV.row(indices(2))[1];
		x(3*indices(2)+2) = TV.row(indices(2))[2];

		x(3*indices(3)) = TV.row(indices(3))[0];
		x(3*indices(3)+1) = TV.row(indices(3))[1];
		x(3*indices(3)+2) = TV.row(indices(3))[2];
	}
	return;
}

void Simulation::calculateElasticForces(VectorXd &f, MatrixXd& TV){
	f.setZero();
	//elastic
	for(unsigned int i=0; i< M.tets.size(); i++){
		Vector4i indices = M.tets[i].verticesIndex;
		MatrixXd F_tet = M.tets[i].computeElasticForces(TV, 1);
		f.segment<3>(3*indices(0)) += F_tet.col(0);
		f.segment<3>(3*indices(1)) += F_tet.col(1);
		f.segment<3>(3*indices(2)) += F_tet.col(2);
		f.segment<3>(3*indices(3)) += F_tet.col(3);
	}
	return;
}

void Simulation::calculateForceGradient(MatrixXd &TVk, SparseMatrix<double>& forceGradient){
	forceGradient.setZero();
	
	vector<Trip> triplets1;
	triplets1.reserve(3*TVk.rows()*3*TVk.rows());	
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

void Simulation::setInitPosition(vector<int> moveVertices, MatrixXd& TV, MatrixXi& TT, int fv){	
	//size of move
	double move_amount = 1;
	int number_of_moves = 10;
	for(int j=0; j<number_of_moves; j++){
		//Move vertices slightly
		for(int i=0; i<moveVertices.size(); i++){
			TV.row(TV.rows()-i-1)[1]+= move_amount/number_of_moves;
		}
		
		//Newtons method static solve for minimum Strain E
		int ignorePastIndex = TV.rows() - moveVertices.size() - fv;
		double strainE;
		SparseMatrix<double> forceGradient;
		forceGradient.resize(3*TV.rows(), 3*TV.rows());
		SparseMatrix<double> forceGradientStaticBlock;
		forceGradientStaticBlock.resize(3*ignorePastIndex, 3*ignorePastIndex);
		VectorXd f, x;
		f.resize(3*TV.rows());
		f.setZero();
		x.resize(3*TV.rows());
		x.setZero();
		setTVtoX(x, TV);
		cout<<TV<<endl;
		cout<<(ignorePastIndex)<<endl;
		cout<<"-----"<<endl;
		// exit(0);
		int NEWTON_MAX = 100, k=0;
		for(k=0; k<NEWTON_MAX; k++){
			// X to TV
			TV.setZero();
			for(unsigned int i=0; i < M.tets.size(); i++){
				Vector4i indices = M.tets[i].verticesIndex;
				TV.row(indices(0)) = Vector3d(x(3*indices(0)), x(3*indices(0)+1), x(3*indices(0) +2));
				TV.row(indices(1)) = Vector3d(x(3*indices(1)), x(3*indices(1)+1), x(3*indices(1) +2));
				TV.row(indices(2)) = Vector3d(x(3*indices(2)), x(3*indices(2)+1), x(3*indices(2) +2));
				TV.row(indices(3)) = Vector3d(x(3*indices(3)), x(3*indices(3)+1), x(3*indices(3) +2)); 
			}
			calculateForceGradient(TV, forceGradient);
			calculateElasticForces(f, TV);
			
			//Block forceGrad and f to exclude the fixed verts
			forceGradientStaticBlock = forceGradient.block(0,0, 3*(ignorePastIndex), 3*ignorePastIndex);
			// cout<<"Force Gradient"<<endl;
			// cout<<forceGradient<<endl<<endl;
			// cout<<"FG Block"<<endl;
			// cout<<forceGradientStaticBlock<<endlIn
			VectorXd fblock = f.head(ignorePastIndex*3);
			//Sparse QR 
			SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> sqr;
			sqr.compute(forceGradientStaticBlock);
			VectorXd deltaX = -1*sqr.solve(fblock);

			x.segment(0,3*(ignorePastIndex))+=deltaX;
			cout<<"X"<<k<<endl;
			if(x != x){
				cout<<"NAN"<<endl;
				exit(0);
			}
			if (fblock.squaredNorm() < 0.00001){
				break;
			}
		}
		if(k== NEWTON_MAX){
			cout<<"ERROR Static Solve: Newton max reached"<<endl;
			cout<<k<<endl;
			exit(0);
		}
		for(unsigned int i=0; i<M.tets.size(); i++){
			strainE += M.tets[i].undeformedVol*M.tets[i].energyDensity;		
		}
		cout<<"Forces"<<endl;
		cout<<f<<endl;
		cout<<"Strain E"<<endl;
		cout<<strainE<<endl;
		cout<<"New pos"<<endl;
		cout<<TV<<endl;
	}
	
}
