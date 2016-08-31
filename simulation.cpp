#include <iostream>
#include <vector>
#include <pthread.h>
#include <fstream>
#include <string>
#include <math.h>
#include <lbfgs.h>
#include <set>
#include <ctime>
#include <fstream>
#include <igl/writeOBJ.h>
#include <igl/barycenter.h>

#include "simulation.h"
#include "globals.h"



using namespace Eigen;
using namespace std;
typedef Eigen::Triplet<double> Trip;
typedef Matrix<double, 12, 1> Vector12d;

Simulation::Simulation(void){}

int Simulation::initializeSimulation(double deltaT, int iterations, char method, MatrixXi& TT, MatrixXd& TV, MatrixXd& B, vector<int>& moveVertices, vector<int> fixVertices, double youngs, double poissons){
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

	if(moveVertices.size()>0 or fixVertices.size()>0){
		MatrixXd newTV;
		newTV.resize(TV.rows(), TV.cols());
		newTV.setZero();
		MatrixXi newTT;
		newTT.resize(TT.rows(), TT.cols());
		newTT.setZero();

		//TODO: Make this shit more efficient
		//Hash maps or something
		vector<int> vertexNewIndices;
		for(int i=0; i<TV.rows(); i++){
			bool flag =false;
			for(unsigned int j=0; j<fixVertices.size(); j++){
				if(i==fixVertices[j]){
					flag = true;
				}
			}
			for(unsigned int j=0; j<moveVertices.size(); j++){
				if(i==moveVertices[j]){
					flag = true;
				}
			}
			// if vertex not fixed or moved, re-index to front
			//[v, v, v, v...., f, f, f...., m, m, m...,m]
			if(!flag){
				vertexNewIndices.push_back(i);
			}
		}
		//re-index fixed verts
		for(unsigned int j=0; j<fixVertices.size(); j++){
			vertexNewIndices.push_back(fixVertices[j]);
		}
		//re-index move verts
		for(unsigned int j=0; j<moveVertices.size(); j++){
			vertexNewIndices.push_back(moveVertices[j]);
		}

		//these are the new indices for the fixed verts
		vector<int> newfixIndices;
		for(unsigned int i= vertexNewIndices.size() - (moveVertices.size() + fixVertices.size()); i<(vertexNewIndices.size()-moveVertices.size()); i++){
			newfixIndices.push_back(i);
		}

		//new indices for the moving verts
		vector<int> newMoveIndices;
		for(unsigned int i= vertexNewIndices.size() - moveVertices.size(); i<vertexNewIndices.size(); i++){
			newMoveIndices.push_back(i);
		}


		reIndexTVandTT(vertexNewIndices, fixVertices.size(), moveVertices.size(), TV, TT, newTV, newTT);


		igl::barycenter(newTV, newTT, B);
		//Initialize Solid Mesh
		M.initializeMesh(newTT, newTV, youngs, poissons);
		if(moveVertices.size() != 0){
			setInitPosition(newMoveIndices, newTV, newTT, fixVertices.size(), B);
		}
		
		integrator->initializeIntegrator(deltaT, M, newTV, newTT);
		integrator->fixVertices(newfixIndices);

	}else{
		igl::barycenter(TV, TT, B);
		M.initializeMesh(TT, TV, youngs, poissons);
		integrator->initializeIntegrator(deltaT, M, TV, TT);
		integrator->fixVertices(fixVertices);
	}

	
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

//TODO: Clean up function params size Fixed and size Move are not needed
void Simulation::reIndexTVandTT(vector<int> newVertsIndices, int sizeFixed, int sizeMove, MatrixXd& TV, MatrixXi& TT, MatrixXd& newTV, MatrixXi& newTT){
	//apply re-index to TV
	for(unsigned int i=0; i<newVertsIndices.size(); i++){
		newTV.row(i) = TV.row(newVertsIndices[i]);
	}

	//create map out of newVertsIndex
	//map keys = newVertsIndex values = old indices in TV
	//map vals = newVertsIndex index = new indices in TV
	map<int, int> oldToNewmap;
	pair<map<int, int>::iterator, bool> err;
	for(unsigned int i=0; i<newVertsIndices.size(); i++){
		err = oldToNewmap.insert(pair<int, int>(newVertsIndices[i], i));
		if(err.second==false){
			cout<<"ERROR::Simulation.cpp::reIndexTVandTT::>>Map already contains this value(";
			cout<< err.first->second <<". Indices should not be repeated"<<endl;
		}
	}

	//Check map, see if its working
	// map<int,int>::iterator it = oldToNewmap.begin();
	// for (it=oldToNewmap.begin(); it!=oldToNewmap.end(); ++it)
	// 	cout << it->first << " => " << it->second << '\n';


	//apply re-index to TT
	for(int i=0; i< TT.rows(); i++){
		for(int j=0; j<4; j++){
			newTT.row(i)[j] = oldToNewmap.find(TT.row(i)[j])->second;
		}
	}
}

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
	triplets1.reserve(12*12*M.tets.size());	
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
			dx(j) = 0;
		}
	}
	forceGradient.setFromTriplets(triplets1.begin(), triplets1.end());
	return;
}

void Simulation::setInitPosition(vector<int> moveVertices, MatrixXd& TV, MatrixXi& TT, int fv, MatrixXd& B){	
	MatrixXd InitialTV = TV;//Use this to reset the vertices every iteration of binary search
							//Use until some type of optimization method is implemented.

	ofstream distvLoadFile;
	distvLoadFile.open("../Scripts/distvLoad.txt");
	ofstream youngsFile;
	youngsFile.open("../Scripts/youngs.txt");

	//REAL VALUES FROM EXPERIMENT
	//dist, load
	vector<pair<double, double>> realLoads = 
	{
		{0.10058, 137.37688},
		{0.20064, 226.74218},
		{0.30028, 312.74922},
		{0.39975, 395.96558},
		{0.49947, 476.31368},
		{0.59996, 556.20089},
		{0.70037, 635.27435},
		{0.80035, 712.96558},
		{0.9005, 788.87},
		{1.00022, 863.39355},
		{1.0997, 935.25},
		{1.19933, 1004.94318},
		{1.29982, 1072.37},
		{1.40039, 1137.45287},
		{1.50063, 1200.1016},
		{1.60052, 1259.40174},
		{1.6999, 1314.902},
		{1.79946, 1366.0726},
		{1.89978, 1412.04},
		{2.00027, 1452.18},
		{2.10042, 1483.09},
		{2.2004, 1500.43219},
		{2.30046, 1501.043},
		{2.39993, 1486.33704},
		{2.49931, 1460.93},
		{2.59955, 1434.53971}
	};
	vector<double> derivedYoungs;

	//size of move
	double move_amount = 2.6;
	int number_of_moves = 100;
	double dist_moved = 0;
	double curr_youngs = 1;

	int count=0;

	for(int j=0; j<realLoads.size(); j++){
		dist_moved = 0;
		double min_youngs = 600000;
		double max_youngs = 4000000;
		double load_scalar = 0;

		//Binary Search Code below
		while(abs(load_scalar-realLoads[j].second)>(realLoads[j].second/10)){
			curr_youngs = (min_youngs+max_youngs)/2; //just a guess
			M.setNewYoungsPoissons(curr_youngs, 0.35);
			dist_moved =0;
			load_scalar =0;
			TV = InitialTV;

			while(dist_moved< realLoads[j].first){
				//Move vertices slightly in x,y,z direction
				// [v, v, v..., f, f, ...(m), (m), (m)...]
				for(unsigned int i=0; i<moveVertices.size(); i++){
					TV.row(TV.rows()-i-1)[0]+= move_amount/number_of_moves;
				}
				dist_moved += move_amount/number_of_moves;
				
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
				// cout<<TV<<endl;
				cout<<"	Move-dist "<<dist_moved<<"--"<<endl;
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
					cout<<"		Newton Iter "<<k<<endl;
					if(x != x){
						cout<<"NAN"<<endl;
						exit(0);
					}
					// cout<<"fblock"<<endl;
					// cout<<fblock.squaredNorm()<<endl;
					if (fblock.squaredNorm()/fblock.size() < 0.00001){
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

				//Calculate Load on moving verts
				Vector3d load(0,0,0);
				for(unsigned int i=f.size() - 3*moveVertices.size(); i<f.size(); i++){
					load+=f.segment<3>(i);
					i++;
					i++;
				}

				load_scalar = abs(load(0)/1000);
				//WRITE TO distvLoad FILE
				// distvLoadFile<<j*move_amount/number_of_moves<<", "<<abs(load(0)/1000)<<endl; //UNITS: Divide by 1000 for plotting purposes. Real data measured in N, I use milimeters for lengths.
				

			}
			cout<<"Binary Search "<< j<<endl;
			cout<<"Calculated Load "<<load_scalar<<endl;
			cout<<"Actual Load "<<realLoads[j].second<<endl;
			cout<<"min_youngs "<<min_youngs<<endl;
			cout<<"curr_youngs "<<curr_youngs<<endl;
			cout<<"max_youngs "<<max_youngs<<endl;
			cout<<"Solve Tet mu, lambda "<<M.tets[0].mu<<", "<<M.tets[0].lambda<<endl<<endl;
			cout<<"----------------"<<endl;
			//PRINT OBJ EACH STEP
			// printObj(count, TV, TT);
			// printObj(count, TV, TT, B);

			count++;

			if((load_scalar - realLoads[j].second)>0){
				max_youngs = curr_youngs;
			}else{
				min_youngs = curr_youngs;
			}
		}
		derivedYoungs.push_back(curr_youngs);
		cout<<endl<<endl;
		youngsFile<<dist_moved<<", "<<curr_youngs<<endl;
		system("( speaker-test -t sine -f 1000 )& pid=$! ; sleep 0.1s ; kill -9 $pid");
	}

	distvLoadFile.close();
	system("( speaker-test -t sine -f 1000 )& pid=$! ; sleep 5s ; kill -9 $pid");
}

void Simulation::printObj(int numberOfPrints, MatrixXd& TV, MatrixXi& TT, MatrixXd& B){

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
		V_temp.row(i*4+0) = TV.row(TT(s[i],0));
		V_temp.row(i*4+1) = TV.row(TT(s[i],1));
		V_temp.row(i*4+2) = TV.row(TT(s[i],2));
		V_temp.row(i*4+3) = TV.row(TT(s[i],3));
		F_temp.row(i*4+0) << (i*4)+0, (i*4)+1, (i*4)+3;
		F_temp.row(i*4+1) << (i*4)+0, (i*4)+2, (i*4)+1;
		F_temp.row(i*4+2) << (i*4)+3, (i*4)+2, (i*4)+0;
		F_temp.row(i*4+3) << (i*4)+1, (i*4)+2, (i*4)+3;
	}

	int q = system("mkdir -p ../TestsResults/NMTests/");
	igl::writeOBJ("../TestsResults/LoadTests/" + to_string(numberOfPrints)+".obj", V_temp, F_temp);

	return;
}
