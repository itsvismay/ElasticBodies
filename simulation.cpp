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
			binarySearchYoungs(newMoveIndices, newTV, newTT, fixVertices.size(), B);
			// syntheticTests(newMoveIndices, newTV, newTT, fixVertices.size(), B);
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

}

void Simulation::staticSolveStep(double move_step, int ignorePastIndex, vector<int>& moveVertices, MatrixXd& TV,  MatrixXi& TT){
	//Move vertices slightly in x,y,z direction
	// [v, v, v..., f, f, ...(m), (m), (m)...]
	for(unsigned int i=0; i<moveVertices.size(); i++){
		TV.row(TV.rows()-i-1)[0]+= move_step;//move step
	}
	
	//Newtons method static solve for minimum Strain E
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
		
		if (fblock.squaredNorm()/fblock.size() < 0.00001){
			break;
		}
	}

	if(k== NEWTON_MAX){
		cout<<"ERROR Static Solve: Newton max reached"<<endl;
		cout<<k<<endl;
		exit(0);
	}				
}

void Simulation::binarySearchYoungs(vector<int> moveVertices, MatrixXd& TV, MatrixXi& TT, int fv, MatrixXd& B){
	cout<<"############Starting Binary Search for Youngs ######################"<<endl;

	ofstream distvLoadFile;
	distvLoadFile.open("../Scripts/distvLoad.txt");

	ofstream youngsFile;
	youngsFile.open("../Scripts/youngsnew.txt");

	//REAL VALUES FROM EXPERIMENT
	//dist, load
	// vector<pair<double, double>> realLoads = 
	// {
	// 	{0.100244,	101.073609},
	// 	{0.20048,	184.592875},
	// 	{0.300347,	265.06366},
	// 	{0.399884,	342.551726},
	// 	{0.499522,	417.940362},
	// 	{0.599679,	492.371951},
	// 	{0.700077,	565.750461},
	// 	{0.800373,	638.291505},
	// 	{0.900439,	709.507457},
	// 	{1.000194,	779.029706},
	// 	{1.099773,	846.406884},
	// 	{1.199565,	911.991954},
	// 	{1.299774,	975.713439},
	// 	{1.400103,	1037.549957},
	// 	{1.500424,	1097.35856},
	// 	{1.600455,	1154.854633},
	// 	{1.700085,	1209.47858},
	// 	{1.799595,	1260.718421},
	// 	{1.8996,	1308.14702},
	// 	{1.999928,	1351.454444},
	// 	{2.100283,	1389.728216},
	// 	{2.200424,	1421.50609},
	// 	{2.300369,	1444.923954},
	// 	{2.39994,	1458.153129},
	// 	{2.499569,	1460.703522},
	// 	{2.599659,	1453.403456},
	// 	{2.699964,	1438.137844},
	// 	{2.800305,	1427.789346},
	// 	{2.900472,	1412.155562},
	// 	{3.000338,	1393.765691}
	// };
	
	vector<pair<double, double>> realLoads = 
	{
		{0.06, 37.0006},
		{0.66, 405.748},
		{1.26, 772.218},
		{1.86, 1136.43},
		{2.46, 1498.42},
		{3.06, 1858.2},
		{3.66, 2215.81},
		{4.26, 2571.27},
		{4.86, 2924.6},
		{5.46, 3275.83},
		{6.06, 3625}

	};

	vector<double> derivedYoungs;

	//size of move
	double move_amount = 6;
	int number_of_moves = 100;
	double dist_moved = 0;
	double curr_youngs = 1;
	double load_scalar = 0;
	int ignorePastIndex = TV.rows() - moveVertices.size() - fv;

	int count=0;
	
	int reals_index =0;
	
	//New Binary Search
	while(dist_moved<move_amount){
		M.setNewYoungsPoissons(1000000, 0.35);

		//Newton Solve for positions
		while(reals_index<realLoads.size() && (dist_moved<realLoads[reals_index].first && (abs(dist_moved-realLoads[reals_index].first)>1e-5))){
			cout<<"	Move next step"<<endl;
			dist_moved += move_amount/number_of_moves;//move step
			cout<<"distance moved"<<dist_moved<<endl;
			cout<<(dist_moved<realLoads[reals_index].first && (abs(dist_moved-realLoads[reals_index].first)>1e-5) )<<endl;
			staticSolveStep(move_amount/number_of_moves, ignorePastIndex, moveVertices, TV, TT);	
			
		}

		//binary search for youngs
		double min_youngs = 800000;
		double max_youngs = 5000000;
		load_scalar = 0;
		while(abs(load_scalar-realLoads[reals_index].second)>(realLoads[reals_index].second/1000) && dist_moved<move_amount){
			load_scalar =0;
			curr_youngs = (min_youngs+max_youngs)/2; //just a guess
			M.setNewYoungsPoissons(curr_youngs, 0.35);
			VectorXd f;
			f.resize(3*TV.rows());
			calculateElasticForces(f, TV);


			//Calculate Load on moving verts
			Vector3d load(0,0,0);
			for(unsigned int i=f.size() - 3*moveVertices.size(); i<f.size(); i++){
				load+=f.segment<3>(i);
				i++;
				i++;
			}

			load_scalar = abs(load(0)/1000);

			count++;
			if((load_scalar - realLoads[reals_index].second)>0){
				cout<<"too high"<<endl;
				max_youngs = curr_youngs;
			}else{
				cout<<"too low"<<endl;
				min_youngs = curr_youngs;
			}
			cout<<"Binary Search "<< reals_index<<endl;
			cout<<"Calculated Load "<<load_scalar<<endl;
			cout<<"Actual Load "<<realLoads[reals_index].second<<endl;
			cout<<"min_youngs "<<min_youngs<<endl;
			cout<<"curr_youngs "<<curr_youngs<<endl;
			cout<<"max_youngs "<<max_youngs<<endl;
			cout<<"Solve Tet mu, lambda "<<M.tets[0].mu<<", "<<M.tets[0].lambda<<endl;
			cout<<"distance moved"<<dist_moved<<endl;
			cout<<"----------------"<<endl;
		}
		youngsFile<<dist_moved<<", "<<curr_youngs<<endl;
		reals_index+=1;
	}


	distvLoadFile.close();
	youngsFile.close();
	//system("( speaker-test -t sine -f 1000 )& pid=$! ; sleep 5s ; kill -9 $pid");
}

void Simulation::syntheticTests(vector<int> moveVertices, MatrixXd& TV, MatrixXi& TT, int fv, MatrixXd& B){
	cout<<"############Starting Synthetic Load Generation######################"<<endl;

	ofstream generateLoadsFile;
	generateLoadsFile.open("../Scripts/syntheticGeneratedLoads.txt");
	
	int setYoungs = 2e6;
	M.setNewYoungsPoissons(setYoungs, 0.35);

	double dist_moved = 0;
	double move_amount = 6;
	double number_of_moves = 100;
	double load_scalar =0;
	double number_of_data_points =10;

	int count = 1;

	while(dist_moved<move_amount){
		double move_step = move_amount/number_of_moves;
		int ignorePastIndex = TV.rows() - moveVertices.size() - fv;
		dist_moved += move_step;
		staticSolveStep(move_step, ignorePastIndex, moveVertices, TV, TT);

		if(count%10 == 0){
			cout<<"	print dist/load to file"<<endl;
			VectorXd f;
			f.resize(3*TV.rows());
			calculateElasticForces(f, TV);
			//Calculate Load on moving verts
			Vector3d load(0,0,0);
			for(unsigned int i=f.size() - 3*moveVertices.size(); i<f.size(); i++){
				load+=f.segment<3>(i);
				i++;
				i++;
			}

			load_scalar = abs(load(0)/1000);

			generateLoadsFile<<dist_moved<<", "<<load_scalar<<endl;
		}
		count++;
	}

	generateLoadsFile.close();
	cout<<"############End Synthetic Load Generation######################"<<endl;
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

