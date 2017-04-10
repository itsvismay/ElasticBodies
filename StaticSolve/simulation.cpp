#include "simulation.h"


Simulation::Simulation(void){}

int staticSolveDirection = 0;

int Simulation::initializeSimulation(double deltaT, int iterations, char method, MatrixXi& TT, MatrixXd& TV, MatrixXd& B, vector<int>& moveVertices, vector<int> fixVertices, double youngs, double poissons){
	iters = iterations;
	if (method =='e'){
		//integrator = new Verlet();
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
	//*********BEAM******************
	// move vertices
	for(int i=0; i<TV.rows(); i++){
	 	if(TV.row(i)[0]>=160){
	 		moveVertices.push_back(i);
	 	}
	}

	//fix vertices
	for(int i=0; i<TV.rows(); i++){
		if(TV.row(i)[0]<=30){
			fixVertices.push_back(i);
		}
	}
	//***************************
	VectorXd force;
	force.resize(3*TV.rows());
	force.setZero();
	TV_k = TV;
	// setInitPosition(force, fixVertices, moveVertices);
	cout<<"Fixed and moving"<<endl;
	cout<<fixVertices.size()<<endl;
	cout<<moveVertices.size()<<endl;
//

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

		VectorXd new_force;
		new_force.resize(3*TV.rows());
		reIndexTVandTT(vertexNewIndices, fixVertices.size(), moveVertices.size(), TV, TT, force, newTV, newTT, new_force);

		igl::barycenter(newTV, newTT, B);
		//Initialize Solid Mesh
		M.initializeMesh(newTT, newTV, youngs, poissons);
		if(moveVertices.size() != 0){
			// cout<<"Move vertices "<<moveVertices.size()<<endl;
			// cout<<"fix verts "<<fixVertices.size()<<endl;
			binarySearchYoungs(newMoveIndices, newTV, newTT, fixVertices.size(), B);
			// syntheticTests(newMoveIndices, newTV, newTT, fixVertices.size(), B);
			// staticSolveInitialPosition(newMoveIndices, newTV, newTT, fixVertices.size(), B);
			exit(0);

		}

		integrator->initializeIntegrator(deltaT, M, newTV, newTT);
		this->external_force = new_force;
		integrator->fixVertices(newfixIndices);

		// staticSolveNewtonsForces(newTV, newTT, B, new_force, ignorePastIndex);



	}else{
		igl::barycenter(TV, TT, B);
		M.initializeMesh(TT, TV, youngs, poissons);
		integrator->initializeIntegrator(deltaT, M, TV, TT);
		this->external_force = force;
		integrator->fixVertices(fixVertices);
	}


	return 1;
}

void Simulation::applyExternalForces(){
	this->external_force.setZero();
}

void Simulation::headless(){
	clock_t begin = clock();

	while(integrator->simTime<iters){
		integrator->render(this->external_force);
		cout<<"Min Displacement (called maxDisp in code)"<<endl;
		double disp =0;
		for(int i=0; i<this->putForceOnTheseVerts.rows(); i++){
			if (integrator->TV.row(this->putForceOnTheseVerts(i))(2) < disp)
				disp = integrator->TV.row(this->putForceOnTheseVerts(i))(2);
		}
		if(disp < maxDisp){
			maxDisp = disp;
		}
		cout<<maxDisp<<"\n";
	}
	optimizationFile<<maxDisp<<endl;

	clock_t end = clock();
	cout<<"TIME ELAPSED"<<endl<<endl;
	cout<<"Seconds Elapsed: "<<double(end-begin)/CLOCKS_PER_SEC<<endl;

}

void Simulation::render(){
	integrator->render(this->external_force);
	cout<<"Max Displacement (called maxDisp in code)"<<endl;
	double disp =0;
	for(int i=0; i<this->putForceOnTheseVerts.rows(); i++){
		if (integrator->TV.row(this->putForceOnTheseVerts(i))(2) < disp)
			disp = integrator->TV.row(this->putForceOnTheseVerts(i))(2);
	}
	if(disp < maxDisp){
		maxDisp = disp;
	}
	cout<<maxDisp<<endl;
}

//TODO: Clean up function params size Fixed and size Move are not needed
void Simulation::reIndexTVandTT(
	vector<int> newVertsIndices,
	int sizeFixed,
	int sizeMove,
	MatrixXd& TV,
	MatrixXi& TT,
	VectorXd& force,
	MatrixXd& newTV,
	MatrixXi& newTT,
	VectorXd& new_force){
	//apply re-index to TV
	for(unsigned int i=0; i<newVertsIndices.size(); i++){
		newTV.row(i) = TV.row(newVertsIndices[i]);
		new_force.segment<3>(3*i) = force.segment<3>(3*newVertsIndices[i]);
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

void Simulation::staticSolveNewtonsForces(MatrixXd& TV, MatrixXi& TT, MatrixXd& B, VectorXd& fixed_forces, int ignorePastIndex){
	cout<<"------------I am here----------------"<<endl;
	cout<<ignorePastIndex<<endl;
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

	int NEWTON_MAX = 10, k=0;
	for(k=0; k<NEWTON_MAX; k++){
		xToTV(x, TV);

		calculateForceGradient(TV, forceGradient);
		calculateElasticForces(f, TV);
		//PLAYGROUND - Check forces in mathematica
		cout<<TV<<endl;
		cout<<f<<endl;
		//--------------
		for(int i=0; i<fixed_forces.rows(); i++){
			if(fabs(fixed_forces(i))>0.00001){
				f(i) = fixed_forces(i);
				//cout<<f(i)<<endl;
			}
		}

		//Block forceGrad and f to exclude the fixed verts
		forceGradientStaticBlock = forceGradient.block(0,0, 3*(ignorePastIndex), 3*ignorePastIndex);
		VectorXd fblock = f.head(ignorePastIndex*3);

		// Conj Grad
		ConjugateGradient<SparseMatrix<double>> cg;
		cg.compute(forceGradientStaticBlock);
		if(cg.info() == Eigen::NumericalIssue){
			cout<<"ConjugateGradient numerical issue"<<endl;
			exit(0);
		}
		VectorXd deltaX = -1*cg.solve(fblock);

		x.segment(0,3*(ignorePastIndex))+=deltaX;
		cout<<"		Newton Iter "<<k<<endl;

		if(x != x){
			cout<<"NAN"<<endl;
			exit(0);
		}
		cout<<"fblock"<<endl;
		cout<<fblock.squaredNorm()/fblock.size()<<endl;

		if (fblock.squaredNorm()/fblock.size() < 0.00001){
			break;
		}

	}
	if(k== NEWTON_MAX){
		cout<<"ERROR Static Solve: Newton max reached"<<endl;
		cout<<k<<endl;
		exit(0);
	}
	double strainE = 0;
	for(int i=0; i< M.tets.size(); i++){
		strainE += M.tets[i].energy;
	}
	cout<<"strain E"<<strainE<<endl;
	cout<<"x[0] "<<x(0)<<endl;
	cout<<"x[1] "<<x(1)<<endl;


	cout<<"----------------------"<<endl;
	cout<<TV<<endl;
	exit(0);
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
		M.tets[i].computeElasticForces(TV, f);

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

void Simulation::setInitPosition(VectorXd& force, vector<int>& fixVertices, vector<int>& moveVertices){
	//TODO: implement this later - with Zack's code
	//hard coded the force file for now
	vector<int> temp;
	//cout<<force.rows()<<endl;
	ifstream forceInputFile (TUTORIAL_SHARED_PATH "shared/"+objectName+".txt");
	cout<<TUTORIAL_SHARED_PATH "shared/"+objectName+".txt"<<endl;
	if(forceInputFile.is_open()){
		string line;
		int index =0;
		int fixedIndex=0;
		while(getline(forceInputFile, line)){
			istringstream iss(line);
			double fx, fy, fz;
			int fixedOrNot; //1 is fixed, 0 not fixed
			if(!(iss >> fx >> fy >> fz >> fixedOrNot)){break;}
			if(fabs(fx + fy + fz)>0){
				moveVertices.push_back(index);
				temp.push_back(index - fixedIndex);
				// cout<<index<<endl;
			}
			force(3*index) = fx;
			force(3*index+1) = fy;
			force(3*index+2) = fz;
			if(fixedOrNot == 1){
				fixVertices.push_back(index);
				// cout<<"fix "<<index<<endl;
				fixedIndex++;
			}
			index+=1;
		}
		this->putForceOnTheseVerts.resize(temp.size());
		for(int i=0; i<temp.size(); i++){
			this->putForceOnTheseVerts(i) = temp[i];
			// cout<<TV_k.row(temp[i])<<endl;
		}
	}else{
		cout<<"Check yo self: Force input error, file not found"<<endl;
	}
	// exit(0);
}

void Simulation::staticSolveInitialPosition(vector<int> moveVertices, MatrixXd& TV, MatrixXi& TT, int fv, MatrixXd& B){
	//Time stuff------
 	time_t now = time(0);
 	string dt = ctime(&now);//local time, replace all spaces and new lines
 	dt.erase('\n');
 	replace(dt.begin(), dt.end(), ' ', '-');
 	//-----------------
 	string saveTestToHere = OUTPUT_SAVED_PATH"TestsResults/InitPositionsSolve20moves/"+objectName+"/"+to_string(TT.rows())+"tets@"+tetgen_code+"@"+material_model+"/";
	cout<<"Static Solve initial position"<<endl;
	cout<<saveTestToHere<<endl;
	system(("mkdir -p "+saveTestToHere).c_str());
	system(("(git log -1; echo ''; echo 'Ran Test On:'; date;) >>"+saveTestToHere+"log.txt").c_str());

	double move_amount = 3.5;
	int number_of_moves = 20;
	double dist_moved = 0;
	double curr_youngs = 6.9e6;
	int ignorePastIndex = TV.rows() - moveVertices.size() - fv;

	int count=0;

	int reals_index =0;
	printObj(saveTestToHere, count, TV, TT, B);
	//New Binary Search
	M.setNewYoungsPoissons(curr_youngs, 0.35);

	//Newton Solve for positions
	while(dist_moved<3.5){
		cout<<"	Move next step"<<endl;
		dist_moved += move_amount/number_of_moves;//move step
		cout<<"distance moved"<<dist_moved<<endl;
		staticSolveStepNewtonsMethod(move_amount/number_of_moves, ignorePastIndex, moveVertices, TV, TT);
		printObj(saveTestToHere,  count, TV, TT, B);
		count++;
	}


	VectorXd f;
	f.resize(3*TV.rows());
	calculateElasticForces(f, TV);
	double load_scalar = 0.0;
	Vector3d load(0,0,0);
	//Calculate Load on moving verts
	load.setZero();
	for(unsigned int i=f.size() - 3*moveVertices.size(); i<f.size(); i++){
		load+=f.segment<3>(i);
		i++;
		i++;
	}

	//BELOW: dividing loads by 1000 because we measure forces in Newtons(real data), but calculate in mm based force
	load_scalar = fabs(load(staticSolveDirection)/1000);
	cout<<fabs(load(0)/1000)<<", "<<fabs(load(1)/1000)<<", "<<fabs(load(2)/1000)<<endl;

	ofstream new_force_file;
	new_force_file.open("forces.txt");
	for(int i=0; i<f.rows(); i++){
		new_force_file<<f(i)<<endl;
	}
	ofstream new_init_file;
	new_init_file.open("mesh.off");
	new_init_file<<
    "OFF\n"<<TV.rows()<<" "<<TT.rows()<<" 0\n"<<
    TV.format(IOFormat(FullPrecision,DontAlignCols," ","\n","","","","\n"))<<
    (TT.array()).format(IOFormat(FullPrecision,DontAlignCols," ","\n","4 ","","","\n"));
    new_init_file.close();
}

// static lbfgsfloatval_t evaluateStaticSolveLBFGS(void *s, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step){
// 	Simulation* sim = (Simulation*) s;
// 	unsigned int i=0;
// 	//from x to x_k
// 	for(i=0; i<n; i++){
// 		sim->x_k(i) = x[i];
// 	}
// 	//cout<<"first i "<<i<<endl;

// 	sim->xToTV(sim->x_k, sim->TV_k);
// 	sim->calculateElasticForces(sim->f_k, sim->TV_k);
// 	VectorXd fblock = sim->f_k.head(sim->ignorePastIndex*3);

// 	for(i=0; i<fblock.rows(); i++){
// 		g[i] = -1*fblock(i);
// 	}


// 	//cout<<"second i "<<i<<endl;

// 	double strainE = 0;
// 	for(i=0; i< sim->M.tets.size(); i++){
// 		strainE += sim->M.tets[i].undeformedVol*sim->M.tets[i].energyDensity;
// 	}

// 	double fdstrainE = 0;
// 	double diffVal = 4e-8;
// 	sim->x_k(4) += diffVal;
// 	sim->xToTV(sim->x_k, sim->TV_k);
// 	sim->calculateElasticForces(sim->f_k, sim->TV_k);

// 	for(i=0; i< sim->M.tets.size(); i++){
// 		fdstrainE += sim->M.tets[i].undeformedVol*sim->M.tets[i].energyDensity;
// 	}
// 	cout<<"finite diff E "<<(fdstrainE-strainE)/diffVal<<endl;
// 	cout<<"real g"<<g[4]<<endl;

// 	//cout<<"third i "<<i<<endl;

// 	lbfgsfloatval_t fx = strainE;
// 	//cout<<"here"<<endl;
// 	return fx;
// }

// static int progressStaticSolveLBFGS(void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls){

// 	printf("Iteration %d:\n", k);
//     printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
//     printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
//     printf("\n");
//     return 0;
// }

// void Simulation::staticSolveStepLBFGS(double move_step, int ignorePastIndex, vector<int>& moveVertices, MatrixXd& TV,  MatrixXi& TT){
// 	//Move vertices slightly in x,y,z direction
// 	// [v, v, v..., f, f, ...(m), (m), (m)...]
// 	for(unsigned int i=0; i<moveVertices.size(); i++){
// 		TV.row(TV.rows()-i-1)[staticSolveDirection] += move_step;//move step
// 	}

// 	int N = ignorePastIndex*3;
// 	cout<<"N "<<N<<endl;
// 	cout<<"TV "<<TV.rows()<<endl;
// 	cout<<"x_k "<<x_k.rows()<<endl;
// 	cout<<"igP "<<ignorePastIndex<<endl;

// 	int i, ret =0;
// 	lbfgsfloatval_t fx;
// 	lbfgsfloatval_t *x = lbfgs_malloc(N);

// 	lbfgs_parameter_t param;

// 	if(x ==NULL){
// 		printf("ERROR: Failed to allocate a memory block for variables.\n");
// 	}

// 	//Initialize variables
// 	this->ignorePastIndex = ignorePastIndex;
// 	this->TV_k = TV;
// 	this->f_k.setZero();
// 	this->x_k.setZero();
// 	setTVtoX(x_k, TV_k);

// 	for(i=0; i<N; i++){
// 		x[i] = x_k(i);
// 	}

// 	lbfgs_parameter_init(&param);
// 	// param.gtol = 0.0001;
// 	param.ftol = 0.0001;

// 	ret = lbfgs(N, x, &fx, evaluateStaticSolveLBFGS, progressStaticSolveLBFGS, this, &param);
// 	if(ret<0){
// 		cout<<"ERROR: liblbfgs did not converge in static solve -- code: "<<ret<<endl;
// 		exit(0);
// 	}

// 	// X to TV
// 	xToTV(x_k, TV);

// 	lbfgs_free(x);
// }

void Simulation::xToTV(VectorXd& x, MatrixXd& TV){
	TV.setZero();
	for(unsigned int i=0; i < M.tets.size(); i++){
		Vector4i indices = M.tets[i].verticesIndex;
		TV.row(indices(0)) = Vector3d(x(3*indices(0)), x(3*indices(0)+1), x(3*indices(0) +2));
		TV.row(indices(1)) = Vector3d(x(3*indices(1)), x(3*indices(1)+1), x(3*indices(1) +2));
		TV.row(indices(2)) = Vector3d(x(3*indices(2)), x(3*indices(2)+1), x(3*indices(2) +2));
		TV.row(indices(3)) = Vector3d(x(3*indices(3)), x(3*indices(3)+1), x(3*indices(3) +2));
	}
}

void Simulation::staticSolveStepNewtonsMethod(double move_step, int ignorePastIndex, vector<int>& moveVertices, MatrixXd& TV,  MatrixXi& TT){
	//Move vertices slightly in x,y,z direction
	// [v, v, v..., f, f, ...(m), (m), (m)...]
	for(unsigned int i=0; i<moveVertices.size(); i++){
		TV.row(TV.rows()-i-1)[staticSolveDirection] += move_step;//move step
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

		// cout<<TV.rows()<<endl;
		// cout<<ignorePastIndex<<endl;
		// cout<<forceGradientStaticBlock.rows()<<endl;
		// cout<<forceGradientStaticBlock.cols()<<endl;
		// SparseMatrix<double> forceGradientStaticBlockTranspose = forceGradientStaticBlock.transpose();
		// cout<<(forceGradientStaticBlockTranspose - forceGradientStaticBlock).norm()<<endl;

		//Sparse QR
		SPQR<SparseMatrix<double>> sqr;
		sqr.compute(forceGradientStaticBlock);
		VectorXd deltaX = -1*sqr.solve(fblock);

		// Conj Grad
		// ConjugateGradient<SparseMatrix<double>> cg;
		// cg.compute(forceGradientStaticBlock);
		// if(cg.info() == Eigen::NumericalIssue){
		// 	cout<<"ConjugateGradient numerical issue"<<endl;
		// 	exit(0);
		// }
		// VectorXd deltaX = -1*cg.solve(fblock);

		// // Sparse Cholesky LL^T
		// SimplicialLLT<SparseMatrix<double>> llt;
		// llt.compute(forceGradientStaticBlock);
		// if(llt.info() == Eigen::NumericalIssue){
		// 	cout<<"Possibly using a non- pos def matrix in the LLT method"<<endl;
		// 	exit(0);
		// }
		// VectorXd deltaX = -1*llt.solve(fblock);

		// cout<< (fblock - forceGradientStaticBlock*deltaX).squaredNorm()<<endl;

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
	double strainE = 0;
	for(int i=0; i< M.tets.size(); i++){
		strainE += M.tets[i].energy;
	}
	cout<<"strain E"<<strainE<<endl;
	cout<<"x[0] "<<x(0)<<endl;
	cout<<"x[1] "<<x(1)<<endl;
}

void Simulation::binarySearchYoungs(vector<int> moveVertices, MatrixXd& TV, MatrixXi& TT, int fv, MatrixXd& B){
	cout<<"############Starting Binary Search for Youngs ######################"<<endl;
	this->f_k.resize(3*TV.rows());
	this->x_k.resize(3*TV.rows());


	//REAL VALUES FROM EXPERIMENT
	//dist, load
	vector<pair<double, double>> realLoads =
	{
		{0.100244,	101.073609},
		{0.20048,	184.592875},
		{0.300347,	265.06366},
		{0.399884,	342.551726},
		{0.499522,	417.940362},
		{0.599679,	492.371951},
		{0.700077,	565.750461},
		{0.800373,	638.291505},
		{0.900439,	709.507457},
		{1.000194,	779.029706},
		{1.099773,	846.406884},
		{1.199565,	911.991954},
		{1.299774,	975.713439},
		{1.400103,	1037.549957},
		{1.500424,	1097.35856},
		{1.600455,	1154.854633},
		{1.700085,	1209.47858},
		{1.799595,	1260.718421},
		{1.8996,	1308.14702},
		{1.999928,	1351.454444},
		{2.100283,	1389.728216},
		{2.200424,	1421.50609},
		{2.300369,	1444.923954},
		{2.39994,	1458.153129},
		{2.499569,	1460.703522},
		{2.599659,	1453.403456},
		{2.699964,	1438.137844},
		{2.800305,	1427.789346},
		{2.900472,	1412.155562},
		{3.000338,	1393.765691}
	};

	// vector<pair<double, double>> realLoads =
	// {
	// 	{0.06, 37.0006},
	// 	{0.66, 405.748},
	// 	{1.26, 772.218},
	// 	{1.86, 1136.43},
	// 	{2.46, 1498.42},
	// 	{3.06, 1858.2},
	// 	{3.66, 2215.81},
	// 	{4.26, 2571.27},
	// 	{4.86, 2924.6},
	// 	{5.46, 3275.83},
	// 	{6.06, 3625}

	// };



	//#############################

	//Time stuff------
 	time_t now = time(0);
 	string dt = ctime(&now);//local time, replace all spaces and new lines
 	dt.erase('\n');
 	replace(dt.begin(), dt.end(), ' ', '-');
 	//-----------------

	string saveTestToHere = OUTPUT_SAVED_PATH"TestsResults/StaticSolve/"+objectName+"/"+to_string(TT.rows())+"tets@"+tetgen_code+"@"+material_model+"/";
	system(("mkdir -p "+saveTestToHere).c_str());
	system(("(git log -1; echo ''; echo 'Ran Test On:'; date;) >>"+saveTestToHere+"log.txt").c_str());

	ofstream distvLoadFile;
	distvLoadFile.open(saveTestToHere+"#distVLoad.txt");

	ofstream youngsFile;
	youngsFile.open(saveTestToHere+"#Youngs.txt");
	vector<double> derivedYoungs;

	//size of move
	double move_amount = 2.5;
	int number_of_moves = 50;
	double dist_moved = 0;
	double curr_youngs = 1;
	double load_scalar = 0;
	int ignorePastIndex = TV.rows() - moveVertices.size() - fv;

	int count=0;

	int reals_index =0;
	printObj(saveTestToHere, count, TV, TT, B);
	//New Binary Search
	while(dist_moved<move_amount && reals_index<realLoads.size()){
		M.setNewYoungsPoissons(1000000, 0.35);

		//Newton Solve for positions
		while(reals_index<realLoads.size() && (dist_moved<realLoads[reals_index].first && (fabs(dist_moved-realLoads[reals_index].first)>1e-3))){
			cout<<"	Move next step"<<endl;
			dist_moved += move_amount/number_of_moves;//move step
			cout<<"distance moved"<<dist_moved<<endl;
			cout<<(dist_moved<realLoads[reals_index].first && (fabs(dist_moved-realLoads[reals_index].first)>1e-3) )<<endl;
			staticSolveStepNewtonsMethod(move_amount/number_of_moves, ignorePastIndex, moveVertices, TV, TT);
			// staticSolveStepLBFGS(move_amount/number_of_moves, ignorePastIndex, moveVertices, TV, TT);
			count++;
		}
		printObj(saveTestToHere,  count, TV, TT, B);

		//binary search for youngs
		double min_youngs = 8e7;
		double max_youngs = 6e10;
		load_scalar = 0;
		Vector3d load(0,0,0);
		while(fabs(load_scalar-realLoads[reals_index].second)>(realLoads[reals_index].second/1000) && dist_moved<move_amount){
			load_scalar =0;
			curr_youngs = (min_youngs+max_youngs)/2; //just a guess

			if(fabs(curr_youngs - max_youngs)<0.0000001 || fabs(curr_youngs - min_youngs)<0.0000001){
				cout<<"Error: Young's Modulus did not converge."<<endl;
				cout<<"Current y: "<<curr_youngs<<endl;
				cout<<"Max y: "<<max_youngs<<endl;
				cout<<"Min y: "<<min_youngs<<endl;
				cout<<"Message: Check that the correct vertices are being moved and held"<<endl;
				cout<<"Message: Check that the shape is stretched in the correct direction"<<endl;
				exit(0);
			}
			M.setNewYoungsPoissons(curr_youngs, 0.35);
			VectorXd f;
			f.resize(3*TV.rows());
			calculateElasticForces(f, TV);


			//Calculate Load on moving verts
			load.setZero();
			for(unsigned int i=f.size() - 3*moveVertices.size(); i<f.size(); i++){
				load+=f.segment<3>(i);
				i++;
				i++;
			}

			//BELOW: dividing loads by 1000 because we measure forces in Newtons(real data), but calculate in mm based force
			load_scalar = fabs(load(staticSolveDirection)/1e6);

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
		distvLoadFile<<dist_moved<<", "<<fabs(load(0)/1e6)<<", "<<fabs(load(1)/1e6)<<", "<<fabs(load(2)/1e6)<<endl;
		youngsFile<<dist_moved<<", "<<curr_youngs<<", "<<min_youngs<<", "<<max_youngs <<endl;
		reals_index+=1;
	}


	distvLoadFile.close();
	youngsFile.close();
	system(("(echo 'Finished Test On:'; date;)>>"+saveTestToHere+"log.txt").c_str());
	//system("( speaker-test -t sine -f 1000 )& pid=$! ; sleep 5s ; kill -9 $pid");
	exit(0);
}

void Simulation::syntheticTests(vector<int> moveVertices, MatrixXd& TV, MatrixXi& TT, int fv, MatrixXd& B){
	//TODO: move this into the tests section
	cout<<"############Starting Synthetic Load Generation######################"<<endl;
	string saveTestToHere = OUTPUT_SAVED_PATH"TestsResults/SyntheticStaticSolve/";
	ofstream generateLoadsFile;
	generateLoadsFile.open(saveTestToHere+"syntheticGeneratedLoads.txt");

	int setYoungs = 2e6;
	M.setNewYoungsPoissons(setYoungs, 0.35);

	double dist_moved = 0;
	double move_amount = 5;
	double number_of_moves = 100;
	double load_scalar =0;
	double number_of_data_points =10;

	int count = 1;

	while(dist_moved<move_amount){
		double move_step = move_amount/number_of_moves;
		int ignorePastIndex = TV.rows() - moveVertices.size() - fv;
		dist_moved += move_step;

		staticSolveStepNewtonsMethod(move_step, ignorePastIndex, moveVertices, TV, TT);

		if(true /*count%10==0*/){
			cout<<"	print dist/load to file"<<endl;
			//printObj(count/10, TV, TT, B);

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

			load_scalar = fabs(load(staticSolveDirection)/1000);

			cout<<dist_moved<<", "<<load_scalar<<endl;
			generateLoadsFile<<dist_moved<<", "<<load_scalar<<endl;
		}
		cout<<count<<endl;
		count++;
	}

	generateLoadsFile.close();
	cout<<"############End Synthetic Load Generation######################"<<endl;
}

void Simulation::printObj(string printToHere, int numberOfPrints, MatrixXd& TV, MatrixXi& TT, MatrixXd& B){

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

	system(("mkdir -p " +printToHere).c_str());
	igl::writeOBJ(printToHere + to_string(numberOfPrints)+".obj", V_temp, F_temp);

	return;
}
