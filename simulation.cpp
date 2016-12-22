#include "simulation.h"

Simulation::Simulation(void){}

int staticSolveDirection = 0;

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
	VectorXd force;
	force.resize(3*TV.rows());
	force.setZero();

	//setInitPosition(force, fixVertices, moveVertices);

	if(moveVertices.size()>0 or fixVertices.size()>0){
		//cout << "DOING STUFFS" << endl;
		MatrixXd newTV;
		newTV.resize(TV.rows(), TV.cols());
		newTV.setZero();
		MatrixXi newTT;
		newTT.resize(TT.rows(), TT.cols());
		newTT.setZero();
		//cout << "MoveVertsSize :: " << moveVertices.size() << endl;

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

		//cout << "NewMoveIndicesSize :: " << newMoveIndices.size() << endl;

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
		
		}
		
		integrator->initializeIntegrator(deltaT, M, newTV, newTT);
		this->external_force = new_force;
		integrator->fixVertices(newfixIndices);
		// int ignorePastIndex = newTV.rows() - newfixIndices.size();
		// staticSolveNewtonsForces(newTV, newTT, B, new_force, ignorePastIndex);


	}else{
		//cout << "Doing Other Stuffs" << endl;
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

	//int tempRows = integrator->TV.rows();
	
	//cout << "BASELINE ::" << endl;
	//for (int i = 0; i < this->putForceOnTheseVerts.rows(); i++)
	//	cout << integrator->TV.row(integrator->TV.rows()-i-1)(2) << endl;

	//cout << "INDEXES ::" << endl;
	//for (int i = 0; i < this->putForceOnTheseVerts.rows(); i++)
	//	cout << this->putForceOnTheseVerts(i) << endl;

        //cout << "POINTS ::" << endl;
        //for (int i = 0; i < this->putForceOnTheseVerts.rows(); i++)
	//	cout << "X:" << integrator->TV.row(this->putForceOnTheseVerts(i))(0) << " Y:" << integrator->TV.row(this->putForceOnTheseVerts(i))(1) << " Z:" << integrator->TV.row(this->putForceOnTheseVerts(i))(2) << endl;

	//cout << "TT ::" << endl;
        //for (int i = 0; i < this->putForceOnTheseVerts.rows(); i++)
	//	cout << integrator->TT.row(this->putForceOnTheseVerts(i)(0)) << endl;

	//cout << endl;

	while(integrator->simTime<iters){
		integrator->render(this->external_force);
		cout<<"Min Displacement (called maxDisp in code)"<<endl;
		double disp =0;
		for(int i=0; i<this->putForceOnTheseVerts.rows(); i++){
			//disp += integrator->TV.row(this->putForceOnTheseVerts(i))(2);
			if (integrator->TV.row(this->putForceOnTheseVerts(i))(2) < disp)
				disp = integrator->TV.row(this->putForceOnTheseVerts(i))(2);
		}
		if(disp < maxDisp){
			maxDisp = disp;
		}
		cout<<"Current Disp :: "<<disp<<"\n";
		cout<<"Max Disp ::"<<maxDisp<<"\n";
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
	int vertRows = integrator->TV.rows();
	for(int i=0; i<this->putForceOnTheseVerts.rows(); i++){
		//disp += integrator->TV.row(this->putForceOnTheseVerts(i))(2);
		if (integrator->TV.row(this->putForceOnTheseVerts(i))(2) < disp)
			disp = integrator->TV.row(this->putForceOnTheseVerts(i))(2);
	}
	if(disp < maxDisp){
		maxDisp = disp;
	}
	cout<<"Current Disp ::"<<disp<<endl;
	cout<<"Max Disp ::"<<maxDisp<<endl;
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

void Simulation::staticSolveNewtonsForces(MatrixXd& TV, MatrixXi& TT, MatrixXd& B, VectorXd& fixed_forces, int ignorePastIndex){
	cout<<"----------------Static Solve Newtons, Fix Forces-------------"<<endl;
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
	cout<<x<<endl;
	int NEWTON_MAX = 100, k=0;
	for(k=0; k<NEWTON_MAX; k++){
		xToTV(x, TV);

		calculateForceGradient(TV, forceGradient);
		calculateElasticForces(f, TV);
		for(unsigned int i=0; i<fixed_forces.rows(); i++){
			if(abs(fixed_forces(i))>0.000001){
				if(i>3*ignorePastIndex){
					cout<<"Problem Check simulation.cpp file"<<endl;
					cout<<ignorePastIndex<<endl;
					cout<<i<<" - "<<fixed_forces(i)<<endl;
					exit(0);
				}
				f(i) = fixed_forces(i);
			}
		}
		
		//Block forceGrad and f to exclude the fixed verts
		forceGradientStaticBlock = forceGradient.block(0,0, 3*(ignorePastIndex), 3*ignorePastIndex);
		VectorXd fblock = f.head(ignorePastIndex*3);

		//Sparse QR 
		// SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> sqr;
		// sqr.compute(forceGradientStaticBlock);
		// VectorXd deltaX = -1*sqr.solve(fblock);

		// Conj Grad
		ConjugateGradient<SparseMatrix<double>> cg;
		cg.compute(forceGradientStaticBlock);
		if(cg.info() == Eigen::NumericalIssue){
			cout<<"ConjugateGradient numerical issue"<<endl;
			exit(0);
		}
		VectorXd deltaX = -1*cg.solve(fblock);

		// // Sparse Cholesky LL^T
		// SimplicialLLT<SparseMatrix<double>> llt;
		// llt.compute(forceGradientStaticBlock);
		// if(llt.info() == Eigen::NumericalIssue){
		// 	cout<<"Possibly using a non- pos def matrix in the LLT method"<<endl;
		// 	exit(0);
		// }
		// VectorXd deltaX = -1*llt.solve(fblock);

		x.segment(0,3*(ignorePastIndex))+=deltaX;
		cout<<"		Newton Iter "<<k<<endl;

		if(x != x){
			cout<<"NAN"<<endl;
			exit(0);
		}
		cout<<"fblock square norm"<<endl;
		cout<<fblock.squaredNorm()/fblock.size()<<endl;
		double strainE = 0;
		for(int i=0; i< M.tets.size(); i++){
			strainE += M.tets[i].undeformedVol*M.tets[i].energyDensity;
		}
		cout<<strainE<<endl;
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
		strainE += M.tets[i].undeformedVol*M.tets[i].energyDensity;
	}
	cout<<"strain E"<<strainE<<endl;
	cout<<"x[0] "<<x(0)<<endl;
	cout<<"x[1] "<<x(1)<<endl;
	exit(0);				

	cout<<"-------------------"<<endl;
}

void Simulation::setInitPosition(VectorXd& force, vector<int>& fixVertices, vector<int>& moveVertices){
	vector<int> temp;
	ifstream forceInputFile (TUTORIAL_SHARED_PATH "shared/"+objectName+".txt");
	if(forceInputFile.is_open()){
		string line;
		int index =0;
		int fixedIndex = 0;
		while(getline(forceInputFile, line)){
			istringstream iss(line);
			double fx, fy, fz;
			int fixedOrNot; //1 is fixed, 0 not fixed
			if(!(iss >> fx >> fy >> fz >> fixedOrNot)){break;}
			if(abs(fx + fy + fz)>0){
				temp.push_back(index - fixedIndex);
				cout<<index<<endl;
			}
			force(3*index) = fx;
			force(3*index+1) = fy;
			force(3*index+2) = fz;
			if(fixedOrNot == 1){
				// cout<<fx<<" "<<fy<<" "<<fz<<" "<<fixedOrNot<<endl;
				// cout<<index<<endl;
				fixVertices.push_back(index);
				fixedIndex++;
			}
			index+=1;
		}
		this->putForceOnTheseVerts.resize(temp.size());
		for(int i=0; i<temp.size(); i++){
			this->putForceOnTheseVerts(i) = temp[i];
			//moveVertices.push_back(temp[i]);
			cout << "VERTIND :: " << temp[i] << endl;
			//cout<<TV_k.row(temp[i])<<endl;
		}
	}else{
		cout<<"Check yo self: Force input error, file not found"<<endl;
	}
	// cout<<"fixed verts"<<endl;
	// for(int i=0; i<fixVertices.size(); i++){
	// 	cout<<fixVertices[i]<<endl;
	// }
	// cout<<"forces####"<<endl;
	// for(int i=0; i<force.rows(); i+=3){
	// 	cout<<force(i)<<" "<<force(i+1)<<" "<<force(i+2)<<" "<<endl;
	// }
}


static lbfgsfloatval_t evaluateStaticSolveLBFGS(void *s, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step){
	Simulation* sim = (Simulation*) s;
	unsigned int i=0;
	//from x to x_k
	for(i=0; i<n; i++){
		sim->x_k(i) = x[i];
	}
	//cout<<"first i "<<i<<endl;

	sim->xToTV(sim->x_k, sim->TV_k);
	sim->calculateElasticForces(sim->f_k, sim->TV_k);
	VectorXd fblock = sim->f_k.head(sim->ignorePastIndex*3);

	for(i=0; i<fblock.rows(); i++){
		g[i] = -1*fblock(i);
	}


	//cout<<"second i "<<i<<endl;

	double strainE = 0;
	for(i=0; i< sim->M.tets.size(); i++){
		strainE += sim->M.tets[i].undeformedVol*sim->M.tets[i].energyDensity;
	}

	double fdstrainE = 0;
	double diffVal = 4e-8;
	sim->x_k(4) += diffVal;
	sim->xToTV(sim->x_k, sim->TV_k);
	sim->calculateElasticForces(sim->f_k, sim->TV_k);

	for(i=0; i< sim->M.tets.size(); i++){
		fdstrainE += sim->M.tets[i].undeformedVol*sim->M.tets[i].energyDensity;
	}
	cout<<"finite diff E "<<(fdstrainE-strainE)/diffVal<<endl;
	cout<<"real g"<<g[4]<<endl;

	//cout<<"third i "<<i<<endl;

	lbfgsfloatval_t fx = strainE;
	//cout<<"here"<<endl;
	return fx;
}

static int progressStaticSolveLBFGS(void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls){	

	printf("Iteration %d:\n", k);
    printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
    printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    printf("\n");
    return 0;	
}

void Simulation::staticSolveStepLBFGS(double move_step, int ignorePastIndex, vector<int>& moveVertices, MatrixXd& TV,  MatrixXi& TT){
	//Move vertices slightly in x,y,z direction
	// [v, v, v..., f, f, ...(m), (m), (m)...]
	for(unsigned int i=0; i<moveVertices.size(); i++){
		TV.row(TV.rows()-i-1)[staticSolveDirection] += move_step;//move step
	}
	

	int N = ignorePastIndex*3;
	cout<<"N "<<N<<endl;
	cout<<"TV "<<TV.rows()<<endl;
	cout<<"x_k "<<x_k.rows()<<endl;
	cout<<"igP "<<ignorePastIndex<<endl;

	int i, ret =0;
	lbfgsfloatval_t fx;
	lbfgsfloatval_t *x = lbfgs_malloc(N);

	lbfgs_parameter_t param;

	if(x ==NULL){
		printf("ERROR: Failed to allocate a memory block for variables.\n");
	}

	//Initialize variables
	this->ignorePastIndex = ignorePastIndex;
	this->TV_k = TV;
	this->f_k.setZero();
	this->x_k.setZero();
	setTVtoX(x_k, TV_k);

	for(i=0; i<N; i++){
		x[i] = x_k(i);
	}

	lbfgs_parameter_init(&param);
	// param.gtol = 0.0001;
	param.ftol = 0.0001;

	ret = lbfgs(N, x, &fx, evaluateStaticSolveLBFGS, progressStaticSolveLBFGS, this, &param);
	if(ret<0){
		cout<<"ERROR: liblbfgs did not converge in static solve -- code: "<<ret<<endl;
		exit(0);
	}

	// X to TV
	xToTV(x_k, TV);

	lbfgs_free(x);
}

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
		// SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> sqr;
		// sqr.compute(forceGradientStaticBlock);
		// VectorXd deltaX = -1*sqr.solve(fblock);

		// Conj Grad
		ConjugateGradient<SparseMatrix<double>> cg;
		cg.compute(forceGradientStaticBlock);
		if(cg.info() == Eigen::NumericalIssue){
			cout<<"ConjugateGradient numerical issue"<<endl;
			exit(0);
		}
		VectorXd deltaX = -1*cg.solve(fblock);

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
		strainE += M.tets[i].undeformedVol*M.tets[i].energyDensity;
	}
	cout<<"strain E"<<strainE<<endl;
	cout<<"x[0] "<<x(0)<<endl;
	cout<<"x[1] "<<x(1)<<endl;
	exit(0);				
}

void Simulation::binarySearchYoungs(vector<int> moveVertices, MatrixXd& TV, MatrixXi& TT, int fv, MatrixXd& B){
	cout<<"############Starting Binary Search for Youngs ######################"<<endl;
	this->f_k.resize(3*TV.rows());
	this->x_k.resize(3*TV.rows());


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

	//############Spring synth test
	// vector<pair<double, double>> realLoads =
	// {
	// 	{0.5, 5.56078},
	// 	{1, 11.1308},
	// 	{1.5, 16.7229},
	// 	{2, 22.3498},
	// 	{2.5, 28.024},
	// 	{3, 33.7579},
	// 	{3.5, 39.5637}

	// };

	// vector<pair<double, double>> realLoads = 
	// {
	// 	{0.1001330556,	1.685690278},
	// 	//{0.2004038889,	3.482108056},
	// 	{0.3003497222,	5.150964444},
	// 	//{0.4000427778,	6.832918333},
	// 	{0.4998455556,	8.463683333},
	// 	//{0.5998727778,	10.14414472},
	// 	{0.7000202778,	11.77661167},
	// 	//{0.8002011111,	13.39418889},
	// 	{0.9003483333,	15.01527306},
	// 	//{1.000226111,	16.67514917},
	// 	{1.099982222,	18.26396861},
	// 	//{1.199848889,	19.89806639},
	// 	{1.300013056,	21.48296167},
	// 	//{1.400086667,	23.05293972},
	// 	{1.500295	,	24.63098806},
	// 	//{1.600397778,	26.21527389},
	// 	{1.700163333,	27.72697556},
	// 	//{1.799876944,	29.241815},
	// 	{1.899856111,	30.74729222},
	// 	//{1.999977778,	32.20479444},
	// 	{2.100109444,	33.65643083},
	// 	//{2.200305556,	35.04660222},
	// 	{2.300320278,	36.41194389},
	// 	//{2.400083333,	37.76987889},
	// 	{2.499865	,	39.02258972},
	// 	//{2.599960833,	40.309035},
	// 	{2.700056667,	41.50979361},
	// 	//{2.800188611,	42.70812972},
	// 	{2.900375	,	43.87573139},
	// 	//{3.000270833,	44.95603361},
	// 	{3.099961667,	45.98499917},
	// 	//{3.199837222,	46.96560389},
	// 	{3.299879722,	47.94410972},
	// 	//{3.400075556,	48.88907111},
	// 	{3.500241111,	49.72337972},
	// 	//{3.600337778,	50.53218222},
	// 	{3.700195278,	51.32930417},
	// 	//{3.799969444,	52.08023306},
	// 	{3.899893056,	52.78361528},
	// 	//{4.000053056,	53.19902667},
	// 	{4.100091111,	53.94501889},
	// 	//{4.200299722,	54.59250778},
	// 	{4.300336389,	55.13545944}
	// };

	vector<pair<double, double>> realLoads = 
	{
		{0.01708,	0.25298},
		{0.02459,	0.42219},
		{0.03355,	0.61954},
		{0.04124,	0.79588},
		{0.05063,	0.91467},
		{0.05789,	1.06907},
		{0.06702,	1.18365},
		{0.07437,	1.34859},
		{0.08393,	1.48492},
		{0.0911	,	1.59278},
		{0.10041,	1.74828},
		{0.10766,	1.89247},
		{0.11731,	2.04129},
		{0.12457,	2.18455},
		{0.13379,	2.35536},
		{0.14105,	2.50726},
		{0.15061,	2.71388},
		{0.15804,	2.88332},
		{0.16709,	2.97436},
		{0.17452,	3.12331},
		{0.18374,	3.2688},
		{0.19168,	3.40847},
		{0.20047,	3.59617},
		{0.20807,	3.67408},
		{0.21695,	3.85588},
		{0.22515,	3.93631},
		{0.23368,	4.10039},
		{0.24162,	4.20311},
		{0.25016,	4.40673},
		{0.25861,	4.55356},
		{0.2669	,	4.76446},
		{0.27509,	4.88167},
		{0.28329,	5.0596},
		{0.292	,	5.14403},
		{0.30011,	5.33197},
		{0.30865,	5.4899},
		{0.3165	,	5.63556},
		{0.32538,	5.7101},
		{0.33315,	5.86506},
		{0.34212,	5.97343},
		{0.34971,	6.0844},
		{0.35876,	6.27656},
		{0.36636,	6.48709},
		{0.37558,	6.63124},
		{0.38293,	6.7856},
		{0.39206,	6.90234},
		{0.39958,	7.07678},
		{0.40897,	7.20909},
		{0.41622,	7.35599},
		{0.42545,	7.51178},
		{0.43279,	7.56089},
		{0.44227,	7.74383},
		{0.44944,	7.85254},
		{0.45874,	8.02209},
		{0.46609,	8.16599},
		{0.47556,	8.35534},
		{0.48291,	8.48604},
		{0.49213,	8.67024},
		{0.49947,	8.80885},
		{0.50886,	8.95673},
		{0.51637,	9.09092},
		{0.52551,	9.25478},
		{0.53285,	9.33346},
		{0.54207,	9.52713},
		{0.54993,	9.60187},
		{0.55881,	9.75249},
		{0.56641,	9.84894},
		{0.57529,	10.03423},
		{0.58348,	10.2213},
		{0.59211,	10.37955},
		{0.59996,	10.51487},
		{0.60841,	10.68351},
		{0.61695,	10.7722},
		{0.62523,	10.92284},
		{0.6336	,	11.07693},
		{0.64163,	11.21744},
		{0.65034,	11.34825},
		{0.65828,	11.48543},
		{0.66707,	11.57553},
		{0.67475,	11.71481},
		{0.68372,	11.92454},
		{0.6914	,	12.07428},
		{0.70062,	12.22542},
		{0.70797,	12.38059},
		{0.71702,	12.55692},
		{0.72453,	12.68159},
		{0.73401,	12.83822},
		{0.74126,	12.9322},
		{0.7504	,	13.09948},
		{0.75774,	13.21157},
		{0.76739,	13.36231},
		{0.77456,	13.43525},
		{0.7837	,	13.61628},
		{0.79104,	13.70739},
		{0.8006	,	13.91004},
		{0.80803,	14.04445},
		{0.817	,	14.26001},
		{0.82442,	14.3797},
		{0.83382,	14.53185},
		{0.8415	,	14.67118},
		{0.85046,	14.81958},
		{0.85789,	14.96015},
		{0.86703,	15.06137},
		{0.87505,	15.21407},
		{0.88376,	15.31573},
		{0.89145,	15.42412},
		{0.90024,	15.5607},
		{0.90852,	15.75775},
		{0.91698,	15.91721},
		{0.925	,	16.04398},
		{0.93345,	16.20599},
		{0.94199,	16.3585},
		{0.95019,	16.49678},
		{0.95856,	16.64456},
		{0.96658,	16.81985},
		{0.97537,	16.95093},
		{0.98332,	17.03351},
		{0.99202,	17.15232},
		{0.99979,	17.30681},
		{1.00876,	17.44692},
		{1.01644,	17.55606},
		{1.02558,	17.7202},
		{1.03301,	17.88259},
		{1.04206,	18.01454},
		{1.04966,	18.13639},
		{1.05896,	18.33371},
		{1.06622,	18.50289},
		{1.07544,	18.63088},
		{1.08287,	18.7809},
		{1.09234,	18.94007},
		{1.09952,	19.02288},
		{1.10874,	19.123},
		{1.11608,	19.2554},
		{1.12564,	19.41167},
		{1.13281,	19.5554},
		{1.14212,	19.68067},
		{1.14938,	19.81684},
		{1.15886,	20.00899},
		{1.1662	,	20.10608},
		{1.1755	,	20.2971},
		{1.18276,	20.44939},
		{1.19215,	20.61627},
		{1.19975,	20.77633},
		{1.20889,	20.8633},
		{1.21623,	20.96754},
		{1.22537,	21.09533},
		{1.23339,	21.2446},
		{1.24219,	21.39535},
		{1.24987,	21.53676},
		{1.25858,	21.67618},
		{1.26686,	21.78997},
		{1.2754	,	21.92999},
		{1.28334,	22.09986},
		{1.29179,	22.28079},
		{1.30024,	22.41031},
		{1.30853,	22.55981},
		{1.31698,	22.70258},
		{1.32492,	22.78282},
		{1.33371,	22.91223},
		{1.34157,	23.04497},
		{1.35045,	23.18036},
		{1.35805,	23.38344},
		{1.3671	,	23.52233},
		{1.37469,	23.66914},
		{1.384	,	23.8488},
		{1.39126,	23.96028},
		{1.40039,	24.05413},
		{1.40782,	24.26275},
		{1.41738,	24.38892},
		{1.42456,	24.49606},
		{1.43378,	24.5844},
		{1.44103,	24.742},
		{1.4506	,	24.86861},
		{1.45794,	25.03281},
		{1.46708,	25.21182},
		{1.47442,	25.32112},
		{1.4839	,	25.50219},
		{1.49141,	25.65068},
		{1.50046,	25.77975},
		{1.5078	,	25.84671},
		{1.51711,	26.04051},
		{1.52496,	26.15532},
		{1.53376,	26.29304},
		{1.54136,	26.37214},
		{1.55032,	26.52069},
		{1.55852,	26.62709},
		{1.56714,	26.79148},
		{1.57482,	26.96355},
		{1.58353,	27.18522},
		{1.59199,	27.26762},
		{1.60027,	27.43705},
		{1.60846,	27.47519},
		{1.61666,	27.63762},
		{1.62537,	27.74542},
		{1.63348,	27.91935},
		{1.64202,	27.98132},
		{1.64987,	28.15039},
		{1.65875,	28.29656},
		{1.66661,	28.4345},
		{1.6754	,	28.58163},
		{1.68309,	28.76842},
		{1.69205,	28.90943},
		{1.69973,	29.04132},
		{1.70887,	29.1632},
		{1.7163	,	29.28292},
		{1.72543,	29.37723},
		{1.74234,	29.64538},
		{1.74951,	29.78411},
		{1.75873,	29.90338},
		{1.76607,	30.08482},
		{1.77564,	30.21732},
		{1.78281,	30.37564},
		{1.79211,	30.54395},
		{1.79937,	30.66848},
		{1.80902,	30.81138},
		{1.81619,	30.90787},
		{1.8255	,	31.03343},
		{1.83276,	31.12435},
		{1.84215,	31.27966},
		{1.84966,	31.38907},
		{1.85888,	31.52755},
		{1.86614,	31.66244},
		{1.87545,	31.81878},
		{1.88321,	31.95701},
		{1.89218,	32.14753},
		{1.89969,	32.30344},
		{1.90866,	32.4484},
		{1.91677,	32.52985},
		{1.92548,	32.66881},
		{1.93325,	32.76123},
		{1.94187,	32.89001},
		{1.95024,	32.97517},
		{1.95869,	33.14239},
		{1.96672,	33.2887},
		{1.975	,	33.46284},
		{1.98362,	33.59881},
		{1.99173,	33.75123},
		{2.00036,	33.89053}
	};
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
	double move_amount = 2.1;
	int number_of_moves = 500;
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
		while(reals_index<realLoads.size() && (dist_moved<realLoads[reals_index].first && (abs(dist_moved-realLoads[reals_index].first)>1e-5))){
			cout<<"	Move next step"<<endl;
			dist_moved += move_amount/number_of_moves;//move step
			cout<<"distance moved"<<dist_moved<<endl;
			cout<<(dist_moved<realLoads[reals_index].first && (abs(dist_moved-realLoads[reals_index].first)>1e-5) )<<endl;
			// staticSolveStepNewtonsMethod(move_amount/number_of_moves, ignorePastIndex, moveVertices, TV, TT);	
			staticSolveStepLBFGS(move_amount/number_of_moves, ignorePastIndex, moveVertices, TV, TT);
			count++;
		}
		printObj(saveTestToHere,  count, TV, TT, B);

		//binary search for youngs
		double min_youngs = 800000;
		double max_youngs = 60000000;
		load_scalar = 0;
		Vector3d load(0,0,0);
		while(abs(load_scalar-realLoads[reals_index].second)>(realLoads[reals_index].second/1000) && dist_moved<move_amount){
			load_scalar =0;
			curr_youngs = (min_youngs+max_youngs)/2; //just a guess
			
			if(abs(curr_youngs - max_youngs)<0.0000001 || abs(curr_youngs - min_youngs)<0.0000001){
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
			load_scalar = abs(load(staticSolveDirection)/1000);
			
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
		distvLoadFile<<dist_moved<<", "<<abs(load(0)/1000)<<", "<<abs(load(1)/1000)<<", "<<abs(load(2)/1000)<<endl;
		youngsFile<<dist_moved<<", "<<curr_youngs<<", "<<min_youngs<<", "<<max_youngs <<endl;
		reals_index+=1;
	}


	distvLoadFile.close();
	youngsFile.close();
	system(("(echo 'Finished Test On:'; date;)>>"+saveTestToHere+"log.txt").c_str());
	//system("( speaker-test -t sine -f 1000 )& pid=$! ; sleep 5s ; kill -9 $pid");
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

			load_scalar = abs(load(staticSolveDirection)/1000);

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
	cout<<printToHere<<endl;
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

