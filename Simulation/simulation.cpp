#include "simulation.h"


Simulation::Simulation(void){}

int staticSolveDirection = 2;

int Simulation::initializeSimulation(double deltaT, int iterations, char method, MatrixXi& TT, MatrixXd& TV, MatrixXd& B, vector<int>& moveVertices, vector<int> fixVertices, double youngs, double poissons){
	iters = iterations;
	if (method =='e'){
	//	integrator = new Verlet();
		cout<<"Initialized Verlet"<<endl;
		exit(0);
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
	cout<<tetgen_code<<endl;
	VectorXd force;
	force.resize(3*TV.rows());
	force.setZero();
	TV_k = TV;

	// setInitPosition(force, fixVertices, moveVertices);

	// fixVertices.push_back(4);
	// fixVertices.push_back(3);
	// fixVertices.push_back(2);
	// moveVertices.push_back(1);


	//BEAM SPRING FIXING vertices and MOVING VERTICES COMMENTED OUT ^ CAUSE IT DOESN"T WORK
	// for(int i=0; i<TV.rows(); i++){
	// 	if(TV.row(i)[0] < 0.1 && TV.row(i)[2]<-3.0){
	// 		moveVertices.push_back(i);
	// 	}
	// 	if(TV.row(i)[0] > 140.1){
	// 		fixVertices.push_back(i);
	// 	}
	// }

	//Analytical Beam fix and moved
	for(int i=0; i<TV.rows(); i++)
	{
		if(TV.row(i)[1]>=112){
			fixVertices.push_back(i);
		}
		if(TV.row(i)[1]<=-192 && TV.row(i)[2]<=0){
			moveVertices.push_back(i);
		}
	}

	//FIXING vertices and MOVING VERTICES COMMENTED OUT ^ CAUSE IT DOESN"T WORK
	// for(int i=0; i<TV.rows(); i++){
	// 	if(TV.row(i)[1] > -1.5){
	// 		moveVertices.push_back(i);
	// 	}
	// 	if(TV.row(i)[1] < -40.5){
	// 		fixVertices.push_back(i);
	// 	}
	// }
	// double amountofforce = 1e8/moveVertices.size();
	// for(int i=0; i<moveVertices.size(); i++){
	// 	force(3*moveVertices[i] + 1) -= amountofforce;
	// }
	//-------------------

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
			//[v, v, v, v...., m, m, m...., f, f, f...,f]
			if(!flag){
				vertexNewIndices.push_back(i);
			}
		}
		//re-index move verts
		for(unsigned int j=0; j<moveVertices.size(); j++){
			vertexNewIndices.push_back(moveVertices[j]);
		}
		//re-index fixed verts
		for(unsigned int j=0; j<fixVertices.size(); j++){
			vertexNewIndices.push_back(fixVertices[j]);
		}

		//new indices for the moving verts
		vector<int> newMoveIndices;
		for(unsigned int i= vertexNewIndices.size() - (moveVertices.size()+fixVertices.size()); i<(vertexNewIndices.size() - fixVertices.size()); i++){
			newMoveIndices.push_back(i);
		}
		//these are the new indices for the fixed verts
		vector<int> newfixIndices;
		for(unsigned int i= vertexNewIndices.size() - fixVertices.size(); i<vertexNewIndices.size(); i++){
			newfixIndices.push_back(i);
		}

		VectorXd new_force;
		new_force.resize(3*TV.rows());
		reIndexTVandTT(vertexNewIndices, fixVertices.size(), moveVertices.size(), TV, TT, force, newTV, newTT, new_force);

		igl::barycenter(newTV, newTT, B);
		//Initialize Solid Mesh
		M.initializeMesh(newTT, newTV, youngs, poissons);
		if(moveVertices.size() != 0){
			moveVertices = newMoveIndices;
			this->moveVerticesStore = newMoveIndices;
			//RECOMMENT
			ifstream meshFile(OUTPUT_SAVED_PATH "TestsResults/Damping/"+objectName+"@"+tetgen_code+".mesh");
			cout<<OUTPUT_SAVED_PATH "TestsResults/Damping/"+objectName+"@"+tetgen_code+".mesh"<<endl;
			if(meshFile.good()){
				igl::readMESH(OUTPUT_SAVED_PATH "TestsResults/Damping/"+objectName+"@"+tetgen_code+".mesh", newTV, newTT, TF);
			}else{
				cout<<"APPLYING STATIC POSITIONS"<<endl;
				applyStaticPositions(newTV, newTT, B, new_force, newMoveIndices, newfixIndices);
				// applyStaticForces(newTV, newTT, B, new_force, newMoveIndices, newfixIndices);
			}

		}

		integrator->initializeIntegrator(deltaT, M, newTV, newTT);

		double forceToSet = 44480000.0; // units: gmm/(ss)
		cout<<"ext force"<<endl;
		cout<<new_force.rows()<<endl;
		for(unsigned int i= this->external_force.size() - (moveVertices.size()+fixVertices.size()); i<(vertexNewIndices.size() - fixVertices.size());){
			i++;
			i++;
			new_force(i) = 3*forceToSet/(moveVertices.size());
			i++;
		}
		this->external_force = -1*new_force;

		integrator->fixVertices(newfixIndices);
		integrator->moveVertices(this->moveVerticesStore);

	}else{
		igl::barycenter(TV, TT, B);
		M.initializeMesh(TT, TV, youngs, poissons);
		integrator->initializeIntegrator(deltaT, M, TV, TT);
		this->external_force = force;
		integrator->fixVertices(fixVertices);
	}

	sB = &B;
	cout<<"NUMBER OF FIXED: "<< fixVertices.size()<<endl;
	cout<<"NUMBER OF MOVING: "<<this->moveVerticesStore.size()<<endl;
	return 1;
}

void Simulation::applyExternalForces(){


	this->external_force.setZero();
}

void Simulation::headless(){
	int printcount =0;
	ofstream dampingPositionFile;
	//	cout<<OUTPUT_SAVED_PATH"TestsResults/Damping/"<<endl;
	dampingPositionFile.open(OUTPUT_SAVED_PATH"TestsResults/Damping/"+solver+"_Y:"+to_string(youngs)+"@R:"+to_string(rayleighCoeff)+"@step"+to_string(integrator->h)+"@"+to_string(integrator->TT.rows())+"tets@"+tetgen_code+"@"+"position.txt");

	integrator->v_old.setZero();
	integrator->f.setZero();
	printDesigns(printcount, integrator->simTime);
	double maxYVel = 0;
	while(integrator->simTime < iters){
		integrator->render(this->external_force);
		//RECOMMENT
		// double z_pos =0;
		// // z_pos = integrator->TV.row(1121)[2] + 0.5;//For pRa5
		// for(int i=0; i<this->moveVerticesStore.size(); i++){
		// 	z_pos += integrator->x_old(3*this->moveVerticesStore[i]+2); //avg in the y direction
		// }
		// dampingPositionFile<<integrator->simTime<<", "<<z_pos/this->moveVerticesStore.size()<<endl;
		double z_pos = 0;
		for(int i=0; i<this->moveVerticesStore.size(); i++){
			z_pos += integrator->TV.row(this->moveVerticesStore[0])(1);
		}
		z_pos /= this->moveVerticesStore.size();

		dampingPositionFile<<integrator->simTime<<", "<<z_pos<<endl;
		// double yvel = printOptimizationOutput();
		// if(yvel>maxYVel)
		// 	maxYVel = yvel;
		if(integrator->simTime%100==0){
			printDesigns(printcount, integrator->simTime);
			// igl::writeMESH(OUTPUT_SAVED_PATH "TestsResults/temp/"+solver+"/"+to_string(printcount), integrator->TV, integrator->TT, TF);
			printcount += 1;
		}
	}
	dampingPositionFile.close();

}

void Simulation::printDesigns(int printcount, int simTime){
	string saveTestsHere = OUTPUT_SAVED_PATH"TestsResults/Analytical/"+solver+"_Y:"+to_string(youngs)+"@R:"+to_string(rayleighCoeff)+"@step"+to_string(integrator->h)+"@"+to_string(integrator->TT.rows())+"tets@"+tetgen_code+"@"+objectName+"/";
	printObj(saveTestsHere, printcount, integrator->TV, integrator->TT, *sB);
	cout<<printcount<<endl;
}

double Simulation::printOptimizationOutput(){
	// double avgmove = 0.0;
	// for(int i = 0; i < this->moveVerticesStore.size(); i++){
	// 	avgmove += integrator->TV.row(this->moveVerticesStore[i])(1);
	// }
	//
	// optimizationFile<<integrator->simTime <<avgmove/this->moveVerticesStore.size()<<endl;
	// cout<<"Spring has moved to here:"<<endl;
	// cout<<avgmove/this->moveVerticesStore.size()<<"\n";

	double avgV = 0.0;
	for(int i=0; i<this->moveVerticesStore.size(); i++){
		avgV += integrator->v_old(3*this->moveVerticesStore[i]+1); //avg in the y direction
	}
	optimizationFile<< integrator->simTime <<", "<< avgV/this->moveVerticesStore.size()<<"\n";
	cout<<"Moving this fast:"<<endl;
	cout<<avgV/this->moveVerticesStore.size()<<endl;
	return avgV/this->moveVerticesStore.size();
}

bool Simulation::render(){
	//These changes are for the spring
	integrator->render(this->external_force);
	if(integrator->simTime%100 == 0){
		printDesigns(integrator->simTime, integrator->simTime);
	}
	// printOptimizationOutput();

	return true;
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
void Simulation::applyStaticPositions(MatrixXd& TV, MatrixXi& TT, MatrixXd& B, VectorXd& fixed_forces, vector<int>& moveVertices, vector<int>& fixVertices){
	cout<<"***SETTING INITIAL POSITIONS***"<<endl;
	double movePercentOfSpringLength = .1;
	cout<<"		DIRECTION of Static position solve: "<<staticSolveDirection<<endl;
	int direction = staticSolveDirection; //-y
	double designZMax = 0.0;
	double designZMin = 0.0;
	double designXMax = 0.0;
	double designXMin = 0.0;
	double designYMax = 0.0;
	double designYMin = 0.0;
	for(int i=0; i<TV.rows(); ++i){
		if(TV.row(i)[2] > designZMax){
			designZMax = TV.row(i)[2];
		}
		if(TV.row(i)[2] < designZMin){
			designZMin = TV.row(i)[2];
		}
		if(TV.row(i)[1] > designYMax){
			designYMax = TV.row(i)[1];
		}
		if(TV.row(i)[1] < designYMin){
			designYMin = TV.row(i)[1];
		}
		if(TV.row(i)[0] > designXMax){
			designXMax = TV.row(i)[0];
		}
		if(TV.row(i)[0] < designXMin){
			designXMin = TV.row(i)[0];
		}
	}
	Vector3d designMaxes(designXMax, designYMax, designZMax);
	Vector3d designMins(designXMin, designYMin, designZMin);
	cout<<designXMax<<", "<<designXMin<<endl;
	cout<<designYMax<<", "<<designYMin<<endl;
	cout<<designZMax<<", "<<designZMin<<endl;

	// double distance_to_move = (designZMax - designZMin)*movePercentOfSpringLength;
	double distance_to_move = 58.52;
	int number_of_moves = 1000;
	double step_size = distance_to_move/number_of_moves;
	cout<<"STEP SIZE"<<endl;
	cout<<step_size<<endl;
	int ignorePastIndex = TV.rows() - moveVertices.size() - fixVertices.size();
	double amount_moved = 0;
	int c=0;
	while(amount_moved < distance_to_move){
		//Move vertices slightly in x,y,z direction
		// [v, v, v..., m, m, ...f, f, f...]
		for(unsigned int i=0; i<moveVertices.size(); i++){
			TV.row(moveVertices[i])[abs(direction)] += (direction/direction)*step_size;//step
		}
		if(c%10==0){
			printObj(OUTPUT_SAVED_PATH"TestsResults/Analytical/", c, TV, TT, B);
		}
		staticSolveNewtonsPosition(TV, TT, B, moveVertices, ignorePastIndex, c);
		c++;
		amount_moved+=step_size;
	}
	printObj(OUTPUT_SAVED_PATH "TestsResults/Damping/", c, TV, TT, B);
	cout<<"WRITE MESH HERE: "<< OUTPUT_SAVED_PATH "TestsResults/Damping/"+objectName+"@"+tetgen_code+".mesh"<<endl;
	igl::writeMESH(OUTPUT_SAVED_PATH "TestsResults/Damping/"+objectName+"@"+tetgen_code+".mesh", TV, TT, TF);
}
void Simulation::staticSolveNewtonsPosition(MatrixXd& TV, MatrixXi& TT, MatrixXd& B, vector<int>& moveVertices, int ignorePastIndex, int step){
	cout<<"------------Static Solve Newtons POSITION - Iteration"<< step<<"--------------"<<endl;
	//Newtons method static solve for minimum Strain E
	SparseMatrix<double> forceGradient;
	SparseMatrix<double> forceGradientStaticBlock;
	SparseMatrix<double> DeletionMatrix;
	VectorXd fblock, fblocktail;
	forceGradient.resize(3*TV.rows(), 3*TV.rows());

	forceGradientStaticBlock.resize((3*ignorePastIndex + 2*moveVertices.size()), (3*ignorePastIndex + 2*moveVertices.size()));
	DeletionMatrix.resize((3*ignorePastIndex + 2*moveVertices.size()), (3*ignorePastIndex + 3*moveVertices.size()));
	fblock.resize(3*ignorePastIndex + 2*moveVertices.size());
	fblocktail.resize(2*moveVertices.size());

	// fblock.resize(3*ignorePastIndex);
	// forceGradientStaticBlock.resize(3*ignorePastIndex, 3*ignorePastIndex);

	VectorXd f, x;
	f.resize(3*TV.rows());
	f.setZero();
	x.resize(3*TV.rows());
	x.setZero();
	setTVtoX(x, TV);

	int NEWTON_MAX = 300, k=0;
	for(k=0; k<NEWTON_MAX; k++){
			cout<<"		Newton Iter "<<k<<endl;
			xToTV(x, TV);
			calculateForceGradient(TV, forceGradient);
			calculateElasticForces(f, TV);

			if(k==0){
				DeletionMatrix.setZero();//optimize this by setting it only once
				int r=0, c=0;
				for(int i=0; i<3*ignorePastIndex; ++i)
				{
					DeletionMatrix.coeffRef(r, c) = 1;
					++r; ++c;
				}
				//because r, c went 1 step too far in the last iteration ^ for the next part to work
				--r;--c;
				for(int i =0; i<moveVertices.size(); ++i)
				{
					// if(abs(staticSolveDirection)==0)
					// {
					// 	++c;
					// 	DeletionMatrix.coeffRef(r, c) = 1;
					// 	++r;++c;
					// 	DeletionMatrix.coeffRef(r, c) = 1;
					// 	++r;++c;
					// }
					// if(abs(staticSolveDirection)==1)
					// {
					// 	++r; ++c;
					// 	DeletionMatrix.coeffRef(r, c) = 1;
					// 	++c;
					// 	DeletionMatrix.coeffRef(r, c) = 1;
					// 	++r;++c;
					// }
					if(abs(staticSolveDirection)==2)
					{
						++r; ++c;
						DeletionMatrix.coeffRef(r, c) = 1;
						++r;++c;
						DeletionMatrix.coeffRef(r, c) = 1;
						++c;
					}
				}
			}

			//Block forceGrad and f to exclude the fixed verts
			forceGradientStaticBlock = DeletionMatrix*forceGradient.block(0,0, 3*ignorePastIndex + 3*moveVertices.size(), 3*ignorePastIndex + 3*moveVertices.size())*DeletionMatrix.transpose();

			fblock.setZero();
			fblocktail.setZero();
			fblock.segment(0, 3*ignorePastIndex) = f.head(3*ignorePastIndex);
			for(int i=0; i<moveVertices.size(); i++){
				// if(abs(staticSolveDirection) ==0)
				// if(abs(staticSolveDirection) ==1)
				if(abs(staticSolveDirection) ==2)
				{
					fblocktail(2*i) = f(3*ignorePastIndex + 3*i);
					fblocktail(2*i + 1) = f(3*ignorePastIndex + 3*i + 1);
				}
			}
			fblock.segment(3*ignorePastIndex, 2*moveVertices.size()) = fblocktail;

			//Sparse QR
			SPQR<SparseMatrix<double>> solver;
			// SparseQR<SparseMatrix<double>> solver;
			solver.compute(forceGradientStaticBlock);
			//------------- Conj Grad------------------
			//  ConjugateGradient<SparseMatrix<double>> solver;
			//  solver.compute(forceGradient);
			//  if(solver.info() == Eigen::NumericalIssue){
			//  	cout<<"ConjugateGradient numerical issue"<<endl;
			//  	exit(0);
			//  }
			// //-----------------------------------------
			VectorXd deltaX = -1*solver.solve(fblock);

			x.segment(0, 3*ignorePastIndex) += deltaX.head(3*ignorePastIndex);
			for(int i=0; i<moveVertices.size(); i++)
			{
				// if(abs(staticSolveDirection)==0)
				// {
				// 	x(3*ignorePastIndex + 2*i+1) += deltaX(3*ignorePastIndex +i);
				// 	x(3*ignorePastIndex + 2*i+2) += deltaX(3*ignorePastIndex +i +1);
				// }
				// if(abs(staticSolveDirection)==1)
				// {
				// 	x(3*ignorePastIndex + 2*i) += deltaX(3*ignorePastIndex +i);
				// 	x(3*ignorePastIndex + 2*i+2) += deltaX(3*ignorePastIndex +i +1);
				// }
				if(abs(staticSolveDirection)==2)
				{
					x(3*ignorePastIndex + 3*i) += deltaX(3*ignorePastIndex + 2*i);
					x(3*ignorePastIndex + 3*i+1) += deltaX(3*ignorePastIndex + 2*i +1);
				}
			}


			if(x != x){
				cout<<"NAN"<<endl;
				exit(0);
			}

			cout<<"fblock/size"<<endl;
			cout<<fblock.squaredNorm()/fblock.size()<<endl;
			if (fblock.squaredNorm()/fblock.size() < 1){
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

void Simulation::applyStaticForces(MatrixXd& TV, MatrixXi& TT, MatrixXd& B, VectorXd& fixed_forces, vector<int>& moveVertices, vector<int>& fixVertices){
	cout<<"***APPLYING STATIC FORCES****"<<endl;

	int staticSolveSteps = 1;
	int ignorePastIndex = TV.rows() - fixVertices.size();

	while(staticSolveSteps<11){
		VectorXd f = fixed_forces*1e-1*staticSolveSteps;
		staticSolveNewtonsForces(TV, TT, B, f, moveVertices, ignorePastIndex, staticSolveSteps);
		staticSolveSteps += 1;
		// printObj(OUTPUT_SAVED_PATH"Boba/StaticSolve/", staticSolveSteps, TV, TT, B);
	}
	//MAX
	double maxd = 100.0;
	for(int i=0; i<moveVertices.size(); i++){
		if(TV.row(moveVertices[i])(2) < maxd){
			maxd = TV.row(moveVertices[i])(2);
		}
	}
	cout<<"Max D"<<endl;
	cout<<maxd<<endl;

}

void Simulation::staticSolveNewtonsForces(MatrixXd& TV, MatrixXi& TT, MatrixXd& B, VectorXd& fixed_forces, vector<int>& moveVertices, int ignorePastIndex, int step){
	cout<<"------------Static Solve Newtons FORCES - Iteration"<< step<<"--------------"<<endl;

	//Newtons method static solve for minimum Strain E
	SparseMatrix<double> forceGradient;
	forceGradient.resize(3*TV.rows(), 3*TV.rows());
	SparseMatrix<double> forceGradientStaticBlock;
	forceGradientStaticBlock.resize(3*ignorePastIndex, 3*ignorePastIndex);
	VectorXd f, x, f_prev;
	f.resize(3*TV.rows());
	f.setZero();
	x.resize(3*TV.rows());
	x.setZero();
	setTVtoX(x, TV);
	int NEWTON_MAX = 100, k=0;

	for(k=0; k<NEWTON_MAX; k++){
		xToTV(x, TV);

		calculateForceGradient(TV, forceGradient);
		calculateElasticForces(f, TV);
		//--------------
		for(int i=0; i<fixed_forces.rows(); i++){
			if(fabs(fixed_forces(i))>0.00001){
				f(i) += fixed_forces(i);
			}
		}

		//Block forceGrad and f to exclude the fixed verts
		forceGradientStaticBlock = forceGradient.block(0,0, 3*(ignorePastIndex), 3*ignorePastIndex);
		VectorXd fblock = f.head(ignorePastIndex*3);



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


		x.segment(0,3*(ignorePastIndex))+=deltaX;
		cout<<"		Newton Iter "<<k<<endl;
		if(x != x){
			cout<<"NAN"<<endl;
			exit(0);
		}

		cout<<"fblock"<<endl;
		cout<<fblock.squaredNorm()/fblock.size()<<endl;
		if (fblock.squaredNorm()/fblock.size() < 0.0001){
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
	ifstream forceInputFile (HOME_SAVED_PATH "shared/"+objectName+".txt");
	cout<<HOME_SAVED_PATH "shared/"+objectName+".txt"<<endl;
	if(forceInputFile.is_open()){
		string line;
		int index =0;
		int fixedIndex=0;
		int moveIndex =0;
		force.setZero();
		while(getline(forceInputFile, line)){
			istringstream iss(line);
			double fx, fy, fz;
			int fixedOrNot; //1 is fixed, 0 not fixed
			if(!(iss >> fx >> fy >> fz >> fixedOrNot)){break;}
			if(fabs(fx + fy + fz)>0){
				temp.push_back(index);
				// cout<<index-fixedIndex<<endl;
				force(3*(index)) = fx;
				force(3*(index)+1) = fy;
				force(3*(index)+2) = fz;
				moveIndex++;
			}

			if(fixedOrNot == 1){
				fixVertices.push_back(index);
				// cout<<"fix "<<index<<endl;
				fixedIndex++;
			}
			index+=1;
		}
		for(int i=0; i<temp.size(); i++){
			moveVertices.push_back(temp[i]);
		}
	}else{
		cout<<"Check yo self: Force input error, file not found"<<endl;
	}
	cout<<"Fixed Verts"<<endl;
	cout<<fixVertices.size()<<endl;
	cout<<"Moving Verts"<<endl;
	cout<<moveVertices.size()<<endl;
	// cout<<"Forces here"<<endl;
	// cout<<force<<endl;
	// exit(0);
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
	cout<<printToHere<<numberOfPrints<< endl;

	// double hausdorffDist;
	// if(integrator->simTime > 1){
	// 	igl::hausdorff(V, F, V, F, hausdorffDist);
	// }
	// V = V_temp;
	// F = F_temp;
	// cout<<"**** Hausdorff Between this and prev it"<<endl;
	// cout<< hausdorffDist<<endl;
	igl::writeOFF(printToHere + to_string(numberOfPrints)+".off", V_temp, F_temp);

	return;
}
