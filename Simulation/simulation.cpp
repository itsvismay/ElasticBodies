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
	cout<<tetgen_code<<endl;
	VectorXd force;
	force.resize(3*TV.rows());
	force.setZero();
	TV_k = TV;
	cout<<"TV.rows()"<<endl;
	cout<<TV.rows()<<endl;
	setInitPosition(force, fixVertices);
	//TEMPORARY CODE REMOVE BEFORE DOING ANYTHIGN ELSE
	// for(int i=0; i< TV.rows(); i++){
	// 	if(TV.row(i)(1)<-40)
	// 		fixVertices.push_back(i);
	// }

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
		
		}
		
		integrator->initializeIntegrator(deltaT, M, newTV, newTT);
		this->external_force = new_force;
		integrator->fixVertices(newfixIndices);
		// int ignorePastIndex = TV.rows() - newfixIndices.size();
		// staticSolveNewtonsForces(newTV, newTT, B, new_force, ignorePastIndex);



	}else{
		igl::barycenter(TV, TT, B);
		M.initializeMesh(TT, TV, youngs, poissons);
		integrator->initializeIntegrator(deltaT, M, TV, TT);
		this->external_force = force;
		integrator->fixVertices(fixVertices);
	}

	sB = &B;
	return 1;
}

void Simulation::applyExternalForces(){
	this->external_force.setZero();
}

void Simulation::headless(){
	int printcount =0;

	integrator->external_f = this->external_force;
	printDesigns(printcount, integrator->simTime);
	while(integrator->simTime < iters){
		integrator->render(this->external_force);
		if(integrator->simTime%10==0){
			printDesigns(printcount, integrator->simTime);
			// printOptimizationOutput();
			printcount += 1;
		}
		igl::writeMESH(to_string(integrator->simTime)+"test1.mesh", integrator->TV, integrator->TT, TF);
	}

}

void Simulation::printDesigns(int printcount, int simTime){
	string saveTestsHere = OUTPUT_SAVED_PATH"TestsResults/SolverTests/"+solver+"/"+to_string(integrator->TT.rows())+"tets@"+tetgen_code+"@"+objectName+"/"+to_string(integrator->h)+"-test1/";
	printObj(saveTestsHere, printcount, integrator->TV, integrator->TT, *sB);

}

void Simulation::printOptimizationOutput(){
	double disp =0;
	for(int i=0; i<this->putForceOnTheseVerts.rows(); i++){
		if (integrator->TV.row(this->putForceOnTheseVerts(i))(2) < disp)
			disp = integrator->TV.row(this->putForceOnTheseVerts(i))(2);
	}
	if(disp < maxDisp){
		maxDisp = disp;
	}
	optimizationFile<<integrator->simTime <<maxDisp<<endl;
	cout<<maxDisp<<"\n";
}

bool Simulation::render(){
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

void Simulation::staticSolveNewtonsForces(MatrixXd& TV, MatrixXi& TT, MatrixXd& B, VectorXd& fixed_forces, int ignorePastIndex){
	cout<<"------------I am here----------------"<<endl;
	// cout<<ignorePastIndex<<endl;
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
		Vector4i indices = M.tets[i].verticesIndex;
		// MatrixXd F_tet = M.tets[i].computeElasticForces(TV, 1);
		// f.segment<3>(3*indices(0)) += F_tet.col(0);
		// f.segment<3>(3*indices(1)) += F_tet.col(1);
		// f.segment<3>(3*indices(2)) += F_tet.col(2);
		// f.segment<3>(3*indices(3)) += F_tet.col(3);
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

void Simulation::setInitPosition(VectorXd& force, vector<int>& fixVertices){
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
				temp.push_back(index - fixedIndex);
				// cout<<"move"<<endl;
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
	cout<<"Fixed Verts"<<endl;
	cout<<fixVertices.size()<<endl;
	cout<<"Moving Verts"<<endl;
	cout<<putForceOnTheseVerts.size()<<endl;

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
		
		cout<<"Fblock/size"<<endl;
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
	exit(0);				
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
	igl::writeOBJ(printToHere + to_string(numberOfPrints)+".obj", V_temp, F_temp);
	
	return;
}

