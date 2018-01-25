#include "IntegratorsAbstract.h"

typedef Eigen::Triplet<double> Trip;
typedef Matrix<double, 12, 1> Vector12d;
//TODO: Optimize this using hashing
bool IntegratorAbstract::isFixed(int vert){

	for(unsigned int j=0; j<fixedVerts.size(); j++){
		if(vert == fixedVerts[j]){
			return true;
		}
	}
	return false;
}

void IntegratorAbstract::printInfo(){
	////////////////////////////////////

	double TotalEnergy = 0;
	double gravityE =0;
	double kineticE =0;
	double strainE = 0.0;
	for(int i=0; i<vertsNum; i++){
		if(!isFixed(i)){
			int k=3*i;
			// gravityE +=  massVector(k)*-1*gravity*(x_old(k));
			kineticE += 0.5*massVector(k)*v_old(k)*v_old(k);

			k++;
			kineticE += 0.5*massVector(k)*v_old(k)*v_old(k);

			k++;
			// gravityE +=  massVector(k)*-1*gravity*(x_old(k));
			gravityE +=  massVector(k)*gravity*(x_old(k));
			kineticE += 0.5*massVector(k)*v_old(k)*v_old(k);
		}
	}

	for(int i=0; i< M.tets.size(); i++){
		strainE += M.tets[i].energy;
	}
	// cout<<totalVol<<endl;
	// optimizationFile<<totalVol<<endl;

	TotalEnergy+= gravityE + kineticE + strainE;
	// cout<<endl<<"Grav E"<<endl;
	// cout<<strainE<<endl;
	cout<<"Tot E"<<endl;
	cout<<TotalEnergy<<endl;
	cout<<"Strain E"<<endl;
	cout<<strainE<<endl;
	cout<<"Gravity E"<<endl;
	cout<<gravityE<<endl;
	cout<<"Kinetic E"<<endl;
	cout<<kineticE<<endl;
	cout<<endl;


	energyFile<<simTime<<", "<<TotalEnergy<<"\n";
	strainEnergyFile<<simTime<<", "<<strainE<<"\n";
	kineticEnergyFile<<simTime<<", "<<kineticE<<"\n";
	gravityEnergyFile<<simTime<<", "<<gravityE<<"\n";

	////////////////////////////////////

	// //////Momentum Code////////////
	// if(momentumFile.is_open()){
	// 	double xp=0;
	// 	double yp=0;
	// 	double zp=0;
	// 	cout<<"-----------"<<endl;
	// 	for(int i=0; i<v_old.rows(); ){
	// 		xp += v_old(i)*massVector(i);
	// 		cout<<"x"<<endl;
	// 		cout<<v_old(i)<<endl;
	// 		i++;
	// 		yp += v_old(i)*massVector(i);
	// 		cout<<"y"<<endl;
	// 		cout<<v_old(i)<<endl;
	// 		i++;
	// 		zp += v_old(i)*massVector(i);
	// 		cout<<"z"<<endl;
	// 		cout<<v_old(i)<<endl;
	// 		i++;
	// 	}
	// 	cout<<"- - - - -- "<<endl;
	// 	momentumFile<<simTime<<","<<xp<<","<<yp<<","<<zp<<"\n";
	// }else{
	// 	cout<<"no open file"<<endl;
	// }
	/////////////////
}

void IntegratorAbstract::initializeIntegrator(double ph, SolidMesh& pM, MatrixXd& pTV, MatrixXi& pTT){
	//Constants
	vertsNum = pTV.rows();

	h = ph;
	M = pM;
	TV = pTV;
	TT = pTT;

	initVectors();
	initMassMatrices();
	createXFromTet();
}

void IntegratorAbstract::initVectors(){
	x_old.resize(3*vertsNum);
	v_old.resize(3*vertsNum);
	f.resize(3*vertsNum);
	massVector.resize(3*vertsNum);
	forceGradient.resize(3*vertsNum, 3*vertsNum);
	CholeskyAnalyze.resize(3*vertsNum, 3*vertsNum);

	x_old.setZero();
	v_old.setZero();
	f.setZero();
	massVector.setZero();

}

void IntegratorAbstract::analyzeCholeskySetup(){
	if(solver.compare("newton") != 0)
		return; //only do cholesky solve on NM

	CholeskyAnalyze.setZero();
	int ignorePastIndex = TV.rows() - fixedVerts.size();
	cout<<"Ignore past"<<endl;
	cout<<ignorePastIndex<<endl;
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
			kj = j%3;
			//row in order for dfxi/dxi ..dfxi/dzl
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[0], 1));
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[0]+1, 1));
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[0]+2, 1));

			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[1], 1));
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[1]+1, 1));
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[1]+2, 1));

			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[2], 1));
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[2]+1, 1));
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[2]+2, 1));

			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[3], 1));
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[3]+1, 1));
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[3]+2, 1));
		}
	}
	CholeskyAnalyze.setFromTriplets(triplets1.begin(), triplets1.end());

	SparseMatrix<double> CholeskyAnalyzeBlock = CholeskyAnalyze.block(0,0, 3*(ignorePastIndex), 3*ignorePastIndex);
	// cout<<"analyzing pattern******"<<endl;
	// cout<<CholeskyAnalyze<<endl;
	// cout<<"analyzing pattern******"<<endl;
	// cout<<CholeskyAnalyzeStaticBlock<<endl;
	cout<<"analyzing pattern******"<<endl;
	llt_solver.analyzePattern(CholeskyAnalyzeBlock);
}

void IntegratorAbstract::initMassMatrices(){
	InvMass.resize(3*vertsNum, 3*vertsNum);
	RegMass.resize(3*vertsNum, 3*vertsNum);

	for(unsigned int i=0; i<M.tets.size(); i++){
		double vol = (M.tets[i].undeformedVol/4)*1.05e3; //UNITS: kg/cubic m
		Vector4i indices = M.tets[i].verticesIndex;

		massVector(3*indices(0)) += vol;
		massVector(3*indices(0)+1) += vol;
		massVector(3*indices(0)+2) += vol;

		massVector(3*indices(1)) += vol;
		massVector(3*indices(1)+1) += vol;
		massVector(3*indices(1)+2) += vol;

		massVector(3*indices(2)) += vol;
		massVector(3*indices(2)+1) += vol;
		massVector(3*indices(2)+2) += vol;

		massVector(3*indices(3)) += vol;
		massVector(3*indices(3)+1) += vol;
		massVector(3*indices(3)+2) += vol;
	}
	vector<double>tempForMedian;
	for(int i=0; i<3*vertsNum; i++){
		InvMass.coeffRef(i,i) = 1/massVector(i);
		RegMass.coeffRef(i,i) = massVector(i);
		tempForMedian.push_back(massVector(i));
	}
	sort(tempForMedian.begin(), tempForMedian.end());
	if(tempForMedian.size()%2 == 0){
		this->convergence_scaling_paramter = 0.5*(tempForMedian[tempForMedian.size()/2-1]+tempForMedian[tempForMedian.size()/2]);
	}else{
		this->convergence_scaling_paramter = tempForMedian[tempForMedian.size()/2];
	}
	cout<<"MEDIAN"<<endl;
	cout<<this->convergence_scaling_paramter<<endl;
	// cout<<"Mass Vector"<<endl;
	// cout<<"Mass Vector"<<endl;
	// cout<<massVector<<endl;
	// cout<<"INV Mass"<<endl;
	// cout<<InvMass<<endl;
	// cout<<"Reg Mass"<<endl;
	// cout<<RegMass<<endl;
	// cout<<"MASS ENTRY"<<endl;
	// cout<<massVector(4)<<endl;
}

void IntegratorAbstract::fixVertices(vector<int> fixMe){
	fixedVerts.insert(fixedVerts.end(), fixMe.begin(), fixMe.end());

	for(int i=0; i<fixMe.size(); i++){
		massVector(3*fixMe[i]) = 1e20;
		massVector(3*fixMe[i]+1) = 1e20;
		massVector(3*fixMe[i]+2) = 1e20;

		InvMass.coeffRef(3*fixMe[i], 3*fixMe[i]) = 0;
		InvMass.coeffRef(3*fixMe[i]+1, 3*fixMe[i]+1) = 0;
		InvMass.coeffRef(3*fixMe[i]+2, 3*fixMe[i]+2) = 0;

		RegMass.coeffRef(3*fixMe[i], 3*fixMe[i]) = 1e20;
		RegMass.coeffRef(3*fixMe[i]+1, 3*fixMe[i]+1) = 1e20;
		RegMass.coeffRef(3*fixMe[i]+2, 3*fixMe[i]+2) = 1e20;

		//set vels to 0
		v_old.segment<3>(3*fixMe[i])*=0;
	}

	//must recreate whenever fixed
	analyzeCholeskySetup();

}

void IntegratorAbstract::moveVertices(vector<int> moveMe){
	double hard_coded_mass = 0000; //mass in grams
	double factor = hard_coded_mass/(3*moveMe.size());
	double totalMovingMass = 0;
	for(int i=0; i<moveMe.size(); i++){
		// massVector(3*moveMe[i]) = factor;
		// massVector(3*moveMe[i]+1) = factor;
		// massVector(3*moveMe[i]+2) = factor;

		InvMass.coeffRef(3*moveMe[i], 3*moveMe[i]) = 0;
		InvMass.coeffRef(3*moveMe[i] + 1, 3*moveMe[i] + 1) = 0;
		InvMass.coeffRef(3*moveMe[i] + 2, 3*moveMe[i] + 2) = 0;

		RegMass.coeffRef(3*moveMe[i], 3*moveMe[i]) = factor;
		RegMass.coeffRef(3*moveMe[i] + 1, 3*moveMe[i] + 1) = factor;
		RegMass.coeffRef(3*moveMe[i] + 2, 3*moveMe[i] + 2) = factor;
		totalMovingMass += massVector(3*moveMe[i]) + massVector(3*moveMe[i]+1) + massVector(3*moveMe[i]+2);
	}
	cout<<"TOTAL MOVING MASS"<<endl;
	cout<<totalMovingMass<<endl;
	// dampingPositionFile<<totalMovingMass<<", "<< "Spring Mass in Grams"<<endl;

}

void IntegratorAbstract::createXFromTet(){
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
