#include "IntegratorsAbstract.h"

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
	double strainE = 0;
	for(int i=0; i<vertsNum; i++){
		if(!isFixed(i)){
			int k=3*i;
			// gravityE +=  massVector(k)*-1*gravity*(x_old(k));
			kineticE += 0.5*massVector(k)*v_old(k)*v_old(k);
			
			k++;
			gravityE +=  massVector(k)*gravity*(x_old(k));
			kineticE += 0.5*massVector(k)*v_old(k)*v_old(k);
			
			k++;
			// gravityE +=  massVector(k)*-1*gravity*(x_old(k));
			kineticE += 0.5*massVector(k)*v_old(k)*v_old(k);
		}		
	}
	
	for(int i=0; i< M.tets.size(); i++){
		strainE += M.tets[i].undeformedVol*M.tets[i].energyDensity;
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
	
	x_old.setZero();
	v_old.setZero();
	f.setZero();
	massVector.setZero();

	v_old(0) =10;
	v_old(1) =100;
	// // v_old(2) =1;
	// v_old(3) =1;
}

void IntegratorAbstract::initMassMatrices(){
	InvMass.resize(3*vertsNum, 3*vertsNum);
	RegMass.resize(3*vertsNum, 3*vertsNum);

	for(unsigned int i=0; i<M.tets.size(); i++){
		double vol = (M.tets[i].undeformedVol/4)*9.7e-7; //UNITS: THIS ACTUALLY DENSITY
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

	for(int i=0; i<3*vertsNum; i++){
		InvMass.coeffRef(i,i) = 1/massVector(i);
		RegMass.coeffRef(i,i) = massVector(i);
	}
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
		massVector(3*fixMe[i]) = 1000000000000;
		massVector(3*fixMe[i]+1) = 1000000000000;
		massVector(3*fixMe[i]+2) = 1000000000000;

		InvMass.coeffRef(3*fixMe[i], 3*fixMe[i]) = 0;
		InvMass.coeffRef(3*fixMe[i]+1, 3*fixMe[i]+1) = 0;
		InvMass.coeffRef(3*fixMe[i]+2, 3*fixMe[i]+2) = 0;

		RegMass.coeffRef(3*fixMe[i], 3*fixMe[i]) = 1000000000000;
		RegMass.coeffRef(3*fixMe[i]+1, 3*fixMe[i]+1) = 1000000000000;
		RegMass.coeffRef(3*fixMe[i]+2, 3*fixMe[i]+2) = 1000000000000;

		//set vels to 0
		v_old.segment<3>(3*fixMe[i])*=0;
	}

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