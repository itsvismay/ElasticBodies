#include <igl/viewer/Viewer.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/barycenter.h>
#include <iostream>
#include <vector>
#include <pthread.h>
#include <fstream>
#include <math.h>

#ifndef TUTORIAL_SHARED_PATH
#define TUTORIAL_SHARED_PATH "../shared"
#endif

using namespace Eigen;
using namespace std;
typedef Eigen::Triplet<double> Trip;
typedef Matrix<double, 12, 1> Vector12d;

// Input polygon
Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd B;
ofstream dataFile;
// Tetrahedralized interior
Eigen::MatrixXd TV;
Eigen::MatrixXi TT;
Eigen::MatrixXi TF;
// double mu = 675450000.0;
// double lambda = 3449700000.0;
double mu = 1000;
double lambda = 1000;
double rayleighCoeff = 10;

//####################################################################################
//####################################################################################
//####################################################################################
//####################################################################################
class Tetrahedron
{
public:
    Matrix3d DeformedShapeMatrix, ReferenceShapeMatrix, InvRefShapeMatrix;
    VectorXi verticesIndex;
    double undeformedVol;
    Tetrahedron(VectorXi k);
    MatrixXd computeElasticForces(MatrixXd& TV, int e);
    void precomputeForces(MatrixXd& TV);
    VectorXd computeForceDifferentials(MatrixXd& TV, Vector12d& dx);
    Matrix3d computeDeltaDs(const Vector12d& dx);
    Matrix3d computeDs(const Vector12d& x);

};

Tetrahedron::Tetrahedron(VectorXi k){
    verticesIndex = k ;
}

void Tetrahedron::precomputeForces(MatrixXd& TV){
	double xi, yi, zi, xj, yj, zj, xk, yk, zk, xl, yl, zl;
	Vector3d ro1 = TV.row(this->verticesIndex(0));
	xi = ro1(0);
    yi = ro1(1);
    zi = ro1(2);

    Vector3d ro2 = TV.row(this->verticesIndex(1));
    xj = ro2(0);
    yj = ro2(1);
    zj = ro2(2);

    Vector3d ro3 = TV.row(this->verticesIndex(2));
    xk = ro3(0);
    yk = ro3(1);
    zk = ro3(2);

    Vector3d ro4 = TV.row(this->verticesIndex(3));
    xl = ro4(0);
    yl = ro4(1);
    zl = ro4(2);

    Matrix3d Dm;
    Dm << (xi - xl), (xj - xl), (xk - xl),
                    (yi - yl), (yj - yl), (yk - yl),
                    (zi - zl), (zj - zl), (zk - zl);
    this->ReferenceShapeMatrix = Dm;
    this->InvRefShapeMatrix = Dm.inverse();
    this->undeformedVol = (1.0/6)*abs(Dm.determinant());
}

Matrix3d Tetrahedron::computeDeltaDs(const Vector12d& dx){
	double xio, yio, zio, xjo, yjo, zjo, xko, yko, zko, xlo, ylo, zlo;
    //Vector3d ro1 = TV.row(this->verticesIndex(0));
	xio = dx(0);
    yio = dx(1);
    zio = dx(2);

    //Vector3d ro2 = TV.row(this->verticesIndex(1));
    xjo = dx(3);
    yjo = dx(4);
    zjo = dx(5);

    //Vector3d ro3 = TV.row(this->verticesIndex(2));
    xko = dx(6);
    yko = dx(7);
    zko = dx(8);

    //Vector3d ro4 = TV.row(this->verticesIndex(3));
    xlo = dx(9);
    ylo = dx(10);
    zlo = dx(11);

    Matrix3d dDs;
    dDs << (xio - xlo), (xjo - xlo), (xko - xlo),
            (yio - ylo), (yjo - ylo), (yko - ylo),
            (zio - zlo), (zjo - zlo), (zko - zlo);
     return dDs;
}

Matrix3d Tetrahedron::computeDs(const Vector12d& x){
	double xi, yi, zi, xj, yj, zj, xk, yk, zk, xl, yl, zl;
    //Vector3d ro1 = TV.row(this->verticesIndex(0));
	xi = x(0);
    yi = x(1);
    zi = x(2);

    //Vector3d ro2 = x(0)(this->verticesIndex(1));
    xj = x(3);
    yj = x(4);
    zj = x(5);

    //Vector3d ro3 = x(0)(this->verticesIndex(2));
    xk = x(6);
    yk = x(7);
    zk = x(8);

    //Vector3d ro4 = x(0)(this->verticesIndex(3));
    xl = x(9);
    yl = x(10);
    zl = x(11);

    Matrix3d Ds;
    Ds << (xi - xl), (xj - xl), (xk - xl),
            (yi - yl), (yj - yl), (yk - yl),
            (zi - zl), (zj - zl), (zk - zl);

    return Ds;
}

VectorXd Tetrahedron::computeForceDifferentials(MatrixXd& TV, Vector12d& dx){
	Vector12d x;
    x.segment<3>(0) = TV.row(this->verticesIndex(0));
    x.segment<3>(3) = TV.row(this->verticesIndex(1));
    x.segment<3>(6) = TV.row(this->verticesIndex(2));
    x.segment<3>(9) = TV.row(this->verticesIndex(3));


    Matrix3d Ds = computeDs(x);
    Matrix3d dDs = computeDeltaDs(dx);



    ////////////////////TEST dDs correctness///////////
    // double epsilon = 0.000001;
	//f(v+[e,0,0,0...]) - f(v) / e = df/dx
	// cout<<"test dDs"<<endl;
	// Matrix3d leftdDs = computeDeltaDs(dx*epsilon);
	// Matrix3d rightdDs = computeDs(x + dx*epsilon) - Ds;
	// cout<<leftdDs<<endl<<endl;
	/////////////////////////////////////////////////

    Matrix3d F = Ds*this->InvRefShapeMatrix;
    Matrix3d dF = dDs*this->InvRefShapeMatrix;

    ////////////////////TEST dF correctness///////////
	// cout<<"test dF"<<endl;
	// Matrix3d leftdF = (leftdDs)*this->InvRefShapeMatrix;
	// Matrix3d rightF = computeDs(x+dx*epsilon)*this->InvRefShapeMatrix;
	// Matrix3d rightdF =  rightF - Ds*this->InvRefShapeMatrix;
	// cout<<leftdF - rightdF<<endl<<endl;
	// //////////////////////////////////////////////////////
    
    //Neohookean
    Matrix3d P = mu*(F - ((F.inverse()).transpose())) + lambda*log(F.determinant())*((F.inverse()).transpose());
    Matrix3d dP = mu*dF + (mu - lambda*log(F.determinant()))*((F.inverse()).transpose())*dF.transpose()*((F.inverse()).transpose()) + lambda*(F.inverse()*dF).trace()*((F.inverse()).transpose());
    
    ////////////////////TEST dP correctness///////////
    // cout<<"test dP"<<endl;
    // Matrix3d leftdP = mu*leftdF + (mu - lambda*log(rightF.determinant()))*((rightF.inverse()).transpose())*leftdF.transpose()*((rightF.inverse()).transpose()) + lambda*(rightF.inverse()*leftdF).trace()*((rightF.inverse()).transpose());
	// Matrix3d rightP1 = mu*(rightF - ((rightF.inverse()).transpose())) + lambda*log(rightF.determinant())*((rightF.inverse()).transpose());
	// Matrix3d rightP2 = mu*(F - ((F.inverse()).transpose())) + lambda*log(F.determinant())*((F.inverse()).transpose());
	// Matrix3d rightdP = rightP1 - rightP2;
	//cout<< leftdP - (rightP1 - rightP2)<<endl<<endl;
	//////////////////////////////////////////////////////

    Matrix3d dH = -1*this->undeformedVol*dP*((this->InvRefShapeMatrix).transpose());
    
    ////////////////////TEST dH correctness///////////
    // cout<<"test dH"<<endl;
    // Matrix3d leftdH = -1*this->undeformedVol*leftdP*((this->InvRefShapeMatrix).transpose());
	// Matrix3d rightH1 = -1*this->undeformedVol*rightP1*((this->InvRefShapeMatrix).transpose());
	// Matrix3d rightH2 = -1*this->undeformedVol*P*((this->InvRefShapeMatrix).transpose());
	// Matrix3d rightdH = rightH1 - rightH2;
	// cout<< leftdH - rightdH<<endl<<endl;
	/////////

    MatrixXd dForces(3,4);
    dForces.col(0) = dH.col(0);
    dForces.col(1) = dH.col(1);
    dForces.col(2) = dH.col(2);
    dForces.col(3) = -1*dH.col(0) - dH.col(1) - dH.col(2);
    // cout<<endl<<"dForces matrix"<<endl;
    // cout<<dForces<<endl<<endl;;
    // cout<<"Mapped vector"<<endl;
    // cout<<Map<VectorXd>(dForces.data(), dForces.cols()*dForces.rows())<<endl;
    return Map<VectorXd>(dForces.data(), dForces.cols()*dForces.rows());
}

MatrixXd Tetrahedron::computeElasticForces(MatrixXd &TV, int e){
    Vector12d x;
    x.segment<3>(0) = TV.row(this->verticesIndex(0));
    x.segment<3>(3) = TV.row(this->verticesIndex(1));
    x.segment<3>(6) = TV.row(this->verticesIndex(2));
    x.segment<3>(9) = TV.row(this->verticesIndex(3));

    Matrix3d Ds = computeDs(x);

    Matrix3d F = Ds*this->InvRefShapeMatrix;

    //TODO: Spring Constant value
    Matrix3d E = 0.5*((F.transpose()*F) - MatrixXd::Identity(3,3));

    //SVK
    //Matrix3d P = F*(2*mu*E + lambda*E.trace()*MatrixXd::Identity(3,3));//piola kirchoff	
	
    //Neo
	Matrix3d P = mu*(F - ((F.inverse()).transpose())) + lambda*log(F.determinant())*((F.inverse()).transpose());
    
    
    Matrix3d H = -1*this->undeformedVol*P*((this->InvRefShapeMatrix).transpose());
    Matrix<double, 3, 4> Forces;
    Forces.col(0) = H.col(0);
    Forces.col(1) = H.col(1);
    Forces.col(2) = H.col(2);
    Forces.col(3) = -1*H.col(0) - H.col(1) - H.col(2);
    return Forces;
}
//####################################################################################
//####################################################################################
//####################################################################################
//####################################################################################

class SolidMesh
{
    public:
        vector<Tetrahedron> tets;
        SolidMesh();
        SolidMesh(MatrixXi& TT, MatrixXd& TV);
        void initializeMesh(MatrixXi& TT, MatrixXd& TV);
};

SolidMesh::SolidMesh(void){}

SolidMesh::SolidMesh(MatrixXi& TT, MatrixXd& TV){
    for(int i=0; i<TT.rows(); i++){
    	int i1 = TT.row(i)[0];
    	int i2 = TT.row(i)[1];
    	int i3 = TT.row(i)[2];
    	int i4 = TT.row(i)[3];

    	Tetrahedron t(TT.row(i));
    	t.precomputeForces(TV);
    	this->tets.push_back(t);
    }
}

void SolidMesh::initializeMesh(MatrixXi& TT, MatrixXd& TV){
	for(int i=0; i<TT.rows(); i++){
    	//based on Tet indexes, get the vertices of the tet from TV
    	Tetrahedron t(TT.row(i));
    	t.precomputeForces(TV);
    	this->tets.push_back(t);
    }
}

//####################################################################################
//####################################################################################
//####################################################################################
//####################################################################################
class Simulation{

public:
	int t=0;
	double timestep;
	SparseMatrix<double> InvMass;
	SolidMesh M;
	int vertices;
	VectorXd x_old, v_old, vertex_masses;
	VectorXd f;
	vector<int> mapV2TV;
	MatrixXd TV;
	SparseMatrix<double> ZeroMatrix;


	Simulation(void);
	void renderExplicit();
	void renderImplicit();
	void createInvMassMatrix();
	void createXFromTet();
	void createForceVector();
	void createTVFromTet();
	void initializeSimulation(double deltaT, MatrixXi& TT_One, MatrixXd& TV_One, vector<int>& map);
	
	void setXIntoTV(VectorXd& x_new);

	void calculateGravity();
	void calculateElasticForce();

	void fixVertices(int fixed);

	SparseMatrix<double> ImplicitCalculateElasticForceGradient(MatrixXd& TVk);
	VectorXd ImplicitCalculateForces(MatrixXd& TVk, SparseMatrix<double>& forceGradient, VectorXd& v_k);
	MatrixXd ImplicitXtoTV(VectorXd& x_tv);
	void ImplicitXfromTV(VectorXd& x_n1, MatrixXd& TVk);
};


Simulation::Simulation(void){}


void Simulation::initializeSimulation(double deltaT, MatrixXi& TT_One, MatrixXd& TV_One, vector<int>& map){
	timestep = deltaT;
	mapV2TV = map;
	TV = TV_One;
	vertices = TV_One.rows();
	ZeroMatrix.resize(3*vertices, 3*vertices);
	ZeroMatrix.setZero();

	x_old.resize(3*vertices);
	x_old.setZero();
	v_old.resize(3*vertices);
	v_old.setZero();

	f.resize(3*vertices);
	f.setZero();
	
	M.initializeMesh(TT_One, TV);	
	createInvMassMatrix(); //creates InvMass
	createXFromTet(); //creates x_old
}

void Simulation::createInvMassMatrix(){
	vertex_masses.resize(3*vertices);
	vertex_masses.setZero();

	InvMass.resize(3*vertices, 3*vertices);
	InvMass.setZero();

	for(int i=0; i<M.tets.size(); i++){
		double vol = (M.tets[i].undeformedVol/4);// TODO: add density 1/Mass for inv matrix,  assume density is 1, so volume=mass
		Vector4i indices = M.tets[i].verticesIndex;

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

	// cout<<vertex_masses<<endl;

}

void Simulation::fixVertices(int fixed){
	
	vertex_masses(3*fixed) = 1000000000000;
	vertex_masses(3*fixed+1) = 1000000000000;
	vertex_masses(3*fixed+2) = 1000000000000;

	InvMass.coeffRef(3*fixed, 3*fixed) = 0;
	InvMass.coeffRef(3*fixed+1, 3*fixed+1) = 0;
	InvMass.coeffRef(3*fixed+2, 3*fixed+2) = 0;
}

void Simulation::createXFromTet(){
	x_old.setZero();
	for(int i = 0; i < M.tets.size(); i++){
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
	for(int i=0; i < M.tets.size(); i++){
		Vector4i indices = M.tets[i].verticesIndex;
		TV.row(indices(0)) = Vector3d(x_new(3*indices(0)), x_new(3*indices(0)+1), x_new(3*indices(0) +2));
		TV.row(indices(1)) = Vector3d(x_new(3*indices(1)), x_new(3*indices(1)+1), x_new(3*indices(1) +2));
		TV.row(indices(2)) = Vector3d(x_new(3*indices(2)), x_new(3*indices(2)+1), x_new(3*indices(2) +2));
		TV.row(indices(3)) = Vector3d(x_new(3*indices(3)), x_new(3*indices(3)+1), x_new(3*indices(3) +2)); 
	}
}

void Simulation::createForceVector(){
	f.setZero();
	calculateGravity();
	calculateElasticForce();
	// cout<<endl<<"Force sq norm"<<endl;
	//cout<<f.squaredNorm()<<endl;
}

void Simulation::calculateGravity(){
	for(int i=0; i<M.tets.size(); i++){
		double vertex_mass = M.tets[i].undeformedVol/4;//assume const density 1
		Vector4i indices = M.tets[i].verticesIndex;

		f(3*indices(0)+1) += vertex_mass*-9.81;

		f(3*indices(1)+1) += vertex_mass*-9.81; 

		f(3*indices(2)+1) += vertex_mass*-9.81;

		f(3*indices(3)+1) += vertex_mass*-9.81;
	}
}

void Simulation::calculateElasticForce(){
	for(int i=0; i<M.tets.size(); i++){
		Vector4i indices = M.tets[i].verticesIndex;
		MatrixXd F_tet = M.tets[i].computeElasticForces(TV, t%2);

		f(3*indices(0)) += F_tet.col(0)[0];
		f(3*indices(0)+1) += F_tet.col(0)[1];
		f(3*indices(0)+2) += F_tet.col(0)[2];

		f(3*indices(1)) += F_tet.col(1)[0];
		f(3*indices(1)+1) += F_tet.col(1)[1]; 
		f(3*indices(1)+2) += F_tet.col(1)[2];

		f(3*indices(2)) += F_tet.col(2)[0];
		f(3*indices(2)+1) += F_tet.col(2)[1];
		f(3*indices(2)+2) += F_tet.col(2)[2];

		f(3*indices(3)) += F_tet.col(3)[0];
		f(3*indices(3)+1) += F_tet.col(3)[1];
		f(3*indices(3)+2) += F_tet.col(3)[2];
	}
}

MatrixXd Simulation::ImplicitXtoTV(VectorXd& x_tv){
	MatrixXd TVk(vertices, 3);
	TVk.setZero();
	for(int i=0; i < M.tets.size(); i++){
		Vector4i indices = M.tets[i].verticesIndex;
		TVk.row(indices(0)) = Vector3d(x_tv(3*indices(0)), x_tv(3*indices(0)+1), x_tv(3*indices(0) +2));
		TVk.row(indices(1)) = Vector3d(x_tv(3*indices(1)), x_tv(3*indices(1)+1), x_tv(3*indices(1) +2));
		TVk.row(indices(2)) = Vector3d(x_tv(3*indices(2)), x_tv(3*indices(2)+1), x_tv(3*indices(2) +2));
		TVk.row(indices(3)) = Vector3d(x_tv(3*indices(3)), x_tv(3*indices(3)+1), x_tv(3*indices(3) +2)); 
	}
	return TVk;
}

VectorXd Simulation::ImplicitCalculateForces( MatrixXd& TVk, SparseMatrix<double>& forceGradient, VectorXd& v_k){
	//gravity
	VectorXd forces(3*vertices);
	forces.setZero();

	for(int i=0; i<M.tets.size(); i++){
		double vertex_mass = M.tets[i].undeformedVol/4;//assume const density 1
		Vector4i indices = M.tets[i].verticesIndex;

		forces(3*indices(0)+1) += vertex_mass*-9.81;

		forces(3*indices(1)+1) += vertex_mass*-9.81; 

		forces(3*indices(2)+1) += vertex_mass*-9.81;

		forces(3*indices(3)+1) += vertex_mass*-9.81;
	}

	//elastic
	for(int i=0; i<M.tets.size(); i++){
		Vector4i indices = M.tets[i].verticesIndex;
		MatrixXd F_tet = M.tets[i].computeElasticForces(TVk, t%2);
		forces(3*indices(0)) += F_tet.col(0)[0];
		forces(3*indices(0)+1) += F_tet.col(0)[1];
		forces(3*indices(0)+2) += F_tet.col(0)[2];

		forces(3*indices(1)) += F_tet.col(1)[0];
		forces(3*indices(1)+1) += F_tet.col(1)[1]; 
		forces(3*indices(1)+2) += F_tet.col(1)[2];

		forces(3*indices(2)) += F_tet.col(2)[0];
		forces(3*indices(2)+1) += F_tet.col(2)[1];
		forces(3*indices(2)+2) += F_tet.col(2)[2];

		forces(3*indices(3)) += F_tet.col(3)[0];
		forces(3*indices(3)+1) += F_tet.col(3)[1];
		forces(3*indices(3)+2) += F_tet.col(3)[2];
	}

	//damping
	// forces+= rayleighCoeff*forceGradient*v_k;
	cout<<forces.squaredNorm()<<endl;
	return forces;
}

SparseMatrix<double> Simulation::ImplicitCalculateElasticForceGradient(MatrixXd& TVk){
	SparseMatrix<double> global_gradForces(3*vertices, 3*vertices);
	global_gradForces.setZero();

	for(int i=0; i<M.tets.size(); i++){
		//Get P(dxn), dx = [1,0, 0...], then [0,1,0,....], and so on... for all 4 vert's x, y, z
		//P is the compute Force Differentials blackbox fxn
		MatrixXd local_gradForces(12, 12);
		local_gradForces.setZero();

		Vector12d dx(12);
		dx.setZero();
		for(int j=0; j<12; j++){
			dx(j) = 1;
			VectorXd df = M.tets[i].computeForceDifferentials(TVk, dx);
			local_gradForces.col(j) = df;// TODO values are really big. Check if correct
			dx(j) = 0; //ASK check is this efficient?
		}


		//TODO optimize this somehow. Create Triplet list from local grad forces
		vector<Trip> triplets;
		triplets.reserve(9*16);//estimation of entries
		Vector4i indices = M.tets[i].verticesIndex;
		for(int ki=0; ki<4; ki++){
			for(int kj=0; kj<4; kj++){
				triplets.push_back(Trip(3*indices[ki], 3*indices[kj], local_gradForces(3*ki, 3*kj))); //dfxi/dxi
				triplets.push_back(Trip(3*indices[ki], 3*indices[kj]+1, local_gradForces(3*ki, 3*kj+1))); //dfxi/dyi
				triplets.push_back(Trip(3*indices[ki], 3*indices[kj]+2, local_gradForces(3*ki, 3*kj+2))); //dfxi/dzi

				triplets.push_back(Trip(3*indices[ki]+1,  3*indices[kj], local_gradForces(3*ki+1, 3*kj))); //dfyi/dxi
				triplets.push_back(Trip(3*indices[ki]+1,  3*indices[kj]+1, local_gradForces(3*ki+1, 3*kj+1))); //dfyi/dyi
				triplets.push_back(Trip(3*indices[ki]+1,  3*indices[kj]+2, local_gradForces(3*ki+1, 3*kj+2))); //dfyi/dzi

				triplets.push_back(Trip(3*indices[ki]+2,  3*indices[kj], local_gradForces(3*ki+2, 3*kj))); //dfzi/dxi
				triplets.push_back(Trip(3*indices[ki]+2,  3*indices[kj]+1, local_gradForces(3*ki+2, 3*kj+1))); //dfzi/dyi
				triplets.push_back(Trip(3*indices[ki]+2,  3*indices[kj]+2, local_gradForces(3*ki+2, 3*kj+2))); //dfzi/dzi
			}
		}
		//Now insert triplets into global matrix
		global_gradForces.setFromTriplets(triplets.begin(), triplets.end());
	}
	// cout<<"global forces"<<endl;
	// SparseMatrix<double> gFT = global_gradForces.transpose();
	// cout<<global_gradForces-gFT<<endl;
	return global_gradForces; //ASK is the -1 correct?
}
void Simulation::ImplicitXfromTV(VectorXd& x_n1, MatrixXd& TVk){
	x_n1.setZero();
	for(int i = 0; i < M.tets.size(); i++){
		Vector4i indices = M.tets[i].verticesIndex;

		x_n1(3*indices(0)) = TVk.row(indices(0))[0];
		x_n1(3*indices(0)+1) = TVk.row(indices(0))[1];
		x_n1(3*indices(0)+2) = TVk.row(indices(0))[2];

		x_n1(3*indices(1)) = TVk.row(indices(1))[0];
		x_n1(3*indices(1)+1) = TVk.row(indices(1))[1];
		x_n1(3*indices(1)+2) = TVk.row(indices(1))[2];

		x_n1(3*indices(2)) = TVk.row(indices(2))[0];
		x_n1(3*indices(2)+1) = TVk.row(indices(2))[1];
		x_n1(3*indices(2)+2) = TVk.row(indices(2))[2];

		x_n1(3*indices(3)) = TVk.row(indices(3))[0];
		x_n1(3*indices(3)+1) = TVk.row(indices(3))[1];
		x_n1(3*indices(3)+2) = TVk.row(indices(3))[2];
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
	MatrixXd TVk(vertices, 3);
	TVk.setZero();
	MatrixXd TVk_prev = TVk;
	VectorXd v_k = v_old;
	VectorXd x_k = x_old+timestep*v_k;
	SparseMatrix<double> Ident = MatrixXd::Identity(3*vertices, 3*vertices).sparseView();
	SparseMatrix<double> grad_g(3*vertices, 3*vertices);

	
	for(int i=0; i<20; i++){
		TVk = ImplicitXtoTV(x_k);//TVk value changed in function
		grad_g.setZero();
		
		SparseMatrix<double> forceGradient = ImplicitCalculateElasticForceGradient(TVk); 
		
		VectorXd g = v_k - (v_old + timestep*InvMass*ImplicitCalculateForces(TVk, forceGradient, v_k));
		grad_g = Ident - timestep*timestep*InvMass*forceGradient;// - ZeroMatrix -timestep*InvMass*rayleighCoeff*forceGradient;//Zero matrix for the Hessian term of damping

		////////////////////TEST GRAD_F correctness///////////
		//f(v+[e,0,0,0...]) - f(v) / e = df/dx
		// VectorXd x_k_test = x_k;
		// x_k_test(0)+=0.000001;
		// TVk_prev = ImplicitXtoTV(x_k_test);
		// SparseMatrix<double> left = ImplicitCalculateElasticForceGradient(TVk_prev);
		// VectorXd right = (ImplicitCalculateForces(TVk_prev) - ImplicitCalculateForces(TVk))/0.000001;
		
		// cout<<"test Grad F"<<endl;
		// cout<<right.transpose()<<endl;
		// cout<<left.col(0).transpose()<<endl<<endl;
		//////////////////////////////////////////////////////

		////////////////////TEST GRAD_G correctness///////////
		//f(v+[e,0,0,0...]) - f(v) / e = df/dx
		//////////////////////////////////////////////////////
		
		//solve for 
		ConjugateGradient<SparseMatrix<double>> cg;
		cg.compute(grad_g);
		VectorXd deltaV = -1*cg.solve(g);
		
		v_k = v_k + deltaV;
		x_k = x_k + timestep*v_k;

		//cout<<g.squaredNorm()<<endl;
		if(g.squaredNorm()<0.000001){
			//cout<<i<<endl;
			break;
		}
	}
	//cout<<endl;
	setXIntoTV(x_k);
	x_old = x_k;
	v_old = v_k;
	
	////////////////////
	// if(dataFile.is_open()){
	// 	double xp=0;
	// 	double yp=0;
	// 	double zp=0;
	// 	for(int i=0; i<v_old.rows(); ){
	// 		xp += v_old(i)*vertex_masses(i);
	// 		i++;
	// 		yp += v_old(i)*vertex_masses(i);
	// 		i++;
	// 		zp += v_old(i)*vertex_masses(i);
	// 		i++;
	// 	}
	// 	dataFile<<t<<","<<xp<<","<<yp<<","<<zp<<"\n";
	// }else{
	// 	cout<<"no open file"<<endl;
	// }
	///////////////////
}


MatrixXi TT_One_G;
MatrixXd TV_One_G;

Simulation Sim;

bool drawLoopTest(igl::viewer::Viewer& viewer){
	viewer.data.clear();
	cout<<Sim.t<<endl;
	Sim.renderImplicit();
	

	viewer.data.add_points(Sim.TV, RowVector3d(1,0,0));
	viewer.data.add_edges(Sim.TV.row(0), Sim.TV.row(1), RowVector3d(1,0,0));
	viewer.data.add_edges(Sim.TV.row(0), Sim.TV.row(2), RowVector3d(1,0,0));
	viewer.data.add_edges(Sim.TV.row(0), Sim.TV.row(3), RowVector3d(1,0,0));
	viewer.data.add_edges(Sim.TV.row(1), Sim.TV.row(2), RowVector3d(1,0,0));
	viewer.data.add_edges(Sim.TV.row(1), Sim.TV.row(3), RowVector3d(1,0,0));
	viewer.data.add_edges(Sim.TV.row(2), Sim.TV.row(3), RowVector3d(1,0,0));

	viewer.data.add_edges(Sim.TV.row(4), Sim.TV.row(3), RowVector3d(0,1,0));
	viewer.data.add_edges(Sim.TV.row(4), Sim.TV.row(0), RowVector3d(0,0,1));
	viewer.data.add_edges(Sim.TV.row(4), Sim.TV.row(2), RowVector3d(0,0,0));

	viewer.data.add_edges(Sim.TV.row(5), Sim.TV.row(3), RowVector3d(0,1,0));
	viewer.data.add_edges(Sim.TV.row(5), Sim.TV.row(0), RowVector3d(0,0,1));
	viewer.data.add_edges(Sim.TV.row(5), Sim.TV.row(1), RowVector3d(0,0,0));

	viewer.data.add_edges(Sim.TV.row(6), Sim.TV.row(3), RowVector3d(0,1,0));
	viewer.data.add_edges(Sim.TV.row(6), Sim.TV.row(0), RowVector3d(0,0,1));
	viewer.data.add_edges(Sim.TV.row(6), Sim.TV.row(5), RowVector3d(0,0,0));


	viewer.data.add_edges(RowVector3d(0,0,0), RowVector3d(100,0,0), RowVector3d(1,1,1));
	viewer.data.add_edges(RowVector3d(0,0,0), RowVector3d(0,100,0), RowVector3d(1,1,1));
	//viewer.data.add_edges(RowVector3d(0,0,0), RowVector3d(0,0,100), RowVector3d(1,1,1));
	
	return false;
}

bool drawLoop(igl::viewer::Viewer& viewer){
	viewer.data.clear();
	cout<<Sim.t<<endl;
	Sim.renderImplicit();
	Sim.renderImplicit();
	Sim.renderImplicit();
	Sim.renderImplicit();
	Sim.renderImplicit();
	for(int i=0; i<Sim.mapV2TV.size(); i++){
		V.row(i) = Sim.TV.row(Sim.mapV2TV[i]);
	}

	viewer.data.add_edges(RowVector3d(0,0,0), RowVector3d(100,0,0), RowVector3d(1,1,1));
	viewer.data.add_edges(RowVector3d(0,0,0), RowVector3d(0,100,0), RowVector3d(1,1,1));
	//viewer.data.add_edges(RowVector3d(0,0,0), RowVector3d(0,0,100), RowVector3d(1,1,1));
	
	viewer.data.set_mesh(V, F);
	viewer.data.set_face_based(false);
	return false;
}

void useFullObject(){
	vector<int> mapV2TV;
	// Load a surface mesh

	// igl::readOFF(TUTORIAL_SHARED_PATH "/3holes.off",V,F);
	igl::readOBJ(TUTORIAL_SHARED_PATH "/wheel.obj", V, F);
	// Tetrahedralize the interior
	igl::copyleft::tetgen::tetrahedralize(V,F,"pq1.414Y", TV,TT,TF);
	// // Compute barycenters
	igl::barycenter(TV, TT,B);
	// //constructs the map from V to TV
	for(int i=0; i< V.rows(); i++){
		for(int j=0; j<TV.rows(); j++){
			if( (V.row(i)-TV.row(j)).squaredNorm() <= 0.001){
				mapV2TV.push_back(j);
				break;
			}
		}
	}
	Sim.initializeSimulation(0.001, TT, TV, mapV2TV);
	Sim.fixVertices(1);
	igl::viewer::Viewer viewer;
	viewer.callback_pre_draw = &drawLoop;
	viewer.launch();
}

void useMyObject(){
	vector<int> mapV2TV;

	TT_One_G.resize(4, 4);
	TT_One_G<< 0, 2, 1, 3,
				4, 2, 0, 3,
				5, 3, 0, 1,
				6, 3, 0, 5;

	TV_One_G.resize(7, 3);
	TV_One_G << 0, 0, 10, //affect
				0, 10, 0,
				10, 0, 0,
				0, 0, 0,
				10, -10, 0,
				-10, 10, 0,
				-10, 0, 0;
				
	Sim.initializeSimulation(0.001, TT_One_G, TV_One_G, mapV2TV);
	Sim.fixVertices(1);
	// int i=0;
	// while(i<2){
	// 	i++;
	// 	cout<<"--------------------------------"<<endl;
	// 	Sim.render();
	// 	cout<<"--------------------------------"<<endl;
	// }

	igl::viewer::Viewer viewer;
	bool boolVariable = true;
	double timeVariable = 0.001;
	viewer.callback_init = [&](igl::viewer::Viewer& viewer)
	{
		viewer.ngui-> addGroup("My settings");

		viewer.ngui-> addVariable("bool", boolVariable);
		viewer.ngui-> addVariable("double", timeVariable);

		// viewer.ngui-> addButton("Print Hello", [](){cout<<"Hello"<<endl;});

		viewer.screen-> performLayout();
		return false;
	};
	viewer.callback_pre_draw = &drawLoopTest;
	viewer.launch();
}
int main(int argc, char *argv[])
{
	dataFile.open("../PythonScripts/data.txt");
	cout<<"###########################My Code ###################"<<endl;
	useMyObject();
	// useFullObject();
	cout<<"###########################My Code ###################"<<endl;
	dataFile.close();
	return 0;
}