#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/LU>
#include <iostream>
#include <vector>
#include <pthread.h>
#include <fstream>
#include <math.h>
#include "tetrahedron.h"

using namespace Eigen;
using namespace std;

Tetrahedron::Tetrahedron(VectorXi k, double mu, double lambda){
    verticesIndex = k ;
    this->mu = mu;
    this->lambda = lambda;
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
    // cout<<this->undeformedVol<<endl;
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

MatrixXd Tetrahedron::computeForceDifferentials(MatrixXd& TV, Vector12d& dx){
	Vector12d x;
    x.segment<3>(0) = TV.row(this->verticesIndex(0));
    x.segment<3>(3) = TV.row(this->verticesIndex(1));
    x.segment<3>(6) = TV.row(this->verticesIndex(2));
    x.segment<3>(9) = TV.row(this->verticesIndex(3));


    Matrix3d Ds = computeDs(x);
    Matrix3d dDs = computeDeltaDs(dx);



    ////////////////////TEST dDs correctness///////////
    // double epsilon = 0.000001;
	// f(v+[e,0,0,0...]) - f(v) / e = df/dx
	// cout<<"test dDs"<<endl;
	// Matrix3d leftdDs = computeDeltaDs(dx*epsilon);
	// Matrix3d rightdDs = computeDs(x + dx*epsilon) - Ds;
	// cout<<leftdDs-rightdDs<<endl<<endl;
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
    double detF = F.determinant();
    double logdetF = log(detF);
    Matrix3d FInvTransp = (F.inverse()).transpose();
    Matrix3d dP = mu*dF + (mu - lambda*logdetF)*(FInvTransp)*dF.transpose()*(FInvTransp) + lambda*(F.inverse()*dF).trace()*(FInvTransp);
    // // Matrix3d P = mu*(F - (FInvTransp)) + lambda*logdetF*(FInvTransp);
    // // OLD Matrix3d P = mu*(F - ((F.inverse()).transpose())) + lambda*log(F.determinant())*((F.inverse()).transpose());
    // // OLD Matrix3d dP = mu*dF + (mu - lambda*log(F.determinant()))*((F.inverse()).transpose())*dF.transpose()*((F.inverse()).transpose()) + lambda*(F.inverse()*dF).trace()*((F.inverse()).transpose());
    
    //SVK
    // Matrix3d E = 0.5*((F.transpose()*F) - MatrixXd::Identity(3,3));
    // Matrix3d P = F*(2*mu*E + lambda*E.trace()*MatrixXd::Identity(3,3));
    // Matrix3d dE = 0.5*(F.transpose()*dF + dF.transpose()*F);
    // Matrix3d dP = 2*mu*dF*E + 2*mu*F*dE + lambda*E.trace()*dF*MatrixXd::Identity(3,3) + lambda*dE.trace()*F*MatrixXd::Identity(3,3);
    
    // Matrix3d dP = dF*(2*mu*E + lambda*E.trace()*MatrixXd::Identity(3,3));
    // dP+= F*(2*mu*(0.5*((dF.transpose()*F + F.transpose()*dF) - MatrixXd::Identity(3,3))));
    // dP+= F*(lambda*(0.5*((dF.transpose()*F + F.transpose()*dF)-MatrixXd::Identity(3,3))).trace());

    ////////////////////TEST dP correctness///////////
 //    cout<<"test dP"<<endl;
 //    Matrix3d leftdP = mu*leftdF + (mu - lambda*log(rightF.determinant()))*((rightF.inverse()).transpose())*leftdF.transpose()*((rightF.inverse()).transpose()) + lambda*(rightF.inverse()*leftdF).trace()*((rightF.inverse()).transpose());
	// Matrix3d rightP1 = mu*(rightF - ((rightF.inverse()).transpose())) + lambda*log(rightF.determinant())*((rightF.inverse()).transpose());
	// Matrix3d rightP2 = mu*(F - ((F.inverse()).transpose())) + lambda*log(F.determinant())*((F.inverse()).transpose());
	// Matrix3d rightdP = rightP1 - rightP2;
 //    cout<<leftdP<<endl;
	// cout<< leftdP - (rightP1 - rightP2)<<endl<<endl;
	//////////////////////////////////////////////////////
    // cout<<"dP matrix"<<endl;
    // cout<<dP<<endl;

    Matrix3d dH = -1*this->undeformedVol*dP*((this->InvRefShapeMatrix).transpose());
    
    ////////////////////TEST dH correctness///////////
 //    cout<<"test dH"<<endl;
 //    Matrix3d leftdH = -1*this->undeformedVol*leftdP*((this->InvRefShapeMatrix).transpose());
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
    return dForces;
}

MatrixXd Tetrahedron::computeElasticForces(MatrixXd &TV, int e){
    Vector12d x;
    x.segment<3>(0) = TV.row(this->verticesIndex(0));
    x.segment<3>(3) = TV.row(this->verticesIndex(1));
    x.segment<3>(6) = TV.row(this->verticesIndex(2));
    x.segment<3>(9) = TV.row(this->verticesIndex(3));

    Matrix3d Ds = computeDs(x);

    Matrix3d F = Ds*this->InvRefShapeMatrix;


    //SVK
    //TODO: Spring Constant value
    // Matrix3d E = 0.5*((F.transpose()*F) - MatrixXd::Identity(3,3));
    // Matrix3d P = F*(2*mu*E + lambda*E.trace()*MatrixXd::Identity(3,3));//piola kirchoff	
	// this->energyDensity = mu*(E*E).trace() + (lambda/2)*E.trace()*E.trace();

    //Neo
	Matrix3d P = mu*(F - ((F.inverse()).transpose())) + lambda*log(F.determinant())*((F.inverse()).transpose());
    this->energyDensity = (mu/2.0)*((F.transpose()*F).trace() -3) - mu*log(F.determinant()) + (lambda/2)*log(F.determinant())*log(F.determinant());


    Matrix3d H = -1*this->undeformedVol*P*((this->InvRefShapeMatrix).transpose());

    Matrix<double, 3, 4> Forces;
    Forces.col(0) = H.col(0);
    Forces.col(1) = H.col(1);
    Forces.col(2) = H.col(2);
    Forces.col(3) = -1*H.col(0) - H.col(1) - H.col(2);

    return Forces;
}