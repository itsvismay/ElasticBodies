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
    this->undeformedVol = (1.0/6)*fabs(Dm.determinant());
}

Matrix3d Tetrahedron::computeDeltaDs(const Vector12d& dx){
	// double xio, yio, zio, xjo, yjo, zjo, xko, yko, zko, xlo, ylo, zlo;
 //    //Vector3d ro1 = TV.row(this->verticesIndex(0));
	// xio = dx(0);
 //    yio = dx(1);
 //    zio = dx(2);

 //    //Vector3d ro2 = TV.row(this->verticesIndex(1));
 //    xjo = dx(3);
 //    yjo = dx(4);
 //    zjo = dx(5);

 //    //Vector3d ro3 = TV.row(this->verticesIndex(2));
 //    xko = dx(6);
 //    yko = dx(7);
 //    zko = dx(8);

 //    //Vector3d ro4 = TV.row(this->verticesIndex(3));
 //    xlo = dx(9);
 //    ylo = dx(10);
 //    zlo = dx(11);

    Matrix3d dDs;
    dDs << (dx(0) - dx(9)), (dx(3) - dx(9)), (dx(6) - dx(9)),
            (dx(1) - dx(10)), (dx(4) - dx(10)), (dx(7) - dx(10)),
            (dx(2) - dx(11)), (dx(5) - dx(11)), (dx(8) - dx(11));
     return dDs;
}


MatrixXd Tetrahedron::computeForceDifferentials(MatrixXd& TV, Vector12d& dx){

    Matrix3d Ds;
    for(int i=0; i<3; i++){
        Ds.col(i) = TV.row(verticesIndex(i)) - TV.row(verticesIndex(3));
    }

    Matrix3d dDs;
    dDs <<  (dx(0) - dx(9)),   (dx(3) - dx(9)), (dx(6) - dx(9)),
            (dx(1) - dx(10)), (dx(4) - dx(10)), (dx(7) - dx(10)),
            (dx(2) - dx(11)), (dx(5) - dx(11)), (dx(8) - dx(11));

    Matrix3d F = Ds*this->InvRefShapeMatrix;
    Matrix3d dF = dDs*this->InvRefShapeMatrix;

    Matrix3d dP;
    if(material_model.compare("neo") == 0 ){
        //Neohookean
        double detF = F.determinant();
        double logdetF = log(detF);
        Matrix3d FInvTransp = (F.inverse()).transpose();
        dP = mu*dF + (mu - lambda*logdetF)*(FInvTransp)*dF.transpose()*(FInvTransp) + lambda*(F.inverse()*dF).trace()*(FInvTransp);
    }
    else if(material_model.compare("svk") == 0){
        //SVK
        Matrix3d E = 0.5*((F.transpose()*F) - MatrixXd::Identity(3,3));
        Matrix3d P = F*(2*mu*E + lambda*E.trace()*MatrixXd::Identity(3,3));
        Matrix3d dE = 0.5*(F.transpose()*dF + dF.transpose()*F);
        dP = 2*mu*dF*E + 2*mu*F*dE + lambda*E.trace()*dF*MatrixXd::Identity(3,3) + lambda*dE.trace()*F*MatrixXd::Identity(3,3);

    }
    else{
        cout<<"Material model not specified properly"<<endl;
        exit(0);
    }

    Matrix3d dH = -1*this->undeformedVol*dP*((this->InvRefShapeMatrix).transpose());


    MatrixXd dForces(3,4);
    dForces.col(0) = dH.col(0);
    dForces.col(1) = dH.col(1);
    dForces.col(2) = dH.col(2);
    dForces.col(3) = -1*dH.col(0) - dH.col(1) - dH.col(2);

    return dForces;
}

void Tetrahedron::computeElasticForces(MatrixXd &TV, VectorXd& f){

    Matrix3d Ds;
    for(int i=0; i<3; i++){
        Ds.col(i) = TV.row(verticesIndex(i)) - TV.row(verticesIndex(3));
    }

    Matrix3d F = Ds*this->InvRefShapeMatrix;

    Matrix3d P;

    this->currentVol = (1.0/6)*fabs(Ds.determinant());
    if(material_model.compare("neo") == 0){
        //Neo
        double J = F.determinant();
        double I1 = (F.transpose()*F).trace();
        double powj = pow(J, -2.0/3.0);
        double I1bar = powj*I1;

        P = mu*(powj * F) +
        (- mu/3.0 * I1 * powj + lambda*(J-1.0)*J)*F.inverse().transpose();

        this->energy = this->undeformedVol*(mu/2.0 * (I1bar - 3) + lambda/2.0 * (J-1.0) * (J-1.0));

    }else if(material_model.compare("svk") == 0){
        //SVK
        //TODO: Spring Constant value
        Matrix3d E = 0.5*((F.transpose()*F) - MatrixXd::Identity(3,3));
        P = F*(2*mu*E + lambda*E.trace()*MatrixXd::Identity(3,3));//piola kirchoff
        this->energy = (mu*(E*E).trace() + (lambda/2)*E.trace()*E.trace())*this->undeformedVol;

    }else{
        cout<<"Material model not specified properly"<<endl;
        exit(0);
    }
    if(F.determinant()<0){
        this->energy = 1e40;
        cout<<"ERROR: F determinant is 0"<<endl;
        cout<<"Decrease timestep maybe - instantaneous force is too much with this timestep"<<endl;
        exit(0);
    }
    if(this->energy != this->energy){
        //NANS
        cout<<"ENERGY nans"<<endl;
        exit(0);
    }

    Matrix3d H = -1*this->undeformedVol*P*((this->InvRefShapeMatrix).transpose());

    f.segment<3>(3*verticesIndex(0)) += H.col(0);
    f.segment<3>(3*verticesIndex(1)) += H.col(1);
    f.segment<3>(3*verticesIndex(2)) += H.col(2);
    f.segment<3>(3*verticesIndex(3)) += -1*H.col(0) - H.col(1) - H.col(2);

}

MatrixXd Tetrahedron::oldComputeElasticForces(MatrixXd &TV, int e){

    Matrix3d Ds;
    for(int i=0; i<3; i++){
        Ds.col(i) = TV.row(verticesIndex(i)) - TV.row(verticesIndex(3));
    }

    Matrix3d F = Ds*this->InvRefShapeMatrix;

    Matrix3d P;

    this->currentVol = (1.0/6)*fabs(Ds.determinant());
    if(material_model.compare("neo") == 0){
        //Neo

        P = mu*(F - ((F.inverse()).transpose())) + lambda*log(F.determinant())*((F.inverse()).transpose());
        double firstTerm = ((mu/2.0)*((F.transpose()*F).trace() -3) - mu*log(F.determinant()));
        if(firstTerm<0 ){
            firstTerm = 0;
        }
        this->energy = firstTerm + (lambda/2)*log(F.determinant())*log(F.determinant());
        this->energy *= this->undeformedVol;

    }else if(material_model.compare("svk") == 0){
        //SVK
        //TODO: Spring Constant value
        Matrix3d E = 0.5*((F.transpose()*F) - MatrixXd::Identity(3,3));
        P = F*(2*mu*E + lambda*E.trace()*MatrixXd::Identity(3,3));//piola kirchoff
        this->energy = (mu*(E*E).trace() + (lambda/2)*E.trace()*E.trace())*this->undeformedVol;

    }else{
        cout<<"Material model not specified properly"<<endl;
        exit(0);
    }
    if(F.determinant()<0){
        this->energy = 1e40;
        cout<<"ERROR: F determinant is 0"<<endl;
    }
    if(this->energy != this->energy){
        //NANS
        cout<<"ENERGY nans"<<endl;
        exit(0);
    }

    Matrix3d H = -1*this->undeformedVol*P*((this->InvRefShapeMatrix).transpose());
    Matrix<double, 3, 4> Forces;
    Forces.col(0) = H.col(0);
    Forces.col(1) = H.col(1);
    Forces.col(2) = H.col(2);
    Forces.col(3) = -1*H.col(0) - H.col(1) - H.col(2);

    return Forces;
}
