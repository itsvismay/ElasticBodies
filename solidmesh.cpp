#include <Eigen/Core>
#include <Eigen/Sparse>
#include "solidmesh.h"

using namespace Eigen;
using namespace std;

SolidMesh::SolidMesh(void){}

SolidMesh::SolidMesh(MatrixXi& TT, MatrixXd& TV, double youngs, double poissons){
    // double youngs = 2200000; //UNITS: pascals (but divide by 1000 to get from m to mm)
    // double poissons = 0.35;
    double mu = youngs/(2+ 2*poissons);
    double lambda = youngs*poissons/((1+poissons)*(1-2*poissons));
    for(int i=0; i<TT.rows(); i++){
    	int i1 = TT.row(i)[0];
    	int i2 = TT.row(i)[1];
    	int i3 = TT.row(i)[2];
    	int i4 = TT.row(i)[3];

    	Tetrahedron t(TT.row(i), mu, lambda);
    	t.precomputeForces(TV);
    	this->tets.push_back(t);
    }
}

void SolidMesh::initializeMesh(MatrixXi& TT, MatrixXd& TV, double youngs, double poissons){
    // double youngs = 2200000; //UNITS: pascals (but divide by 1000 to get from m to mm)
    // double poissons = 0.35;
    double mu = youngs/(2+ 2*poissons);
    double lambda = youngs*poissons/((1+poissons)*(1-2*poissons));
    cout<<"Setting Mu and Lambda from Youngs and Poissons"<<endl;
    cout<<"Solidmesh init Youngs, poissons ="<<youngs<<", "<<poissons<<endl;
    cout<<"Solidmesh init Mu, Lambda ="<<mu<<", "<<lambda<<endl<<endl;
	for(int i=0; i<TT.rows(); i++){
    	//based on Tet indexes, get the vertices of the tet from TV
    	Tetrahedron t(TT.row(i), mu, lambda);
    	t.precomputeForces(TV);
    	this->tets.push_back(t);
    }
}

void SolidMesh::setNewYoungsPoissons(double youngs, double poissons){
    double mu = youngs/(2+ 2*poissons);
    double lambda = youngs*poissons/((1+poissons)*(1-2*poissons));
    cout<<"NEW** Mu and Lambda from Youngs and Poissons"<<endl;
    cout<<"Solidmesh Youngs, poissons ="<<youngs<<", "<<poissons<<endl;
    cout<<"Solidmesh Mu, Lambda ="<<mu<<", "<<lambda<<endl<<endl;
    for(int i=0; i<tets.size(); i++){
        tets[i].mu = mu;
        tets[i].lambda = lambda;
    }
}