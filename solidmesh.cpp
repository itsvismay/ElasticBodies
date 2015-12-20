#include <Eigen/Core>
#include <Eigen/Sparse>
#include "solidmesh.h"

using namespace Eigen;
using namespace std;

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