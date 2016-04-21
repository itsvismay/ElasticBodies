#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <pthread.h>
#include <fstream>
#include <math.h>


#include "solidmesh.h"

using namespace Eigen;
using namespace std;

class Integrator{

protected:
	int t=0;
	double timestep;
	SparseMatrix<double> InvMass;
	vector<int> fixedVertices;
	int vertices; //number of vertices

	SolidMesh M; //NEW IDEA: Pass in from simulation
	MatrixXd TV, TVk;

	VectorXd x_old, v_old, x_k, v_k, vertex_masses, f;

public:
	virtual void render() = 0; //pure virtual render class
	void initializeIntegrator()	
}