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

class Integrators{

public:
	int t=0;
	double timestep;
	SparseMatrix<double> InvMass;
	vector<int> fixedVertices;
	int vertices;

	SolidMesh M;

	VectorXd x_old, v_old, x_k, v_k, vertex_masses, f;
}