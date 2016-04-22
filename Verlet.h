#ifndef VERLET__H
#define VERLET__H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <pthread.h>
#include <fstream>
#include <math.h>

#include "IntegratorsAbstract.h"

using namespace Eigen;
using namespace std;

class Verlet: public IntegratorAbstract{

public:
	void initializeIntegrator(int ph, SolidMesh& pM, MatrixXd& pTV);
	void setXIntoTV(VectorXd& x_new);
	void createForceVector();
	void calculateGravity();
	void calculateElasticForce();
	void render();
};
#endif