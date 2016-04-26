#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <pthread.h>
#include <fstream>
#include <math.h>


#include "Verlet.h"
#include "ImplicitEuler.h"
#include "ImplicitNewmark.h"


using namespace Eigen;
using namespace std;

class Simulation{

public:
	SolidMesh M;
	IntegratorAbstract* integrator;
	vector<int> mapV2TV;

	Simulation(void);
	void render();
	void initializeSimulation(double deltaT, char method, MatrixXi& TT, MatrixXd& TV, vector<int> map);
};