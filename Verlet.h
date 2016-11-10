#ifndef verlet_h
#define verlet_h

#include "IntegratorsAbstract.h"

class Verlet: public IntegratorAbstract{

public:
	void initializeIntegrator(double ph, SolidMesh& pM, MatrixXd& pTV, MatrixXi& pTT);
	void setXIntoTV(VectorXd& x_new);
	void createForceVector();
	void calculateGravity();
	void calculateElasticForce();
	void render(VectorXd& ext_force);
};

#endif