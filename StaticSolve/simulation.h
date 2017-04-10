#ifndef simulation_h
#define simulation_h

#include "../ImplicitEuler.h"
#include "../ImplicitNewmark.h"

#include <igl/writeOBJ.h>
#include <igl/barycenter.h>
#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/read_triangle_mesh.h>
#include <Eigen/SPQRSupport>

class Simulation{

public:
	SolidMesh M;
	IntegratorAbstract* integrator;
	vector<int> mapV2TV;
	int iters;
	MatrixXd sB, TV_k;
	VectorXd x_k, f_k, external_force;
	int ignorePastIndex;
	VectorXi putForceOnTheseVerts;
	double maxDisp = 100;

	Simulation(void);
	int initializeSimulation(double deltaT, int iterations, char method, MatrixXi& TT, MatrixXd& TV, MatrixXd& B, vector<int>& moveVertices, vector<int> fixVertices, double youngs, double poissons);

	void staticSolveNewtonsForces(MatrixXd& TV, MatrixXi& TT, MatrixXd& B, VectorXd& fixed_forces, int ignorePastIndex);
	void binarySearchYoungs(vector<int> moveVertices, MatrixXd& TV, MatrixXi& TT, int fv, MatrixXd& B);
	void staticSolveStepNewtonsMethod(double move_step, int ignorePastIndex, vector<int>& moveVertices, MatrixXd& TV,  MatrixXi& TT);
	void syntheticTests(vector<int> moveVertices, MatrixXd& TV, MatrixXi& TT, int fv, MatrixXd& B);
	void reIndexTVandTT(vector<int> newVertsIndices,
						int sizeFixed,
						int sizeMove,
						MatrixXd& TV,
						MatrixXi& TT,
						VectorXd& force,
						MatrixXd& newTV,
						MatrixXi& newTT,
						VectorXd& new_force);

	void staticSolveStepLBFGS(double move_step, int ignorePastIndex, vector<int>& moveVertices, MatrixXd& TV,  MatrixXi& TT);
	void staticSolveInitialPosition(vector<int> moveVertices, MatrixXd& TV, MatrixXi& TT, int fv, MatrixXd& B);

	void setInitPosition(VectorXd& force, vector<int>& fixVertices,  vector<int>& moveVertices);
	void printObj(string printToHere, int numberOfPrints, MatrixXd& TV, MatrixXi& TT, MatrixXd& B);
	void setTVtoX(VectorXd &x, MatrixXd &TV);
	void xToTV(VectorXd& x, MatrixXd& TV);
	void applyExternalForces();
	void calculateElasticForces(VectorXd &f, MatrixXd &TV);
	void calculateForceGradient(MatrixXd &TVk, SparseMatrix<double>& forceGradient);
	void headless();
	void render();

};
#endif
