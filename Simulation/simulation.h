#ifndef simulation_h
#define simulation_h

#include "../ImplicitEuler.h"
#include "../ImplicitNewmark.h"

#include <igl/writeOBJ.h>
#include <igl/writeMESH.h>
#include <igl/barycenter.h>
#include <igl/readOFF.h>
#include <igl/writeOFF.h>
#include <igl/hausdorff.h>
#include <igl/readOBJ.h>
#include <igl/readMESH.h>
#include <Eigen/SPQRSupport>
#include <Eigen/SparseQR>

class Simulation{

public:
	SolidMesh M;
	IntegratorAbstract* integrator;
	vector<int> moveVerticesStore;
	int iters;
	MatrixXd TV_k;
	MatrixXd* sB;
	VectorXd x_k, f_k, external_force;
	int ignorePastIndex;
	VectorXi putForceOnTheseVerts;
	double maxDisp = 100;

	MatrixXd V;
	MatrixXi F;

	Simulation(void);
	int initializeSimulation(double deltaT, int iterations, char method, MatrixXi& TT, MatrixXd& TV, MatrixXd& B, vector<int>& moveVertices, vector<int> fixVertices, double youngs, double poissons);

	void staticSolveNewtonsForces(MatrixXd& TV, MatrixXi& TT, MatrixXd& B, VectorXd& fixed_forces, vector<int>& moveVertices, int ignorePastIndex, int step);

	void staticSolveNewtonsPosition(MatrixXd& TV, MatrixXi& TT, MatrixXd& B, vector<int>& moveVertices, int ignorePastIndex, int step);
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

	void applyStaticForces(MatrixXd& TV, MatrixXi& TT, MatrixXd& B, VectorXd& fixed_forces, vector<int>& moveVertices, vector<int>& fixVertices);
	void applyStaticPositions(MatrixXd& TV, MatrixXi& TT, MatrixXd& B, VectorXd& fixed_forces, vector<int>& moveVertices, vector<int>& fixVertices);

	void setInitPosition(VectorXd& force, vector<int>& fixVertices, vector<int>& moveVertices);

	void printObj(string printToHere, int numberOfPrints, MatrixXd& TV, MatrixXi& TT, MatrixXd& B);
	void printDesigns(int printcount, int simTime);
	double printOptimizationOutput();

	void setTVtoX(VectorXd &x, MatrixXd &TV);
	void xToTV(VectorXd& x, MatrixXd& TV);
	void applyExternalForces();
	void calculateElasticForces(VectorXd &f, MatrixXd &TV);
	void calculateForceGradient(MatrixXd &TVk, SparseMatrix<double>& forceGradient);
	void headless();
	bool render();

};
#endif
