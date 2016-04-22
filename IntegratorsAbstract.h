#ifndef INTEGRATORS__H
#define INTEGRATORS__H

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

class IntegratorAbstract{

public:
	int simTime =0;
	double h; //timestep
	SparseMatrix<double> InvMass;
	SparseMatrix<double> RegMass;
	SparseMatrix<double> ZeroMatrix;
	SparseMatrix<double> Ident;


	vector<int> fixedVerts;
	int vertsNum; //number of vertices

	SolidMesh M; //NEW IDEA: Pass in from simulation
	MatrixXd TV;

	VectorXd x_old, v_old, f, massVector;
	int width;
	int height;

	bool isFixed(int vert);
	virtual void render()=0; //pure virtual render class
	void initializeIntegrator(int ph, SolidMesh& pM, MatrixXd& pTV);
	void initVectors();
	void initMassMatrices();
	void fixVertices(int fixMe);
	void createXFromTet();
};
#endif