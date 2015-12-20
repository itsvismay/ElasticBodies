#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <pthread.h>
#include <fstream>
#include <math.h>

using namespace Eigen;
using namespace std;

typedef Matrix<double, 12, 1> Vector12d;

class Tetrahedron
{
public:
    Matrix3d DeformedShapeMatrix, ReferenceShapeMatrix, InvRefShapeMatrix;
    VectorXi verticesIndex;
    double undeformedVol;
    Tetrahedron(VectorXi k);
    MatrixXd computeElasticForces(MatrixXd& TV, int e);
    void precomputeForces(MatrixXd& TV);
    VectorXd computeForceDifferentials(MatrixXd& TV, Vector12d& dx);
    Matrix3d computeDeltaDs(const Vector12d& dx);
    Matrix3d computeDs(const Vector12d& x);

};
