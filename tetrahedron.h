#ifndef tetrahedron_h
#define tetrahedron_h

#include "globals.h"

class Tetrahedron
{
public:
    Matrix3d DeformedShapeMatrix, ReferenceShapeMatrix, InvRefShapeMatrix;
    VectorXi verticesIndex;
    double undeformedVol, energyDensity;
    double mu, lambda;
    Tetrahedron(VectorXi k, double mu, double lambda);
    MatrixXd computeElasticForces(MatrixXd& TV, int e);
    void precomputeForces(MatrixXd& TV);
    MatrixXd computeForceDifferentials(MatrixXd& TV, Vector12d& dx);
    Matrix3d computeDeltaDs(const Vector12d& dx);
    Matrix3d computeDs(const Vector12d& x);

};

#endif