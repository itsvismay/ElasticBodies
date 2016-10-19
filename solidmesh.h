#ifndef solidmesh_h
#define solidmesh_h

#include "tetrahedron.h"

class SolidMesh
{
    public:
        vector<Tetrahedron> tets;
        SolidMesh();
        SolidMesh(MatrixXi& TT, MatrixXd& TV, double youngs, double poissons);
        void initializeMesh(MatrixXi& TT, MatrixXd& TV, double youngs, double poissons);

        void setNewYoungsPoissons(double youngs, double poissons);
};

#endif