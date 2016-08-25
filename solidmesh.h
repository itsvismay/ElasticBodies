#include <Eigen/Core>
#include <Eigen/Sparse>
#include "tetrahedron.h"

using namespace Eigen;
using namespace std;

class SolidMesh
{
    public:
        vector<Tetrahedron> tets;
        SolidMesh();
        SolidMesh(MatrixXi& TT, MatrixXd& TV, double youngs, double poissons);
        void initializeMesh(MatrixXi& TT, MatrixXd& TV, double youngs, double poissons);
};