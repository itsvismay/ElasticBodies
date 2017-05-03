#ifndef globals_h
#define globals_h

#include "externals.h"

//MACROS
#define HOME_SAVED_PATH "/home/vismay/ElasticBodies/"
#define OUTPUT_SAVED_PATH "/home/vismay/ElasticBodies/"

// #define HOME_SAVED_PATH "/u/vismay/ElasticBodies/"
// #define OUTPUT_SAVED_PATH "/scratch/cluster/vismay/"


// #define TUTORIAL_SHARED_PATH "/u/vismay/ElasticBodies/"
// #define CONSISTENCY_TEST_SAVE_PATH "/scratch/cluster/vismay/"

// #define TUTORIAL_SHARED_PATH "/home/vismay/ElasticBodies/"
// #define CONSISTENCY_TEST_SAVE_PATH "/home/vismay/ElasticBodies/"

using namespace Eigen;
using namespace std;

typedef Eigen::Triplet<double> Trip;
typedef Matrix<double, 12, 1> Vector12d;

extern ofstream momentumFile;
extern ofstream energyFile;
extern ofstream strainEnergyFile;
extern ofstream kineticEnergyFile;
extern ofstream gravityEnergyFile;
extern ofstream optimizationFile;

extern double rayleighCoeff;
extern double gravity;
extern bool headless;
extern string material_model;
extern string tetgen_code;
extern string solver;
extern string objectName;
extern double youngs;

extern MatrixXi TF;
#endif
