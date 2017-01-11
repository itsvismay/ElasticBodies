#include "../globals.h"
#include <igl/readMESH.h>
#include <igl/writeMESH.h>
#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>

using namespace Eigen;
int main(int argc, char *argv[])
{

	Eigen::MatrixXd TV;
	Eigen::MatrixXi TT;
	Eigen::MatrixXi TF;
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	string filepath = argv[1];
	std::cout<< filepath+argv[2]<<std::endl;
	igl::readOFF(filepath+argv[2], V, F);
	igl::copyleft::tetgen::tetrahedralize(V,F, "pR", TV,TT,TF);

	igl::writeMESH(filepath+argv[3]+".mesh", TV, TT, TF);

}
