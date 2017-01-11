#include <igl/read_triangle_mesh.h>
#include <igl/copyleft/cgal/CSGTree.h>
#include <igl/viewer/Viewer.h>
#include <igl/jet.h>
#include <Eigen/Core>
#include <vector>
#include <igl/writeOBJ.h>

#include "tutorial_shared_path.h"

#define PATH_TO_OBJ_FILES "/home/firal/Documents/Research/ElasticBodies/Pipeline/"

using namespace Eigen;
using namespace igl::copyleft::cgal;
using namespace std;
using namespace igl;

// Hella memory inefficient, 
// make this better in the future
// use pointer or some shit to keep track of splits
CSGTree* binaryMergeTreeCreation(vector<MatrixXd>& Vs, vector<MatrixXi>& Fs){
  if(Vs.size()==2 && Fs.size()==2){
    CSGTree* M;
    M = new CSGTree(CSGTree(Vs[0], Fs[0]),CSGTree(Vs[1], Fs[1]), "u");
    return M;
  }
  else if(Vs.size() ==1 && Fs.size()==1){
    CSGTree* M;
    M = new CSGTree(Vs[0], Fs[0]);
    return M;
  }
  else{
    vector<MatrixXd> Vs2(make_move_iterator(Vs.begin() + Vs.size()/2),make_move_iterator(Vs.end()));
    Vs.erase(Vs.begin() + Vs.size()/2, Vs.end());
    vector<MatrixXi> Fs2(make_move_iterator(Fs.begin() + Fs.size()/2),make_move_iterator(Fs.end()));
    Fs.erase(Fs.begin() + Fs.size()/2, Fs.end());
    CSGTree* M;
    CSGTree* M1 = binaryMergeTreeCreation(Vs, Fs);
    CSGTree* M2 = binaryMergeTreeCreation(Vs2, Fs2);
    M = new CSGTree(CSGTree(M1->V(), M1->F()),CSGTree(M2->V(), M2->F()), "u");
    delete M1;
    delete M2;
    return M;
  }
}



int main(int argc, char * argv[])
{
  vector<MatrixXi> Fs;
  vector<MatrixXd> Vs;
  // Read arguements
  int num_of_layers = atoi(argv[2]);
  string name(argv[1]);
  // Read in inputs as double precision floating point meshes
  if (argc == 3) {
    for(int i=0; i<num_of_layers; i++){
      MatrixXi F;
      MatrixXd V;
      read_triangle_mesh(PATH_TO_OBJ_FILES +name+"_layer_" + to_string(i)+".obj",V,F);
      Fs.push_back(F);
      Vs.push_back(V);  
    }
  } else if (argc == 4){
    string path(argv[3]);
    for(int i=0; i<num_of_layers; i++) {
      MatrixXi F;
      MatrixXd V;
      read_triangle_mesh(path + name + "_layer_" + to_string(i) + ".obj", V, F);
      Fs.push_back(F);
      Vs.push_back(V);
    }
  }

  CSGTree* M;

  M = binaryMergeTreeCreation( Vs, Fs);
  
  writeOBJ(PATH_TO_OBJ_FILES"unioned.obj", M->cast_V<MatrixXd>(), M->F());
}

