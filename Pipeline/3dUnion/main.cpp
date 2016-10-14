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
CSGTree binaryMergeTreeCreation(vector<MatrixXd>& Vs, vector<MatrixXi>& Fs){
  if(Vs.size()==2 && Fs.size()==2){
    CSGTree M;
    M = {{Vs[0], Fs[0]},{Vs[1], Fs[1]}, "u"};
    return M;
  }
  else if(Vs.size() ==1 && Fs.size()==1){
    CSGTree M;
    M = {Vs[0], Fs[0]};
    return M;
  }
  else{

    vector<MatrixXd> Vs2(make_move_iterator(Vs.begin() + Vs.size()/2),make_move_iterator(Vs.end()));
    Vs.erase(Vs.begin() + Vs.size()/2, Vs.end());

    vector<MatrixXi> Fs2(make_move_iterator(Fs.begin() + Fs.size()/2),make_move_iterator(Fs.end()));
    Fs.erase(Fs.begin() + Fs.size()/2, Fs.end());

    CSGTree M;
    CSGTree M1 = binaryMergeTreeCreation(Vs, Fs);
    CSGTree M2 = binaryMergeTreeCreation(Vs2, Fs2);
    M = {{M1.V(), M1.F()},{M2.V(), M2.F()}, "u"};
    return M;
  }
}



int main(int argc, char * argv[])
{
  vector<MatrixXi> Fs;
  vector<MatrixXd> Vs;
  // TODO :: Read in number of layers as input
  // TODO :: Read in input name as input
  // TODO :: Read in output name as input
  int num_of_layers = 11;
  // Read in inputs as double precision floating point meshes
  cout << "WHY IS THIS NOT BEING DISPLAYED" << endl;
  for(int i=0; i<num_of_layers; i++){
    MatrixXi F;
    MatrixXd V;
    //cout << PATH_TO_OBJ_FILES << endl;
    //read_triangle_mesh(PATH_TO_OBJ_FILES "layer"+to_string(i)+".stl.obj",V,F);
    read_triangle_mesh(PATH_TO_OBJ_FILES "test_layer_" + to_string(i)+".obj",V,F);
    Fs.push_back(F);
    Vs.push_back(V);  
  }

  if(Fs.size()!=Vs.size() || Fs.size()<2){
    cout<<"Nope"<<endl;
    exit(0);
  }

  cout<<Fs.size()<<endl;
  cout<<Vs.size()<<endl;
  
  CSGTree M;

  M = binaryMergeTreeCreation( Vs, Fs);
  
  writeOBJ(PATH_TO_OBJ_FILES"unioned.obj", M.cast_V<MatrixXd>(), M.F());

  cout<<"Done"<<endl;
}

