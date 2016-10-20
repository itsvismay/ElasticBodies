// input arguements

// logistical arguements
// --mesh + name             -> the name of the mesh to set up
// --config + name           -> the name of the file that contains the config info
// --output + name           -> the name of the output force file

// alignment arguements
// --corner + #              -> the corner of bounding shape to use as origin
// --cubeAlign + height      -> create bounding shape of cubes of various height
// --topCubeAlign + depth    -> create bounding shape of only the top
// --botCubeAlign + depth    -> create bounding shape of only the bot
// --rightCubeAlign + depth  -> create bounding shape of only the right
// --leftCubeAlign + depth   -> create bounding shape of only the left
// --frontCubeAlign + depth  -> create bounding shape of only the front
// --backCubeAlign + depth   -> create bounding shape of only the back

// force arguements
// --maxForce + #            -> maximum amoung of force applied
// --constant + #            -> constant force to apply
// --impulse + #             -> single impulse to apply
// --bezierFile + file.txt   -> name of bezier representation of continuous force

// force area arguements
// --topCubeForce + depth    -> apply the force to the top cube of verts
// --botCubeForce + depth    -> apply the force to the bot cube of verts
// --rightCubeForce + depth  -> apply the force to the right cube of verts
// --leftCubeForce + depth   -> apply the force to the left cube of verts
// --frontCubeForce + depth  -> apply the force to the front cube of verts
// --backCubeVerts + depth   -> apply the force to the back cube of verts

#include <iostream>
#include <fstream>
#include <string>
#include <glm/glm.hpp>

using namespace std;
using namespace glm;

int main(int argc, char* argv[]) {
  cout << "HELLO WORLD" << endl;
  if (argc == 2) return -1;

  // stages of the program
  // read in all the verts and faces
  // calculate the bounding shape based on the verts
  // translate the object to the origin
  // calculate the force bounds and find all forces within it
  // print all verts with forces into the data file
}
