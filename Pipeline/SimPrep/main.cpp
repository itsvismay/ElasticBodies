// input arguements

// logistical arguements - 3
// --mesh + name             -> the name of the mesh to set up
// --config + name           -> the name of the file that contains the config info
// --output + name           -> the name of the output force file

// alignment arguements - 8
// --corner + #              -> the corner of bounding shape to use as origin
// --cubeAlign + height      -> create bounding shape of cubes of various height
// --topCubeAlign + depth    -> create bounding shape of only the top
// --botCubeAlign + depth    -> create bounding shape of only the bot
// --rightCubeAlign + depth  -> create bounding shape of only the right
// --leftCubeAlign + depth   -> create bounding shape of only the left
// --frontCubeAlign + depth  -> create bounding shape of only the front
// --backCubeAlign + depth   -> create bounding shape of only the back

// force arguements - 4
// --maxForce + #            -> maximum amoung of force applied
// --constant + #            -> constant force to apply
// --impulse + #             -> single impulse to apply
// --bezierFile + file.txt   -> name of bezier representation of continuous force

// force area arguements - 6
// --topCubeForce + depth    -> apply the force to the top cube of verts
// --botCubeForce + depth    -> apply the force to the bot cube of verts
// --rightCubeForce + depth  -> apply the force to the right cube of verts
// --leftCubeForce + depth   -> apply the force to the left cube of verts
// --frontCubeForce + depth  -> apply the force to the front cube of verts
// --backCubeForce + depth   -> apply the force to the back cube of verts

// 10
// force domain arguements - will be included later also will be slower
// to be used to define a custom domain, then detect verts in that domain
// two of the following three must be defined, all three is an error
// --domainMinX + #          -> define min x domain bounds
// --domainMaxX + #          -> define max x domain bounds
// --domainMinY + #          -> define min y domain bounds
// --domainMaxY + #          -> define max y domain bounds
// --domainMinZ + #          -> define min z domain bounds
// --domainMaxZ + #          -> define max z domain bounds
// --domainBezier + file.txt -> define bezier surface domain bounds
// --domainDirX + #          -> define domain direction to be pos or neg x
// --domainDirY + #          -> define domain direction to be pos or neg y
// --domainDirZ + #          -> define domain direction to be pos or neg z

#include <iostream>
#include <fstream>
#include <string.h>
#include <glm/glm.hpp>
#include "mesh.h"

using namespace std;
using namespace glm;

// function declarations
void empty();

// TODO :: CHANGE ALL OF THIS

int main(int argc, char* argv[]) {
  cout << "HELLO WORLD" << endl;
  //if (argc == 2) return -1;

  int currentArg = 1;

  ProgramSettings* settings = new ProgramSettings();

  // take in arguements and settings
  while (currentArg != argc) {
    char* arg = argv[currentArg];

    if (strcmp(arg, "--in") == 0) {
      settings->mesh = string(argv[++currentArg]);
    } else if (strcmp(arg, "--out") == 0) {
      cout << "FOUND OUTPUT:: " << endl;
      settings->outputMesh = string(argv[++currentArg]);
    } else if (strcmp(arg, "--config") == 0) {
      settings->config = string(argv[++currentArg]);
      cout << "FOUND FORCE FILE" << endl;
    } else if (strcmp(arg, "--force") == 0) {
      settings->outputForce = string(argv[++currentArg]);
    } else if (strcmp(arg, "--corner") == 0) {
      settings->corner = atoi(argv[++currentArg]);
    } else if (strcmp(arg, "--cubeAlign") == 0) {
      settings->needsAlignment = true;
      settings->cubeAlign = true;
      settings->depth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--topCubeAlign") == 0) {
      settings->needsAlignment = true;
      settings->topCubeAlign = true;
      settings->depth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--botCubeAlign") == 0) {
      settings->needsAlignment = true;
      settings->botCubeAlign = true;
      settings->depth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--rightCubeAlign") == 0) {
      settings->needsAlignment = true;
      settings->rightCubeAlign = true;
      settings->depth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--leftCubeAlign") == 0) {
      settings->needsAlignment = true;
      settings->leftCubeAlign = true;
      settings->depth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--frontCubeAlign") == 0) {
      settings->needsAlignment = true;
      settings->frontCubeAlign = true;
      settings->depth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--backCubeAlign") == 0) {
      settings->needsAlignment = true;
      settings->backCubeAlign = true;
      settings->depth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--maxForce") == 0) {
      settings->needsForce = true;
      cout << "FOUND MAX FORCE" << endl;
      settings->maxForce = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--constant") == 0) {
      settings->needsForce = true;
      settings->constantForce = true;
      settings->maxForce = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--impulse") == 0) {
      settings->needsForce = true;
      settings->impulseForce = true;
      settings->maxForce = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--bezierFile") == 0) {
      settings->bezierForce = true;
      settings->bezierForceFile = string(argv[++currentArg]);
    } else if (strcmp(arg, "--topCubeForce") == 0) {
      settings->topCubeForce = true;
      settings->forceDepth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--botCubeForce") == 0) {
      settings->botCubeForce = true;
      settings->forceDepth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--rightCubeForce") == 0) {
      settings->rightCubeForce = true;
      settings->forceDepth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--leftCubeForce") == 0) {
      settings->leftCubeForce = true;
      settings->forceDepth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--frontCubeForce") == 0) {
      settings->frontCubeForce = true;
      settings->forceDepth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--backCubeForce") == 0) {
      settings->backCubeForce = true;
      settings->forceDepth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--domainMinX") == 0) {
      settings->domainForce = true;
      settings->domainX[0] = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--domainMaxX") == 0) {
      settings->domainForce = true;
      settings->domainX[1] = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--domainMinY") == 0) {
      settings->domainForce = true;
      settings->domainY[0] = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--domainMaxY") == 0) {
      settings->domainForce = true;
      settings->domainY[1] = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--domainMinZ") == 0) {
      settings->domainForce = true;
      settings->domainZ[0] = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--domainMaxZ") == 0) {
      settings->domainForce = true;
      settings->domainZ[1] = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--domainBezier") == 0) {
      settings->domainBezier = true;
    } else if (strcmp(arg, "--domainDirX") == 0) {
      settings->domainForce = true;
      settings->direction[0] = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--domainDirY") == 0) {
      settings->domainForce = true;
      settings->direction[1] = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--domainDirZ") == 0) {
      settings->domainForce = true;
      settings->direction[2] = atof(argv[++currentArg]);
    }

    currentArg++;
  }

  // stages of the program
  // read in all the verts and faces
  cout << "Making Mesh" << endl;
  Mesh* mesh = new Mesh(settings->mesh);
  cout << "Made Mesh" << endl;
  // calculate the bounding shape based on the verts
  BoundingVolume* volume = 0x0;
  //if (settings->needsAlignment) {
  if (true) {
    // initialize bounding volume
    //if (settings->cubeAlign) volume = mesh->createCubeBound(settings->depth);
    //else if (settings->topCubeAlign) volume = mesh->createTopBound(settings->depth);
    //else if (settings->botCubeAlign) volume = mesh->createBotBound(settings->depth);
    //else if (settings->rightCubeAlign) volume = mesh->createRightBound(settings->depth);
    //else if (settings->leftCubeAlign) volume = mesh->createLeftBound(settings->depth);
    //else if (settings->frontCubeAlign) volume = mesh->createFrontBound(settings->depth);
    //else if (settings->backCubeAlign) volume = mesh->createBackBound(settings->depth);
    cout << "Creating Fixed Volume" << endl;
    volume = mesh->createLeftBound(0.0);
    cout << "Created Fixed Volume" << endl;
    // translate the object to the origin
    //vec3 newOrigin = volume->getNewOrigin(settings->corner);
    cout << "Getting Origin" << endl;
    dvec3 newOrigin = volume->getNewOrigin(1);
    cout << "ORIGIN:: " << newOrigin[0] << " " << newOrigin[1] << " " << newOrigin[2] << endl;
    //mesh->translate(-newOrigin[0], -newOrigin[1], -newOrigin[2]);
    mesh->translate(-mesh->xExtremes[0], -mesh->yExtremes[1], -mesh->zExtremes[1]);
    cout << "Making Fixed" << endl;
    volume->makeFixed();
  }
  // calculate the force bounds and calculate all forces within it
  BoundingVolume* forceVolume = 0x0;
  //if (settings->needsForce
  cout << "Calculating Force Volume" << endl;
  if (true)
  {
    // initialize force bounding volume
    //if (settings->topCubeForce) forceVolume = mesh->createTopBound(settings->forceDepth);
    //else if (settings->botCubeForce) forceVolume = mesh->createBotBound(settings->forceDepth);
    //else if (settings->rightCubeForce) forceVolume = mesh->createRightBound(settings->forceDepth);
    //else if (settings->leftCubeForce) forceVolume = mesh->createLeftBound(settings->forceDepth);
    //else if (settings->frontCubeForce) forceVolume = mesh->createFrontBound(settings->forceDepth);
    //else if (settings->backCubeForce) forceVolume = mesh->createBackBound(settings->forceDepth);
    cout << "FIXED BOUNDS :: " << volume->verts.size() << endl;
    cout << "Creating Line Bound" << endl;
    cout << "MESH XBNDS :: " << mesh->xBnds()[0] << " " << mesh->xBnds()[1] << endl;
    cout << "MESH YBNDS :: " << mesh->yBnds()[0] << " " << mesh->yBnds()[1] << endl;
    cout << "MESH ZBNDS :: " << mesh->zBnds()[0] << " " << mesh->zBnds()[1] << endl;
    forceVolume = mesh->createLineBound(mesh->xBnds()[1], 0.0, mesh->zBnds()[1], false, true, false);
    // distribute force per vert in volume
    cout << "Distributing Max Force" << endl;
    cout << "MAX FORCE :: " << settings->maxForce << endl;
    forceVolume->distributeForce(settings->maxForce);
  }
  // print all verts with forces into the data file
  cout << "Writing All Files" << endl;
  mesh->writeToFile(settings);
  cout << "Cleaning Up" << endl;
  // clean up
  if (forceVolume) delete forceVolume;
  if (volume) delete volume;
  delete mesh;
  delete settings;
}

void empty() {
  // does nothing and is a place holder for non implemented things
}
