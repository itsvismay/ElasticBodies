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

    if (strcmp(arg, "--mesh")) {
      settings->mesh = string(argv[++currentArg]);
    } else if (strcmp(arg, "--config")) {
      settings->config = string(argv[++currentArg]);
    } else if (strcmp(arg, "--output")) {
      settings->outputForce = string(argv[++currentArg]);
    } else if (strcmp(arg, "--corner")) {
      settings->corner = atoi(argv[++currentArg]);
    } else if (strcmp(arg, "--cubeAlign")) {
      settings->needsAlignment = true;
      settings->cubeAlign = true;
      settings->depth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--topCubeAlign")) {
      settings->needsAlignment = true;
      settings->topCubeAlign = true;
      settings->depth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--botCubeAlign")) {
      settings->needsAlignment = true;
      settings->botCubeAlign = true;
      settings->depth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--rightCubeAlign")) {
      settings->needsAlignment = true;
      settings->rightCubeAlign = true;
      settings->depth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--leftCubeAlign")) {
      settings->needsAlignment = true;
      settings->leftCubeAlign = true;
      settings->depth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--frontCubeAlign")) {
      settings->needsAlignment = true;
      settings->frontCubeAlign = true;
      settings->depth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--backCubeAlign")) {
      settings->needsAlignment = true;
      settings->backCubeAlign = true;
      settings->depth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--maxForce")) {
      settings->needsForce = true;
      settings->maxForce = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--constant")) {
      settings->needsForce = true;
      settings->constantForce = true;
      settings->maxForce = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--impulse")) {
      settings->needsForce = true;
      settings->impulseForce = true;
      settings->maxForce = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--bezierFile")) {
      settings->bezierForce = true;
      settings->bezierForceFile = string(argv[++currentArg]);
    } else if (strcmp(arg, "--topCubeForce")) {
      settings->topCubeForce = true;
      settings->forceDepth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--botCubeForce")) {
      settings->botCubeForce = true;
      settings->forceDepth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--rightCubeForce")) {
      settings->rightCubeForce = true;
      settings->forceDepth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--leftCubeForce")) {
      settings->leftCubeForce = true;
      settings->forceDepth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--frontCubeForce")) {
      settings->frontCubeForce = true;
      settings->forceDepth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--backCubeForce")) {
      settings->backCubeForce = true;
      settings->forceDepth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--domainMinX")) {
      settings->domainForce = true;
      settings->domainX[0] = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--domainMaxX")) {
      settings->domainForce = true;
      settings->domainX[1] = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--domainMinY")) {
      settings->domainForce = true; 
      settings->domainY[0] = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--domainMaxY")) {
      settings->domainForce = true;
      settings->domainY[1] = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--domainMinZ")) {
      settings->domainForce = true;
      settings->domainZ[0] = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--domainMaxZ")) {
      settings->domainForce = true;
      settings->domainZ[1] = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--domainBezier")) {
      settings->domainBezier = true;
    } else if (strcmp(arg, "--domainDirX")) {
      settings->domainForce = true;
      settings->direction[0] = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--domainDirY")) {
      settings->domainForce = true;
      settings->direction[1] = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--domainDirZ")) {
      settings->domainForce = true;
      settings->direction[2] = atof(argv[++currentArg]);
    }

    currentArg++;
  }

  // stages of the program
  // read in all the verts and faces
  Mesh* mesh = new Mesh(settings->mesh);
  // calculate the bounding shape based on the verts
  BoundingVolume* volume = 0x0;
  if (settings->needsAlignment) {
    // initialize bounding volume
    if (settings->cubeAlign) volume = mesh->createCubeBound(settings->depth);
    else if (settings->topCubeAlign) volume = mesh->createTopBound(settings->depth);
    else if (settings->botCubeAlign) volume = mesh->createBotBound(settings->depth);
    else if (settings->rightCubeAlign) volume = mesh->createRightBound(settings->depth);
    else if (settings->leftCubeAlign) volume = mesh->createLeftBound(settings->depth);
    else if (settings->frontCubeAlign) volume = mesh->createFrontBound(settings->depth);
    else if (settings->backCubeAlign) volume = mesh->createBackBound(settings->depth);
    // translate the object to the origin
    vec3 newOrigin = volume->getNewOrigin(settings->corner);
    mesh->translate(-newOrigin[0], -newOrigin[1], -newOrigin[2]);
  }
  // calculate the force bounds and calculate all forces within it
  BoundingVolume* forceVolume = 0x0;
  if (settings->needsForce)
  {
    // initialize force bounding volume
    if (settings->topCubeForce) forceVolume = mesh->createTopBound(settings->forceDepth);
    else if (settings->botCubeForce) forceVolume = mesh->createBotBound(settings->forceDepth);
    else if (settings->rightCubeForce) forceVolume = mesh->createRightBound(settings->forceDepth);
    else if (settings->leftCubeForce) forceVolume = mesh->createLeftBound(settings->forceDepth);
    else if (settings->frontCubeForce) forceVolume = mesh->createFrontBound(settings->forceDepth);
    else if (settings->backCubeForce) forceVolume = mesh->createBackBound(settings->forceDepth);
    // distribute force per vert in volume
    forceVolume->distributeForce(settings->maxForce);
  }
  // print all verts with forces into the data file
  mesh->writeToFile(settings);
  // clean up
  if (forceVolume) delete forceVolume;
  if (volume) delete volume;
  delete mesh;
  delete settings;
}

void empty() {
  // does nothing and is a place holder for non implemented things
}
