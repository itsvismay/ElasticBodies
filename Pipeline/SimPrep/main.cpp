// input arguements

// logistical arguements - 3
// --mesh + name             -> the name of the mesh to set up
// --meshObj + name          -> the name of the obj mesh to set up
// --meshOff + name          -> the name of the off mesh to set up
// --config + name           -> the name of the file that contains the config info
// --outputF + name          -> the name of the output force file no extension
// --outputM + name          -> the name of the output mesh file no extension

// alignment arguements - 1
// --corner + #              -> the corner of bounding shape to use as origin

// force arguements - 4
// --maxForce + #            -> maximum amoung of force applied
// --constant + #            -> constant force to apply
// --impulse + #             -> single impulse to apply
// --bezierFile + file.txt   -> name of bezier representation of continuous force

// force area arguements - 6
// --topForce + depth    -> apply the force to the top box of verts
// --botForce + depth    -> apply the force to the bot box of verts
// --rightForce + depth  -> apply the force to the right box of verts
// --leftForce + depth   -> apply the force to the left box of verts
// --frontForce + depth  -> apply the force to the front box of verts
// --backForce + depth   -> apply the force to the back box

// fixed area arguements - 6
// --topFixed + depth    -> make all verts in the top box fixed
// --botFixed + depth    -> make all verts in the bot box fixed
// --rightFixed + depth  -> make all verts in the right box fixed
// --leftFixed + depth   -> make all verts in the left box fixed
// --frontFixed + depth  -> make all verts in the front box fixed
// --backFixed + depth   -> make all verts in the back box fixed

// force direction arguements - 6
// --forceXComp + value  -> sets the x component of the force direction
// --forceYComp + value  -> sets the y component of the force direction
// --forceZComp + value  -> sets the z component of the force direction
// --forceXAxis + scale  -> sets the force to be along the x axis
// --forceYAxis + scale  -> sets the force to be along the y axis
// --forceZAxis + scale  -> sets the force to be along the z axis

#include <iostream>
#include <fstream>
#include <string.h>
#include <glm/glm.hpp>
#include "mesh.h"

using namespace std;
using namespace glm;

int main(int argc, char* argv[]) {
  int currentArg = 1;
  ProgramSettings* settings = new ProgramSettings();

  // take in arguements and settings
  while (currentArg != argc) {
    char* arg = argv[currentArg];
    ////////// LEGACY SUPPORT ////////////
    if (strcmp(arg, "--in") == 0) {
      settings->inputMeshObj = string(argv[++currentArg]);
      settings->obj = true;
    } else if (strcmp(arg, "--out") == 0) {
      settings->outputMeshObj = string(argv[++currentArg]);
      settings->obj = true;
    //////////////////////////////////////
    } else if (strcmp(arg, "--mesh") == 0) {
      settings->inputMeshObj = string(argv[++currentArg]);
      settings->obj = true;
    } else if (strcmp(arg, "--meshObj") == 0) {
      settings->inputMeshObj = string(argv[++currentArg]);
      settings->obj = true;
    } else if (strcmp(arg, "--meshOff") == 0) {
      settings->inputMeshOff = string(argv[++currentArg]);
      settings->off = true;
    } else if (strcmp(arg, "--config") == 0) {
      settings->config = string(argv[++currentArg]);
    } else if (strcmp(arg, "--outputF") == 0) {
      settings->outputForce = string(argv[++currentArg]);
    } else if (strcmp(arg, "--outputM") == 0) {
      if (settings->off) settings->outputMeshOff = string(argv[++currentArg]);
      else settings->outputMeshObj = string(argv[++currentArg]);
    } else if (strcmp(arg, "--force") == 0) {
      settings->outputForce = string(argv[++currentArg]);
      // to be removed in legacy
    } else if (strcmp(arg, "--corner") == 0) {
      settings->corner = atoi(argv[++currentArg]);
    } else if (strcmp(arg, "--maxForce") == 0) {
      settings->maxForce = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--constant") == 0) {
      settings->constantForce = true;
      settings->maxForce = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--impulse") == 0) {
      settings->impulseForce = true;
      settings->maxForce = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--bezierFile") == 0) {
      settings->bezierForce = true;
      settings->bezierForceFile = string(argv[++currentArg]);
    } else if (strcmp(arg, "--topForce") == 0) {
      settings->topForce = true;
      settings->forceDepth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--botForce") == 0) {
      settings->botForce = true;
      settings->forceDepth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--rightForce") == 0) {
      settings->rightForce = true;
      settings->forceDepth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--leftForce") == 0) {
      settings->leftForce = true;
      settings->forceDepth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--frontForce") == 0) {
      settings->frontForce = true;
      settings->forceDepth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--backForce") == 0) {
      settings->backForce = true;
      settings->forceDepth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--topFixed") == 0) {
      settings->topFixed = true;
      settings->fixedDepth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--botFixed") == 0) {
      settings->botFixed = true;
      settings->fixedDepth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--rightFixed") == 0) {
      settings->rightFixed = true;
      settings->fixedDepth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--leftFixed") == 0) {
      settings->leftFixed = true;
      settings->fixedDepth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--frontFixed") == 0) {
      settings->frontFixed = true;
      settings->fixedDepth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--backFixed") == 0) {
      settings->backFixed = true;
      settings->fixedDepth = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--forceXComp") == 0) {
      settings->forceDirection[0] = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--forceYComp") == 0) {
      settings->forceDirection[1] = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--forceZComp") == 0) {
      settings->forceDirection[2] = atof(argv[++currentArg]);
    } else if (strcmp(arg, "--forceXAxis") == 0) {
      settings->forceDirection[0] = atof(argv[++currentArg]);
      settings->forceDirection[1] = 0.0;
      settings->forceDirection[2] = 0.0;
    } else if (strcmp(arg, "--forceYAxis") == 0) {
      settings->forceDirection[0] = 0.0;
      settings->forceDirection[1] = atof(argv[++currentArg]);
      settings->forceDirection[2] = 0.0;
    } else if (strcmp(arg, "--forceZAxis") == 0) {
      settings->forceDirection[0] = 0.0;
      settings->forceDirection[1] = 0.0;
      settings->forceDirection[2] = atof(argv[++currentArg]);
    } else {
      cout << "ERROR :: INPUT " << arg << " NOT RECOGNIZED" << endl;
      exit(-1);
    }

    currentArg++;
  }

  // stages of the program
  // read in all the verts and faces
  Mesh* mesh;
  if (settings->off) mesh = new Mesh(settings->inputMeshOff, settings);
  else mesh = new Mesh(settings->inputMeshObj, settings);

  // calculate the bounding shape based on the verts
  BoundingVolume* fixedVolume = 0x0;
  if (settings->topFixed) fixedVolume = mesh->createTopBound(settings->fixedDepth);
  else if (settings->botFixed) fixedVolume = mesh->createBotBound(settings->fixedDepth);
  else if (settings->rightFixed) fixedVolume = mesh->createRightBound(settings->fixedDepth);
  else if (settings->leftFixed) fixedVolume = mesh->createLeftBound(settings->fixedDepth);
  else if (settings->frontFixed) fixedVolume = mesh->createFrontBound(settings->fixedDepth);
  else if (settings->backFixed) fixedVolume = mesh->createBackBound(settings->fixedDepth);

  /////// temp code ///////
  fixedVolume = mesh->createFrontBound(.01);
  /////////////////////////

  // translate the object to the origin
  //if (settings->corner == -1) settings->corner = 1; // corner not implemented yet
  mesh->translate(-mesh->xBounds[0], -mesh->yBounds[1], -mesh->zBounds[1]);

  // fix verts
  fixedVolume->makeFixed();

  // calculate the force bounds and calculate all forces within it
  BoundingVolume* forceVolume = 0x0;
  if (settings->topForce) forceVolume = mesh->createTopBound(settings->forceDepth);
  else if (settings->botForce) forceVolume = mesh->createBotBound(settings->forceDepth);
  else if (settings->rightForce) forceVolume = mesh->createRightBound(settings->forceDepth);
  else if (settings->leftForce) forceVolume = mesh->createLeftBound(settings->forceDepth);
  else if (settings->frontForce) forceVolume = mesh->createFrontBound(settings->forceDepth);
  else if (settings->backForce) forceVolume = mesh->createBackBound(settings->forceDepth);

  /////// temp code ///////
  //forceVolume = mesh->createLineBound(mesh->xBnds()[1], 0.0, mesh->zBnds()[1], false, true, false);
  forceVolume = mesh->createBackBound(0.01);
  /////////////////////////

  // distribute force per vert in volume
  forceVolume->distributeForce(settings->maxForce);

  // print all verts with forces into the data file
  if (settings->obj) mesh->writeObjToFile(settings);
  if (settings->off) mesh->writeOffToFile(settings);

  // clean up
  if (forceVolume) delete forceVolume;
  if (fixedVolume) delete fixedVolume;
  delete mesh;
  delete settings;
}
