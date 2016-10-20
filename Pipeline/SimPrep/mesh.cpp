#include "mesh.h"

#include <iostream>

Mesh::Mesh(string) {
  file = f;
  readFromFile();
}

Mesh::~Mesh() {
  for (int i=0;i<verts.size();i++) {
    delete verts[i];
  }
}

void Mesh::readFromFile() {
  xExtremes[0] = 1000.0f;
  xExtremes[1] = -1000.0f;
  yExtremes[0] = 1000.0f;
  yExtremes[1] = -1000.0f;
  zExtremes[0] = 1000.0f;
  zExtremes[1] = -1000.0f;
  Loader.loadMesh(this, file);
}

void Mesh::writeToFile() {
  ofstream fileWrite(file);

  for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
    fileWrite << "v " << (*it)->vert[0] << " " << (*it)->vert[1] << " " (*it)->vert[2] << endl;
  }
  for (vector<ivec4>::iterator it = faces.begin(); it != faces.end(); ++it) {
    fileWrite << "f " << (*it)[0] << " " << (*it)[1] << " " (*it)[2] << " " << (*it)[3] << endl;
  }

  fileWrite.close();
}

void Mesh::translate(float x, float y, float z) {
  for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
    (*it)->translate(x,y,z);
  }
}

BoundingVolume* Mesh::createCubeBound(float height) {
  // will be implemented in a future iteration if needed
  cout << "Mesh::createCubeBound not implemented" << endl;
  return 0x0;
}

BoundingVolume* Mesh::createTopBound(float depth) {
  BoundingVolume* volume = new BoundingVolume();
  for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
    vec3 vert = (*it)->vert;
    if (vert[2] >= zExtremes[1] - depth)
      volume->addVert(*it);
  }
  return volume;
}

BoundingVolume* Mesh::createBotBound(float depth) {
  BoundingVolume* volume = new BoundingVolume();
  for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
    vec3 vert = (*it)->vert;
    if (vert[2] <= zExtremes[0] + depth)
      volume->addVert(*it);
  }
  return volume;
}

BoundingVolume* Mesh::createRightBound(float depth) {
  BoundingVolume* volume = new BoundingVolume();
  for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
    vec3 vert = (*it)->vert;
    if (vert[0] >= zExtremes[1] - depth)
      volume->addVert(*it);
  }
  return volume;
}

BoundingVolume* Mesh::createLeftBound(float depth) {
  BoundingVolume* volume = new BoundingVolume();
  for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
    vec3 vert = (*it)->vert;
    if (vert[0] <= zExtremes[0] + depth)
      volume->addVert(*it);
  }
  return volume;
}

BoundingVolume* Mesh::createFrontBound(float depth) {
  BoundingVolume* volume = new BoundingVolume();
  for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
    vec3 vert = (*it)->vert;
    if (vert[1] <= zExtremes[0] + depth)
      volume->addVert(*it);
  }
  return volume;
}

BoundingVolume* Mesh::createBackBound(float depth) {
  BoundingVolume* volume = new BoundingVolume();
  for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
    vec3 vert = (*it)->vert;
    if (vert[1] >= zExtremes[0] - depth)
      volume->addVert(*it);
  }
  return volume;
}
