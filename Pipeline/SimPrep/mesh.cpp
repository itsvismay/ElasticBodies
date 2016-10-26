#include "mesh.h"
#include "loader.h"
#include <iostream>

Mesh::Mesh(string f) {
  file = f.c_str();
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
  Loader::loadMesh(this, file);
}

void Mesh::writeToFile(ProgramSettings* settings) {
  ofstream objFileWrite(settings->mesh);
  ofstream forceFileWrite(settings->outputForce);

  // write geometry data
  if (settings->needsAlignment) {
    for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
      objFileWrite << "v " << (*it)->vert[0] << " " << (*it)->vert[1] << " " << (*it)->vert[2] << endl;
    }
    for (vector<ivec4>::iterator it = faces.begin(); it != faces.end(); ++it) {
      objFileWrite << "f " << (*it)[0] << " " << (*it)[1] << " " << (*it)[2] << " " << (*it)[3] << endl;
    }
  }

  // write force data
  if (settings->needsForce) {
    for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
      forceFileWrite << (*it)->force * settings->direction[0] << " " << (*it)->force * settings->direction[1] << " " << (*it)->force * settings->direction[2] << endl;
    }
  }

  objFileWrite.close();
  forceFileWrite.close();
}

void Mesh::translate(float x, float y, float z) {
  for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
    (*it)->translate(x,y,z);
  }
}

void Mesh::addVert(vec3 vert) {
  if (vert[0] < xExtremes[0]) xExtremes[0] = vert[0];
  if (vert[0] > xExtremes[1]) xExtremes[1] = vert[0];
  if (vert[1] < yExtremes[0]) yExtremes[0] = vert[1];
  if (vert[1] > yExtremes[1]) yExtremes[1] = vert[1];
  if (vert[2] < zExtremes[0]) zExtremes[0] = vert[2];
  if (vert[2] > zExtremes[1]) zExtremes[1] = vert[2];
  verts.push_back(new FVert(vert, verts.size()));
}

void Mesh::addFace(ivec4 face) {
  faces.push_back(face);
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
