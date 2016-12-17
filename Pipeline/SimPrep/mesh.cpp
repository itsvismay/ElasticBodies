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
  xExtremes[0] = 1000.0;
  xExtremes[1] = -1000.0;
  yExtremes[0] = 1000.0;
  yExtremes[1] = -1000.0;
  zExtremes[0] = 1000.0;
  zExtremes[1] = -1000.0;
  Loader::loadMesh(this, file);
}

void Mesh::writeToFile(ProgramSettings* settings) {
  ofstream objFileWrite(settings->outputMesh);
  ofstream forceFileWrite(settings->outputForce);

  // write geometry data
  //if (settings->needsAlignment) {
  if (true) {
    for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
      objFileWrite << "v " << (*it)->vert[0] << " " << (*it)->vert[1] << " " << (*it)->vert[2] << endl;
    }
    for (vector<ivec3>::iterator it = faces.begin(); it != faces.end(); ++it) {
      objFileWrite << "f " << (*it)[0] << " " << (*it)[1] << " " << (*it)[2] << endl;
    }
  }

  // write force data
  //if (settings->needsForce) {
  //if (true) {
  //  for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
  //    forceFileWrite << (*it)->force * settings->direction[0] << " " << (*it)->force * settings->direction[1] << " " << (*it)->force * settings->direction[2] << endl;
  //  }
  //}

  //if (settings->needsForce) {
  if (true) {
    cout << "THIS SHOULD BE WORKING" << endl;
    for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
      if ((*it)->fixed)
        forceFileWrite << 0 << " " << 0 << " " << 0 << " " << 1 << endl;
      else
        forceFileWrite << (*it)->force * 0 << " " << (*it)->force * -1 << " " << (*it)->force * 0 << " " << 0 << endl;
    }
  }

  objFileWrite.close();
  forceFileWrite.close();
}

void Mesh::translate(double x, double y, double z) {
  for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
    (*it)->translate(x,y,z);
  }
  xExtremes[0] += x;
  xExtremes[1] += x;
  yExtremes[0] += y;
  yExtremes[1] += y;
  zExtremes[0] += z;
  zExtremes[1] += z;
}

void Mesh::addVert(dvec3 vert) {
  if (vert[0] < xExtremes[0]) xExtremes[0] = vert[0];
  if (vert[0] > xExtremes[1]) xExtremes[1] = vert[0];
  if (vert[1] < yExtremes[0]) yExtremes[0] = vert[1];
  if (vert[1] > yExtremes[1]) yExtremes[1] = vert[1];
  if (vert[2] < zExtremes[0]) zExtremes[0] = vert[2];
  if (vert[2] > zExtremes[1]) zExtremes[1] = vert[2];
  verts.push_back(new FVert(vert, verts.size()));
}

void Mesh::addFace(ivec3 face) {
  faces.push_back(face);
}

BoundingVolume* Mesh::createCubeBound(double height) {
  // will be implemented in a future iteration if needed
  cout << "Mesh::createCubeBound not implemented" << endl;
  return 0x0;
}

BoundingVolume* Mesh::createTopBound(double depth) {
  BoundingVolume* volume = new BoundingVolume();
  for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
    dvec3 vert = (*it)->vert;
    if (vert[2] >= zExtremes[1] - depth)
      volume->addVert(*it);
  }
  return volume;
}

BoundingVolume* Mesh::createBotBound(double depth) {
  BoundingVolume* volume = new BoundingVolume();
  for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
    dvec3 vert = (*it)->vert;
    if (vert[2] <= zExtremes[0] + depth)
      volume->addVert(*it);
  }
  return volume;
}

BoundingVolume* Mesh::createRightBound(double depth) {
  BoundingVolume* volume = new BoundingVolume();
  for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
    dvec3 vert = (*it)->vert;
    if (vert[0] >= xExtremes[1] - depth)
      volume->addVert(*it);
  }
  return volume;
}

BoundingVolume* Mesh::createLeftBound(double depth) {
  cout << "ZExtreemes :: " << zExtremes[0] << endl;
  BoundingVolume* volume = new BoundingVolume();
  cout << "DEPTH: " << depth << endl;
  for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
    dvec3 vert = (*it)->vert;
    if (vert[0] <= xExtremes[0] + depth)
      volume->addVert(*it);
  }
  return volume;
}

BoundingVolume* Mesh::createFrontBound(double depth) {
  BoundingVolume* volume = new BoundingVolume();
  for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
    dvec3 vert = (*it)->vert;
    if (vert[1] <= yExtremes[0] + depth)
      volume->addVert(*it);
  }
  return volume;
}

BoundingVolume* Mesh::createBackBound(double depth) {
  BoundingVolume* volume = new BoundingVolume();
  for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
    dvec3 vert = (*it)->vert;
    if (vert[1] >= yExtremes[1] - depth)
      volume->addVert(*it);
  }
  return volume;
}

BoundingVolume* Mesh::createLineBound(double x, double y, double z, bool alongX, bool alongY, bool alongZ) {
  BoundingVolume* volume = new BoundingVolume();
  cout << "BOUNDING X: " << x << endl;
  cout << "BOUNDING Y: " << y << endl;
  cout << "BOUNDING Z: " << z << endl;
  if (alongX) {
    for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
      dvec3 vert = (*it)->vert;
      if (vert[1] == y && vert[2] == z)
        volume->addVert(*it);
    }
  } else if (alongY) {
    for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
      dvec3 vert = (*it)->vert;
      if (vert[0] == x && vert[2] == z)
        volume->addVert(*it);
    }
  } else if (alongZ) {
    for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
      dvec3 vert = (*it)->vert;
      if (vert[0] == x && vert[1] == y)
        volume->addVert(*it);
    }
  }
  //cout << "SIZE: " volume->verts.size();
  return volume;
}

dvec2 Mesh::xBnds() { return xExtremes; }
dvec2 Mesh::yBnds() { return yExtremes; }
dvec2 Mesh::zBnds() { return zExtremes; }
