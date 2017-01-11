#include "mesh.h"
#include "loader.h"
#include <iostream>

Mesh::Mesh(string f, ProgramSettings* settings) {
  file = f.c_str();
  if (settings->off) readOffFromFile();
  else readObjFromFile();
}

Mesh::~Mesh() {
  for (int i=0;i<verts.size();i++) {
    delete verts[i];
  }
}

void Mesh::readObjFromFile() {
  xBounds[0] = 1000.0;
  xBounds[1] = -1000.0;
  yBounds[0] = 1000.0;
  yBounds[1] = -1000.0;
  zBounds[0] = 1000.0;
  zBounds[1] = -1000.0;
  Loader::loadMeshOBJ(this, file);
}

void Mesh::readOffFromFile() {
  xBounds[0] = 1000.0;
  xBounds[1] = -1000.0;
  yBounds[0] = 1000.0;
  yBounds[1] = -1000.0;
  zBounds[0] = 1000.0;
  zBounds[1] = -1000.0;
  Loader::loadMeshOFF(this, file);
}

void Mesh::writeObjToFile(ProgramSettings* settings) {
  ofstream objFileWrite(settings->outputMeshObj);
  ofstream forceFileWrite(settings->outputForce);

  for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
    objFileWrite << "v " << (*it)->vert[0] << " " << (*it)->vert[1] << " " << (*it)->vert[2] << endl;
  }
  for (vector<ivec3>::iterator it = faces.begin(); it != faces.end(); ++it) {
    objFileWrite << "f " << (*it)[0] << " " << (*it)[1] << " " << (*it)[2] << endl;
  }

  for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
    dvec3 dir = settings->forceDirection;
    if ((*it)->fixed) {
      forceFileWrite << (*it)->force * dir[0] << " " << (*it)->force * dir[1] << " " << (*it)->force * dir[2] << " " << 1 << endl;
    }
    else
      forceFileWrite << (*it)->force * dir[0] << " " << (*it)->force * dir[1] << " " << (*it)->force * dir[2] << " " << 0 << endl;
  }

  objFileWrite.close();
  forceFileWrite.close();
}

void Mesh::writeOffToFile(ProgramSettings* settings) {
  ofstream offFileWrite(settings->outputMeshOff);
  ofstream forceFileWrite(settings->outputForce);

  offFileWrite << "OFF " << verts.size() << " " << faces.size() << " 0" << endl;
  for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
    offFileWrite << (*it)->vert[0] << " " << (*it)->vert[1] << " " << (*it)->vert[2] << endl;
  }
  for (vector<ivec3>::iterator it = faces.begin(); it != faces.end(); ++it) {
    offFileWrite << "3 " << (*it)[0] << " " << (*it)[1] << " " << (*it)[2] << endl;
  }

  for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
    dvec3 dir = settings->forceDirection;
    if ((*it)->fixed) {
      forceFileWrite << (*it)->force * dir[0] << " " << (*it)->force * dir[1] << " " << (*it)->force * dir[2] << " " << 1 << endl;
    }
    else
      forceFileWrite << (*it)->force * dir[0] << " " << (*it)->force * dir[1] << " " << (*it)->force * dir[2] << " " << 0 << endl;
  }

  offFileWrite.close();
  forceFileWrite.close();
}

void Mesh::translate(double x, double y, double z) {
  for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
    (*it)->translate(x,y,z);
  }
  xBounds[0] += x;
  xBounds[1] += x;
  yBounds[0] += y;
  yBounds[1] += y;
  zBounds[0] += z;
  zBounds[1] += z;
}

void Mesh::addVert(dvec3 vert) {
  if (vert[0] < xBounds[0]) xBounds[0] = vert[0];
  if (vert[0] > xBounds[1]) xBounds[1] = vert[0];
  if (vert[1] < yBounds[0]) yBounds[0] = vert[1];
  if (vert[1] > yBounds[1]) yBounds[1] = vert[1];
  if (vert[2] < zBounds[0]) zBounds[0] = vert[2];
  if (vert[2] > zBounds[1]) zBounds[1] = vert[2];
  verts.push_back(new FVert(vert, verts.size()));
}

void Mesh::addFace(ivec3 face) {
  faces.push_back(face);
}

BoundingVolume* Mesh::createTopBound(double depth) {
  BoundingVolume* volume = new BoundingVolume();
  for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
    dvec3 vert = (*it)->vert;
    if (vert[2] >= zBounds[1] - depth)
      volume->addVert(*it);
  }
  return volume;
}

BoundingVolume* Mesh::createBotBound(double depth) {
  BoundingVolume* volume = new BoundingVolume();
  for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
    dvec3 vert = (*it)->vert;
    if (vert[2] <= zBounds[0] + depth)
      volume->addVert(*it);
  }
  return volume;
}

BoundingVolume* Mesh::createRightBound(double depth) {
  BoundingVolume* volume = new BoundingVolume();
  for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
    dvec3 vert = (*it)->vert;
    if (vert[0] >= xBounds[1] - depth)
      volume->addVert(*it);
  }
  return volume;
}

BoundingVolume* Mesh::createLeftBound(double depth) {
  BoundingVolume* volume = new BoundingVolume();
  for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
    dvec3 vert = (*it)->vert;
    if (vert[0] <= xBounds[0] + depth)
      volume->addVert(*it);
  }
  return volume;
}

BoundingVolume* Mesh::createFrontBound(double depth) {
  BoundingVolume* volume = new BoundingVolume();
  for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
    dvec3 vert = (*it)->vert;
    if (vert[1] <= yBounds[0] + depth)
      volume->addVert(*it);
  }
  return volume;
}

BoundingVolume* Mesh::createBackBound(double depth) {
  BoundingVolume* volume = new BoundingVolume();
  for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
    dvec3 vert = (*it)->vert;
    if (vert[1] >= yBounds[1] - depth)
      volume->addVert(*it);
  }
  return volume;
}

BoundingVolume* Mesh::createLineBound(double x, double y, double z, bool alongX, bool alongY, bool alongZ) {
  BoundingVolume* volume = new BoundingVolume();
  if (alongX) {
    for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
      dvec3 vert = (*it)->vert;
      if (vert[1] >= y -.01 && vert[1] <= y + .01 && vert[2] >= z - .01 && vert[2] <= z + .01)
        volume->addVert(*it);
    }
  } else if (alongY) {
    for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
      dvec3 vert = (*it)->vert;
      if (vert[0] >= x -.01 && vert[0] <= x + .01 && vert[2] >= z - .01 && vert[2] <= z + .01)
        volume->addVert(*it);
    }
  } else if (alongZ) {
    for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
      dvec3 vert = (*it)->vert;
      if (vert[0] == x && vert[1] == y)
        volume->addVert(*it);
    }
  }
  return volume;
}

dvec2 Mesh::xBnds() { return xBounds; }
dvec2 Mesh::yBnds() { return yBounds; }
dvec2 Mesh::zBnds() { return zBounds; }
