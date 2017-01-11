#include "boundingVolume.h"
#include <iostream>

using namespace std;

BoundingVolume::BoundingVolume() {
  xBounds = dvec2(10000.0, -10000.0);
  yBounds = dvec2(10000.0, -10000.0);
  zBounds = dvec2(10000.0, -10000.0);
}

void BoundingVolume::addVert(FVert* vert) {
  if (vert->vert[0] < xBounds[0]) xBounds[0] = vert->vert[0];
  if (vert->vert[0] > xBounds[1]) xBounds[1] = vert->vert[0];
  if (vert->vert[1] < yBounds[0]) yBounds[0] = vert->vert[1];
  if (vert->vert[1] > yBounds[1]) yBounds[1] = vert->vert[1];
  if (vert->vert[2] < zBounds[0]) zBounds[0] = vert->vert[2];
  if (vert->vert[2] > zBounds[1]) zBounds[1] = vert->vert[2];
  verts.push_back(vert);
}

void BoundingVolume::distributeForce(double maxForce) {
  if (verts.size() == 0) {
    cout << "ERROR :: No Verts to Apply Forces To" << endl;
    return;
  }
  double forcePerVert = maxForce / verts.size();
  for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
    (*it)->force += forcePerVert;
  }
}

void BoundingVolume::makeFixed() {
  for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
    (*it)->fixed = true;
  }
}

// Note :: z is up, x is left, y is back
// corner values:
// 0 = +x,+y,+z
// 1 = -x,+y,+z
// 2 = -x,-y,+z
// 3 = +x,-y,+z
// 4 = +x,+y,-z
// 5 = -x,+y,-z
// 6 = -x,-y,-z
// 7 = +x,-y,-z
dvec3 BoundingVolume::getNewOrigin(int corner) {
  if (corner == 0) return vec3(xBounds[1], yBounds[1], zBounds[1]);
  if (corner == 1) return vec3(xBounds[0], yBounds[1], zBounds[1]);
  if (corner == 2) return vec3(xBounds[0], yBounds[0], zBounds[1]);
  if (corner == 3) return vec3(xBounds[1], yBounds[0], zBounds[1]);
  if (corner == 4) return vec3(xBounds[1], yBounds[1], zBounds[0]);
  if (corner == 5) return vec3(xBounds[0], yBounds[1], zBounds[0]);
  if (corner == 6) return vec3(xBounds[0], yBounds[0], zBounds[0]);
  if (corner == 7) return vec3(xBounds[1], yBounds[0], zBounds[0]);
  return dvec3(0.0, 0.0, 0.0);
}
