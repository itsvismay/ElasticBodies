#include "boundingVolume.h"
#include <iostream>

using namespace std;

BoundingVolume::BoundingVolume() {
  xExtremes = vec2(10000.0, -10000.0);
  yExtremes = vec2(10000.0, -10000.0);
  zExtremes = vec2(10000.0, -10000.0);
}

BoundingVolume::~BoundingVolume() {
  for (int i=0;i<containingVolumes.size();i++) {
    delete containingVolumes[i];
  }
}

void BoundingVolume::addVert(FVert* vert) {
  if (vert->vert[0] < xExtremes[0]) xExtremes[0] = vert->vert[0];
  if (vert->vert[0] > xExtremes[1]) xExtremes[1] = vert->vert[0];
  if (vert->vert[1] < yExtremes[0]) yExtremes[0] = vert->vert[1];
  if (vert->vert[1] > yExtremes[1]) yExtremes[1] = vert->vert[1];
  if (vert->vert[2] < zExtremes[0]) zExtremes[0] = vert->vert[2];
  if (vert->vert[2] > zExtremes[1]) zExtremes[1] = vert->vert[2];
  verts.push_back(vert);
}

void BoundingVolume::addBoundingVolume(BoundingVolume* volume) {
  containingVolumes.push_back(volume);
}

void BoundingVolume::distributeForce(double maxForce) {
  if (verts.size() == 0) {
    cout << "VERTS SIZE IS ZERO FOR FORCE" << endl;
    return;
  }
  cout << "VERTS SIZE :: " << verts.size() << endl;
  float forcePerVert = maxForce / verts.size();
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
vec3 BoundingVolume::getNewOrigin(int corner) {
  if (corner == 0) return vec3(xExtremes[1], yExtremes[1], zExtremes[1]);
  if (corner == 1) return vec3(xExtremes[0], yExtremes[1], zExtremes[1]);
  if (corner == 2) return vec3(xExtremes[0], yExtremes[0], zExtremes[1]);
  if (corner == 3) return vec3(xExtremes[1], yExtremes[0], zExtremes[1]);
  if (corner == 4) return vec3(xExtremes[1], yExtremes[1], zExtremes[0]);
  if (corner == 5) return vec3(xExtremes[0], yExtremes[1], zExtremes[0]);
  if (corner == 6) return vec3(xExtremes[0], yExtremes[0], zExtremes[0]);
  if (corner == 7) return vec3(xExtremes[1], yExtremes[0], zExtremes[0]);
  return vec3(0.0, 0.0, 0.0);
}
