#include "boundingVolume.h"

BoundingVolume::BoundingVolume() {
  // to be implemented
}

BoundingVolume::~BoundingVolume() {
  for (int i=0;i<containingVolumes.size();i++) {
    delete containingVolumes[i];
  }
  for (int i=0;i<verts.size();i++) {
    delete verts[i];
  }
}

void BoundingVolume::addVert(FVert* vert) {
  verts.push_back(vert);
}

void BoundingVolume::addBoundingVolume(BoundingVolume* volume) {
  containingVolumes.push_back(volume);
}

void BoundingVolume::distributeForce(float maxForce) {
  float forcePerVert = maxForce / verts.size();
  for (vector<FVert*>::iterator it = verts.begin(); it != verts.end(); ++it) {
    (*it).force += forcePerVert;
  }
}
