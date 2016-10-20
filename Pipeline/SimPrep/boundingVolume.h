#pragma once

#include <vector>
#include <glm/glm.hpp>
#include "fvert.h"

using namespace std;
using namespace glm;

class BoundingVolume {
private:
  // either containingVolumes or verts will be empty
  vector<BoundingVolume*> containingVolumes;
  vector<FVert*> verts;
public:
  BoundingVolume();
  ~BoundingVolume();
  void addVert(FVert* vert);
  void addBoundingVolume(BoundingVolume* volume);
  void distributeForce(float maxForce);
};
