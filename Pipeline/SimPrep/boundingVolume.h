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
  vec2 xExtremes;
  vec2 yExtremes;
  vec2 zExtremes;
public:
  BoundingVolume();
  ~BoundingVolume();
  void addVert(FVert* vert);
  void addBoundingVolume(BoundingVolume* volume);
  void distributeForce(float maxForce);
  vec3 getNewOrigin(int corner);
};
