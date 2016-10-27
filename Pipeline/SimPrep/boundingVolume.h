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
  dvec2 xExtremes;
  dvec2 yExtremes;
  dvec2 zExtremes;
public:
  vector<FVert*> verts;
  BoundingVolume();
  ~BoundingVolume();
  void addVert(FVert* vert);
  void addBoundingVolume(BoundingVolume* volume);
  void distributeForce(double maxForce);
  void makeFixed();
  vec3 getNewOrigin(int corner);
};
