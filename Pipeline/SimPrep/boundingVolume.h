#pragma once

#include <vector>
#include <glm/glm.hpp>
#include "fvert.h"

using namespace std;
using namespace glm;

class BoundingVolume {
private:
  // either containingVolumes or verts will be empty
  dvec2 xBounds;
  dvec2 yBounds;
  dvec2 zBounds;
public:
  vector<FVert*> verts;
  BoundingVolume();
  void addVert(FVert* vert);
  void distributeForce(double maxForce);
  void makeFixed();
  dvec3 getNewOrigin(int corner);
};
