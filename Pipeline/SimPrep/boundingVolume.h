#pragma once

#include <vector>
#include <glm/vec4.hpp>

using namespace std;
using namespace glm;

class BoundingVolume {
private:
  // either containingVolumes or verts will be empty
  vector<BoundingVolume*> containingVolumes;
  vector<vec4> verts;
public:
  BoundingVolume();
  ~BoundingVolume();
};
