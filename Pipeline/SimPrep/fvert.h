#pragma once

#include <glm/glm.hpp>

using namespace glm;

class FVert {
public:
  vec3 vert;
  int index;
  float force;

  FVert();
  FVert(vec3 v, int i);
  void translate(float x, float y, float z);
};