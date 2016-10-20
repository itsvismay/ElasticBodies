#pragma once

#include <vector>
#include <glm/glm.hpp>
#include "boundingVolume.h"

using namespace std;
using namespace glm;

class Mesh {
private:
  vector<ivec4> faces;
  vector<vec4> verts;
  vec2 xExtremes;
  vec2 yExtremes;
  vec2 zExtremes;

  void readFromFile();
public:
  Mesh(char* file);
  BoundingVolume* createTopBound(float depth);
  BoundingVolume* createBotBound(float depth);
  BoundingVolume* createRightBound(float depth);
  BoundingVolume* createLeftBound(float depth);
  BoundingVolume* createFrontBound(float depth);
  BoundingVolume* createBackBound(float depth);
};
