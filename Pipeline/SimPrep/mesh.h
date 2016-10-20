#pragma once

#include <vector>
#include <string>
#include <glm/glm.hpp>
#include "boundingVolume.h"
#include "fvert.h"

using namespace std;
using namespace glm;

class Mesh {
private:
  char* file;
  vector<ivec4> faces;
  vector<FVert*> verts;
  // should these be properties of mesh or bounding volumes
  vec2 xExtremes;
  vec2 yExtremes;
  vec2 zExtremes;

  void readFromFile();
  void writeToFile();
public:
  Mesh(string f);
  ~Mesh();
  void translate(float x, float y, float z);
  BoundingVolume* createCubeBound(float height);
  BoundingVolume* createTopBound(float depth);   // +z
  BoundingVolume* createBotBound(float depth);   // -Z
  BoundingVolume* createRightBound(float depth); // +x
  BoundingVolume* createLeftBound(float depth);  // -x
  BoundingVolume* createFrontBound(float depth); // -y
  BoundingVolume* createBackBound(float depth);  // +y
};
