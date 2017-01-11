#pragma once

#include <vector>
#include <string>
#include <glm/glm.hpp>
#include "boundingVolume.h"
#include "programSettings.h"
#include "fvert.h"

using namespace std;
using namespace glm;

class Mesh {
private:
  const char* file;
  vector<ivec3> faces;
  vector<FVert*> verts;
  // should these be properties of mesh or bounding volumes

  void readObjFromFile();
  void readOffFromFile();
public:
  dvec2 xBounds;
  dvec2 yBounds;
  dvec2 zBounds;

  Mesh(string f, ProgramSettings* settings);
  ~Mesh();
  void translate(double x, double y, double z);
  void addVert(dvec3 vert);
  void addFace(ivec3 face);
  void writeObjToFile(ProgramSettings* settings);
  void writeOffToFile(ProgramSettings* settings);
  BoundingVolume* createTopBound(double depth);   // +z
  BoundingVolume* createBotBound(double depth);   // -z
  BoundingVolume* createRightBound(double depth); // +x
  BoundingVolume* createLeftBound(double depth);  // -x
  BoundingVolume* createFrontBound(double depth); // -y
  BoundingVolume* createBackBound(double depth);  // +y
  BoundingVolume* createLineBound(double x, double y, double z, bool alongX, bool alongY, bool alongZ);
  dvec2 xBnds();
  dvec2 yBnds();
  dvec2 zBnds();
};
