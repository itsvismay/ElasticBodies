#pragma once

// I am representing bezier as its combination form

#include <vector>
#include <glm/glm.hpp>

using namespace std;
using namespace glm;

#define VERTS_PER_CTRL 8

class BezierOne {
private:
  vector<vec2> ctrlPoints;
  vector<vec2> evalPoints;
public:
  BezierOne();
  void addCtrl(vec2 point);
  void removeCtrl(vec2 point);
  void evaluateCtrls();
  void drawBezier();
  void drawBezierLines();
  void drawBezierPoints();
};
