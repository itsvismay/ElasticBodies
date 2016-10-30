#pragma once

// I am representing bezier as its combination form

#include <vector>
#include <glm/glm.hpp>

using namespace std;
using namespace glm;

#define VERTS_PER_CTRL 8

class BezierOne {
private:
  vector<vec2> evalPoints;
public:
  vector<vec2*> ctrlPoints;
  BezierOne();
  void addCtrl(vec2* point);
  void removeCtrl(vec2* point);
  void evaluateCtrls();
  void drawBezier(double pointSize, double lineSize, vec2* selected);
  void drawBezierLines();
  void drawBezierPoints();
};
