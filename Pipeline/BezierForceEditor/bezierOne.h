#pragma once

// I am representing bezier as its combination form

#include <vector>
#include <glm/glm.hpp>
#include "point.h"

using namespace std;
using namespace glm;

#define VERTS_PER_CTRL 8

class BezierOne {
private:
  vector<vec2> evalPoints;
public:
  vector<Point*> ctrlPoints;
  BezierOne();
  void addCtrl(vec2 point);
  void removeCtrl(Point* point);
  void evaluateCtrls();
  void drawBezier(double pointSize, double lineSize, Point* selected);
  void drawBezierLines();
  void drawBezierPoints();
};
