#pragma once

#include <glm/glm.hpp>

using namespace glm;

class Point {
public:
  vec2 pos;
  Point* next;
  Point* prev;

  Point();
  Point(Point* n, Point* p);
  Point(vec2 o);

  void move(double x, double y);
};
