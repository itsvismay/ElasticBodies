#include "point.h"

Point::Point() {
  pos = vec2(0.0f, 0.0f);
  next = 0x0;
  prev = 0x0;
}

Point::Point(Point* n, Point* p) {
  pos = vec2(0.0f, 0.0f);
  next = n;
  prev = p;
}

Point::Point(vec2 o) {
  pos = o;
  next = 0x0;
  prev = 0x0;
}

void Point::move(double x, double y) {
  if (x > 1.0) x = 1.0;
  if (x < 0.0) x = 0.0;
  if (y > 1.0) y = 1.0;
  if (y < 0.0) y = 0.0;
  if (!next || !prev) {
    pos[1] = (float)y;
  } else {
    pos[1] = (float)y;
    if (x <= prev->pos[0])
      pos[0] = prev->pos[0] + .001f;
    else if (x >= next->pos[0])
      pos[0] = next->pos[0] - .001f;
    else
      pos[0] = x;
  }
}
