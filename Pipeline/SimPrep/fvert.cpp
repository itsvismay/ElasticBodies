#include "fvert.h"

FVert::FVert() {
  vert = dvec3(0.0, 0.0, 0.0);
  index = -1;
  force = 0.0;
  fixed = false;
}

FVert::FVert(dvec3 v, int i) {
  vert = v;
  index = i;
  force = 0.0;
  fixed = false;
}

void FVert::translate(double x, double y, double z) {
  vert[0] += x;
  vert[1] += y;
  vert[2] += z;
}
