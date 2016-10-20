#include "fvert.h"

FVert::FVert() {
  vert = vec3(0.0f, 0.0f, 0.0f);
  index = -1;
  force = 0.0f;
}

FVert::FVert(vec3 v, int i) {
  vert = v;
  index = i;
  force = 0.0f;
}

void FVert::translate(float x, float y, float z) {
  vert[0] += x;
  vert[1] += y;
  vert[2] += z;
}
