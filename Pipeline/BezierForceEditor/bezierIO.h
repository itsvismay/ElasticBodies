#pragma once

#include "bezierOne.h"

class BezierIO {
public:
  static void writeBezierToFile(string file, BezierOne* curve);
  static BezierOne* readBezierFromFile(string file);
}
