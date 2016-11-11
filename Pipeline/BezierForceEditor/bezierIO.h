#pragma once

#include "bezierOne.h"
#include <string>

using namespace std;

class BezierIO {
public:
  static void writeBezierToFile(string file, BezierOne* curve);
  static BezierOne* readBezierFromFile(string file);
};
