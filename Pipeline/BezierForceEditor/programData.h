#pragma once

#include "bezierOne.h"
#include "selectionControl.h"

class ProgramData {
public:
  BezierOne* curve;
  SelectionControl* selectionControl;
  double width;
  double height;
  double mouseX;
  double mouseY;
  double pointSize;
  double lineSize;
  bool drag;

  ProgramData();
  ~ProgramData();

  void handleMouseDown();
  void handleMouseUp();
  void handleMouseMove(double x,double y);
};
