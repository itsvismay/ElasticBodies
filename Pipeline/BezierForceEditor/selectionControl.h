#pragma once

#include "point.h"

class SelectionControl {
public:
  Point* selectedCtrl;
  bool hasSelection;
  
  SelectionControl();
  void setSelected(Point* point);
  void clearSelected();
};
