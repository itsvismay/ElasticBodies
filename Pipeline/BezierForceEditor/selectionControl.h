#pragma once

#include <glm/glm.hpp>

using namespace glm;

class SelectionControl {
public:
  vec2* selectedCtrl;
  bool hasSelection;
  
  SelectionControl();
  void setSelected(vec2* point);
  void clearSelected();
};
