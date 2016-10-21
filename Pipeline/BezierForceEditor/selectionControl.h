#pragma once

#include <glm/glm.hpp>

using namespace glm;

class SelectionControl {
private:
  vec2 selectedCtrl;
  bool hasSelection;
public:
  SelectionControl();
  void setSelected(vec2 point);
  void clearSelected();
};
