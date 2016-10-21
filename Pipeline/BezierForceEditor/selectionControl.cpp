#include "selectionControl.h"

SelectionControl::SelectionControl() {
  selectedCtrl = vec2(0.0f, 0.0f);
  hasSelection = false;
}

void SelectionControl::setSelected(vec2 point) {
  selectedCtrl = point;
  hasSelection = true;
}

void SelectionControl::clearSelected() {
  hasSelection = false;
}
