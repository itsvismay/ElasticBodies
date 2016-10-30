#include "selectionControl.h"

SelectionControl::SelectionControl() {
  selectedCtrl = 0x0;
  hasSelection = false;
}

void SelectionControl::setSelected(vec2* point) {
  selectedCtrl = point;
  hasSelection = true;
}

void SelectionControl::clearSelected() {
  hasSelection = false;
}
