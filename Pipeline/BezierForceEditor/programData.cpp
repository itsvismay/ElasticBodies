#include "programData.h"
#include <iostream>

using namespace std;

ProgramData::ProgramData() {
  curve = new BezierOne();
  selectionControl = new SelectionControl();
  selectionControl->selectedCtrl = 0x0;
  mouseX = 0.0;
  mouseY = 0.0;
  width = 0.0;
  height = 0.0;
  lineSize = 4.0;
  pointSize = 14.0;
  drag = false;
  selectionControl->hasSelection = false;
}

ProgramData::~ProgramData() {
  delete curve;
  delete selectionControl;
}

void ProgramData::handleMouseDown() {
  for (vector<vec2*>::iterator it = curve->ctrlPoints.begin(); it != curve->ctrlPoints.end(); ++it) {
    if (mouseX >= (**it)[0] - pointSize / width && mouseX <= (**it)[0] + pointSize / width)
      if (mouseY >= (**it)[1] - pointSize / height && mouseY <= (**it)[1] + pointSize / height) {
        selectionControl->setSelected(*it);
        drag = true;
      }
  }
}

void ProgramData::handleMouseUp() {
  drag = false;
}

void ProgramData::handleMouseMove(double x,double y) {
  mouseX = x / width;
  mouseY = 1.0 - (y / height);
  //cout << x << "," << y << endl;
  if (drag && selectionControl->hasSelection) {
    (*(selectionControl->selectedCtrl))[0] = mouseX;
    (*(selectionControl->selectedCtrl))[1] = mouseY;
    curve->evaluateCtrls();
  }
}
