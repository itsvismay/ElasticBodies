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
  for (vector<Point*>::iterator it = curve->ctrlPoints.begin(); it != curve->ctrlPoints.end(); ++it) {
    if (mouseX >= (*it)->pos[0] - pointSize / width && mouseX <= (*it)->pos[0] + pointSize / width)
      if (mouseY >= (*it)->pos[1] - pointSize / height && mouseY <= (*it)->pos[1] + pointSize / height) {
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
    selectionControl->selectedCtrl->move(mouseX, mouseY);
    curve->evaluateCtrls();
  }
}

void ProgramData::createCtrl() {
  if (!drag) {
    curve->addCtrl(vec2((float)mouseX, (float)mouseY));
    curve->evaluateCtrls();
  }
}

void ProgramData::deleteCtrl() {
  if (!drag && selectionControl->selectedCtrl) {
    curve->removeCtrl(selectionControl->selectedCtrl);
    curve->evaluateCtrls();
    selectionControl->selectedCtrl = 0x0;
  }
}
