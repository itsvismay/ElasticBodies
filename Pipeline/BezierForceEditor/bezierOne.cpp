#include "bezierOne.h"

// this class is a 1D Bezier

BezierOne::BezierOne() {
  ctrlPoints.push_back(vec2(0.0f, 0.5f));
  ctrlPoints.push_back(vec2(0.33f, 0.5f));
  ctrlPoints.push_back(vec2(0.67f, 0.5f));
  ctrlPoints.push_back(vec2(1.0f, 0.5f));
}

void BezierOne::addCtrl(vec2 point) {
  for (vector<vec2>::iterator it = ctrlPoints.begin(), int i = 0; it != ctrlPoints.end(); ++it, ++i) {
    if (point[0] < (*it)[0]) {
      ctrlPoints.insert(ctrlPoints.begin() + i, point);
      it = ctrlPoints.end();
    }
  }
}

void BezierOne::removeCtrl(vec2 point) {
  if (ctrlPoints.size() > 4 && point[0] != 0.0f && point[0] != 1.0f) {
    for (vector<vec2>::iterator it = ctrlPoints.begin(), int i = 0; it != ctrlPoints.end(); ++it, ++i) {
      if ((*it)[0] == point[0] && (*it)[1] == point[1]) {
        ctrlPoints.erase(ctrlPoints.begin() + i);
        it = ctrlPoints.end();
      }
    }
  }
}

void BezierTwo::evaluateCtrls() {
  evalPoints.clear();
  // to be implemented
}
