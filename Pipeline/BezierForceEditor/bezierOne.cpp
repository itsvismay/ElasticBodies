#include "bezierOne.h"
#include "combinationCache.h"
#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>

// this class is a 1D Bezier

BezierOne::BezierOne() {
  ctrlPoints.push_back(vec2(0.0f, 0.25f));
  ctrlPoints.push_back(vec2(0.33f, 0.1f));
  ctrlPoints.push_back(vec2(0.67f, 0.75f));
  ctrlPoints.push_back(vec2(1.0f, 0.25f));
  evaluateCtrls();
}

void BezierOne::addCtrl(vec2 point) {
  int i = 0;
  for (vector<vec2>::iterator it = ctrlPoints.begin(); it != ctrlPoints.end(); ++it) {
    if (point[0] < (*it)[0]) {
      ctrlPoints.insert(ctrlPoints.begin() + i, point);
      i++;
      it = ctrlPoints.end();
    }
  }
}

void BezierOne::removeCtrl(vec2 point) {
  int i = 0;
  if (ctrlPoints.size() > 4 && point[0] != 0.0f && point[0] != 1.0f) {
    for (vector<vec2>::iterator it = ctrlPoints.begin(); it != ctrlPoints.end(); ++it) {
      if ((*it)[0] == point[0] && (*it)[1] == point[1]) {
        ctrlPoints.erase(ctrlPoints.begin() + i);
        i++;
        it = ctrlPoints.end();
      }
    }
  }
}

void BezierOne::evaluateCtrls() {
  evalPoints.clear();
  int numCtlPts = ctrlPoints.size();
  int numEvalPoints = VERTS_PER_CTRL * numCtlPts; //numCtlPts * 4;
  //float length = ctrlPoints[numCtlPts-1][0] - ctrlPoints[0][0];
  float length = 1.0f;
  float uMult = length / numEvalPoints;
  CombinationCache* cc = CombinationCache::getInstance();
  for(int i=0;i<numEvalPoints+1;i++) {
    float pos = ((float)i)/((float)numEvalPoints);
    float ompos = 1 - pos;
    float ypos = 0.0f;
    float xpos = 0.0f;
    float n = numCtlPts - 1;
    for(int j = 0;j<numCtlPts;j++) {
      float tmp = powf(pos,(float)j) * powf(ompos,(float)n-j);
      tmp *= cc->getCombination(n,j);
      ypos += tmp * ctrlPoints[j][1];
      xpos += tmp * ctrlPoints[j][0];
    }
    evalPoints.push_back(vec2(xpos,ypos));
  }
}

void BezierOne::drawBezier() {
  glLineWidth(4.0);
  glColor4f(0.0f,0.0f,1.0f,1.0f);
  glBegin(GL_LINES);
  drawBezierLines();
  glEnd();
  glPointSize(14.0);
  glColor4f(1.0f,0.0f,0.0f,1.0f);
  glBegin(GL_POINTS);
  drawBezierPoints();
  glEnd();
}

void BezierOne::drawBezierLines() {
  glVertex2f(evalPoints[0][0], evalPoints[0][1]);
  for (vector<vec2>::iterator it = evalPoints.begin() + 1; it != evalPoints.end() - 1; ++it) {
    glVertex2f((*it)[0], (*it)[1]);
    glVertex2f((*it)[0], (*it)[1]);
  }
  glVertex2f(evalPoints[evalPoints.size()-1][0], evalPoints[evalPoints.size()-1][1]);
}

void BezierOne::drawBezierPoints() {
  for (vector<vec2>::iterator it = ctrlPoints.begin(); it != ctrlPoints.end(); ++it) {
    glVertex2f((*it)[0], (*it)[1]);
  }
}
