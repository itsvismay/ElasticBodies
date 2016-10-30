#include "bezierOne.h"
#include "combinationCache.h"
#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>
#include <iostream>

using namespace std;

// this class is a 1D Bezier

BezierOne::BezierOne() {
  Point* one = new Point(vec2(0.0f, 0.25f));
  Point* two = new Point(vec2(0.33f, 0.1f));
  Point* three = new Point(vec2(0.67f, 0.75f));
  Point* four = new Point(vec2(1.0f, 0.25f));
  one->next = two;
  two->prev = one;
  two->next = three;
  three->prev = two;
  three->next = four;
  four->prev = three;
  ctrlPoints.push_back(one);
  ctrlPoints.push_back(two);
  ctrlPoints.push_back(three);
  ctrlPoints.push_back(four);
  evaluateCtrls();
}

void BezierOne::addCtrl(vec2 point) {
  int i = 0;
  for (vector<Point*>::iterator it = ctrlPoints.begin(); it != ctrlPoints.end(); ++it) {
    if (point[0] < (*it)->pos[0] && point[0] > 0.0f && point[0] < 1.0f) {
      Point* newPt = new Point(point);
      newPt->prev = *(it - 1);
      newPt->next = *(it);
      (*(it - 1))->next = newPt;
      (*(it))->prev = newPt;
      ctrlPoints.insert(ctrlPoints.begin() + i, newPt);
      it = ctrlPoints.end() - 1;
    }
    i++;
  }
}

void BezierOne::removeCtrl(Point* point) {
  int i = 0;
  if (ctrlPoints.size() > 4 && point->pos[0] != 0.0f && point->pos[0] != 1.0f) {
    for (vector<Point*>::iterator it = ctrlPoints.begin(); it != ctrlPoints.end(); ++it) {
      if (*it == point) {
        (*it)->prev->next = (*it)->next;
        (*it)->next->prev = (*it)->prev;
        ctrlPoints.erase(ctrlPoints.begin() + i);
        it = ctrlPoints.end() - 1;
      }
      i++;
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
      ypos += tmp * ctrlPoints[j]->pos[1];
      xpos += tmp * ctrlPoints[j]->pos[0];
    }
    evalPoints.push_back(vec2(xpos,ypos));
  }
}

void BezierOne::drawBezier(double pointSize, double lineSize, Point* selected) {
  glLineWidth(lineSize);
  glColor4f(0.0f,0.0f,1.0f,1.0f);
  glBegin(GL_LINES);
  drawBezierLines();
  glEnd();
  glPointSize(pointSize);
  glColor4f(1.0f,0.0f,0.0f,1.0f);
  glBegin(GL_POINTS);
  drawBezierPoints();
  glEnd();
  if (selected) {
    glColor4f(0.0f,1.0f,0.0f,1.0f);
    glBegin(GL_POINTS);
    glVertex2f(selected->pos[0], selected->pos[1]);
    glEnd();
  }
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
  for (vector<Point*>::iterator it = ctrlPoints.begin(); it != ctrlPoints.end(); ++it) {
    glVertex2f((*it)->pos[0], (*it)->pos[1]);
  }
}
