#include "bezierOne.h"
#include "combinationCache.h"
#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>
#include <glm/gtx/spline.hpp>
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
  for (int i=0;i<ctrlPoints.size();i++)
    cout << "XPOS :: " << ctrlPoints[i]->pos[0] << endl;
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

void BezierOne::evaluateCtrlsTwo() {
  evalPoints.clear();
  int numCtrlPts = ctrlPoints.size();
  int numEvalPoints = VERTS_PER_CTRL * numCtrlPts;
  for (int i = 0; i < numEvalPoints+1; i++) {
    float pointPos = (float)i / (float)numEvalPoints;
    vec2 point = getValueAtPoint(pointPos, (int)((float)(i-1) / (float)VERTS_PER_CTRL / 4.0f));
    cout << "POIHNT :: " << pointPos << endl;
    evalPoints.push_back(point);
    //int ind = i / VERTS_PER_CTRL;
    //int i0 = ind == 0 ? 0 : ind-1;
    //int i1 = ind;
    //int i2 = ind == numCtrlPoints ? numCtrlPoints : ind+1;
    //int i3 = ind == numCtrlPoints ? numCtrlPoints : ind+2;
    //float fractPoint = (ctrlPoints[i3]->pos[0] - ctrlPoints[i0]->pos[0])
    //vec2 point = glm::gtx::spline::cubic(ctrlPoints[i0], ctrlPoints[i1], ctrlPoints[i2], ctrlPoints[i3], )
    
  }
}

void BezierOne::subdivide() {
  // to be implemented
  cout << "WHAT" << endl;
  int numCtrls = ctrlPoints.size();
  float currentSpacing = 1.0f / numCtrls;
  float newSpacing = currentSpacing / 4.0f;
  vector<vec2> newPts = vector<vec2>();
  int ind = 0;
  int numMade = 0;
  for(float f = 0.0f; f < 1.0f; f += newSpacing) {
    if (fract(f / currentSpacing) != 0.0f)
    {
      cout << "MAKING ONE" << endl;
      newPts.push_back(getValueAtPoint(f, ind));
      cout << "MADE ONE" << endl;
      numMade++;
    }
    else {
      //ind++;
    }
  }
  cout << "afasdfasfas" << endl;
  for (int i = 0; i < newPts.size(); i++) {
    addCtrl(newPts[i]);
  }
  cout << "Subdivision Created :: " << numMade << " Points " << endl;
  evaluateCtrls();
}

vec2 BezierOne::getValueAtPoint(float point, int ind) {
  //int i0 = clamp<int>(point - 1, 0, ctrlPoints.size() - 1);
  //int i1 = clamp<int>(point, 0, ctrlPoints.size() - 1);
  //int i2 = clamp<int>(point + 1, 0, ctrlPoints.size() - 1);
  //int i3 = clamp<int>(point + 2, 0, ctrlPoints.size() - 1);
  int i0 = ind == 0 ? 0 : ind * 4 - 1;
  float partial = fract((point * VERTS_PER_CTRL) / 4.0f);
  

  double pos = point;
  double ompos = 1.0f - point;
  double xpos = 0.0f;
  double ypos = 0.0f;
  int n = 3;
  //cout << "IND :: " << ind << endl;
  CombinationCache* cc = CombinationCache::getInstance();
  for(int j = 0;j<4;j++) {
    double tmp = pow(pos,(double)j) * powf(ompos,(double)n-j);
    tmp *= cc->getCombination(n,j);
    //cout << "j + i0 :: " << j + i0 << endl;
    ypos += tmp * ctrlPoints[j+i0]->pos[1];
    xpos += tmp * ctrlPoints[j+i0]->pos[0];
  }
  return vec2((float)xpos, (float)ypos);
  //  evalPoints.push_back(vec2(xpos,ypos));



  //cout << "PARTIAL :: " << partial << endl;
  //vec2 result = cubic(ctrlPoints[i0]->pos, ctrlPoints[i1]->pos, ctrlPoints[i2]->pos, ctrlPoints[i3]->pos, point); 
  //cout << "RESULT :: " << result[0] << "," << result[1] << endl;
  //result[0] -= 1;
  //result[1] /= 4;
  //return result;
  //return (((ctrlPoints[i3]->pos * point) + ctrlPoints[i2]->pos) * point + ctrlPoints[i1]->pos) * point + ctrlPoints[i0]->pos;
}

void BezierOne::evaluateCtrls() {
  evalPoints.clear();
  int numCtlPts = ctrlPoints.size();
  int numEvalPoints = VERTS_PER_CTRL * numCtlPts; //numCtlPts * 4;
  //float length = ctrlPoints[numCtlPts-1][0] - ctrlPoints[0][0];
  //float length = 1.0f;
  //float uMult = length / numEvalPoints;
  CombinationCache* cc = CombinationCache::getInstance();
  for(int i=0;i<numEvalPoints+1;i++) {
    double pos = ((float)i)/((float)numEvalPoints);
    double ompos = 1 - pos;
    double ypos = 0.0f;
    double xpos = 0.0f;
    double n = numCtlPts - 1;
    for(int j = 0;j<numCtlPts;j++) {
      double tmp = powf(pos,(double)j) * powf(ompos,(double)n-j);
      tmp *= cc->getCombination(n,j);
      ypos += tmp * ctrlPoints[j]->pos[1];
      xpos += tmp * ctrlPoints[j]->pos[0];
    }
    evalPoints.push_back(vec2((float)xpos,(float)ypos));
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
    glColor4f(0.0f, 0.4f, 0.0f, 1.0f);
    glBegin(GL_LINES);
    glVertex2f(0.0f, selected->pos[1]);
    glVertex2f(selected->pos[0], selected->pos[1]);
    glVertex2f(selected->pos[0], 0.0f);
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
