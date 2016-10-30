// input arguements
// some shit will go here

#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>
#include <iostream>
#include "combinationCache.h"
#include "programData.h"

using namespace std;

// openGL methods
void display(BezierOne *curve);
void reshape(int w,int h);
void keyboard(GLFWwindow* window,int key,int scancode,int action,int mods);
void mouseMove(GLFWwindow* window,double x,double y);
void mouseClick(GLFWwindow* window,int button,int action,int mods);
void error(int error, const char* description);

ProgramData* data;

int main(int argc, char* argv[]) {
  // initialization of program data
  cout << "HELLO WORLD" << endl;

  CombinationCache::initialize();
  data = new ProgramData();

  // initialize glfw
  if(!glfwInit())
    exit(EXIT_FAILURE);
  GLFWwindow* window = glfwCreateWindow(1000,1000,"Continuous Function Editor",NULL,NULL);
  if (!window) {
    glfwTerminate();
    exit(EXIT_FAILURE);
  }
  int width;
  int height;
  glfwSetErrorCallback(error);
  glfwMakeContextCurrent(window);
  glfwSetKeyCallback(window,keyboard);
  glfwSetCursorPosCallback(window,mouseMove);
  glfwSetMouseButtonCallback(window,mouseClick);
  glfwGetFramebufferSize(window, &width, &height);
  glfwSwapInterval(1);

  data->height = height;
  data->width = width;

  // main loop
  while (!glfwWindowShouldClose(window))
  {
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.0,1.0,0.0,1.0,-1.0,1.0);
    glMatrixMode(GL_MODELVIEW);
    glClearColor(0.0f,0.0f,0.0f,1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    display(data->curve);
    glfwSwapBuffers(window);
    glfwPollEvents();
  }
  glfwDestroyWindow(window);
  glfwTerminate();

  delete data;
  return 0;
}

void display(BezierOne *curve) {

  glLineWidth(0.5);
  glBegin(GL_LINES);
  glColor4f(0.7f,0.7f,0.7f,0.5f);
  // draw a box
  glVertex2f(0.999,0.999);
  glVertex2f(0.999,0.001);
  glVertex2f(0.999,0.999);
  glVertex2f(0.001,0.999);
  glVertex2f(0.001,0.001);
  glVertex2f(0.001,0.999);
  glVertex2f(0.001,0.001);
  glVertex2f(0.999,0.001);
  // draw hashes
  for (float i=0.1;i<1.0;i+=0.1) {
    glVertex2f(0.001, i);
    glVertex2f(0.05, i);
    glVertex2f(i, 0.001);
    glVertex2f(i, 0.05);
  }
  for (float i=0.05;i<1.0;i+=0.1) {
    glVertex2f(0.001, i);
    glVertex2f(0.025, i);
    glVertex2f(i, 0.001);
    glVertex2f(i, 0.025);
  }
  for (float i=0.01;i<1.0;i+=0.01) {
    glVertex2f(0.001, i);
    glVertex2f(0.010, i);
    glVertex2f(i, 0.001);
    glVertex2f(i, 0.010);
  }

  curve->drawBezier(data->pointSize, data->lineSize, data->selectionControl->selectedCtrl);

  glEnd();
}

void reshape(int w,int h) {
  glViewport(0,0,w,h);
  data->width = (double)w;
  data->height = (double)h;
}

void keyboard(GLFWwindow* window,int key,int scancode,int action,int mods) {
  if (key == GLFW_KEY_C && action == GLFW_PRESS)
    data->createCtrl();
  if (key == GLFW_KEY_D && action == GLFW_PRESS)
    data->deleteCtrl();
}

void mouseMove(GLFWwindow* window,double x,double y) {
  data->handleMouseMove(x,y);
}

void mouseClick(GLFWwindow* window,int button,int action,int mods) {
  if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
    data->handleMouseDown();
  if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE)
    data->handleMouseUp();
}

void error(int error, const char* description) {
  // to be implemented
}
