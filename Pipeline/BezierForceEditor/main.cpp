#pragma once

// input arguements
// some shit will go here

#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>
#include "bezierOne.h"

// openGL methods
void display();
void reshape(int w,int h);
void keyboard(GLFWwindow* window,int key,int scancode,int action,int mods);
void mouseMove(GLFWwindow* window,double x,double y);
void mouseClick(GLFWwindow* window,int button,int action,int mods);
void error(int error, const char* description);

int main(int argc, char* argv[]) {
  // initialization of program data
  BezierOne *curve = new BezierOne();

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
  glViewport(0,0,1.0f,1.0f);
  glfwSwapInterval(1);
  transitionToMapView();

  // main loop
  while (!glfwWindowShouldClose(window))
  {
    dt+=.1;
    display(curve);
    glfwSwapBuffers(window);
    glfwPollEvents();
  }
  glfwDestroyWindow(window);
  glfwTerminate();
  return 0;
}

void display(BezierOne *curve) {
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glMatrixMode(GL_MODELVIEW);
  glClearColor(0.0f,0.0f,0.0f,1.0f);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();

  glLineWidth(0.5);
  glBegin(GL_LINES);
  glColor4f(0.7f,0.7f,0.7f,0.5f);
  // draw a box
  glVertex2f(0.5,0.5);
  glVertex2f(0.5,-0.5);
  glVertex2f(0.5,0.5);
  glVertex2f(-0.5,0.5);
  glVertex2f(-0.5,-0.5);
  glVertex2f(-0.5,0.5);
  glVertex2f(-0.5,-0.5);
  glVertex2f(0.5,-0.5);

  //curve->drawBezier();

  glEnd();
}

void reshape(int w,int h) {
  // to be implemented
}

void keyboard(GLFWwindow* window,int key,int scancode,int action,int mods) {
  // to be implemented
}

void mouseMove(GLFWwindow* window,double x,double y) {
  // to be implemented
}

void mouseClick(GLFWwindow* window,int button,int action,int mods) {
  // to be implemented
}

void error(int error, const char* description) {
  // to be implemented
}
