#pragma once

#include <string>
#include <glm/glm.hpp>

using namespace std;
using namespace glm;

// container for arguements basically

class ProgramSettings {
public:
  string inputMeshObj;
  string inputMeshOff;
  string outputMeshObj;
  string outputMeshOff;
  string config;
  string outputForce;
  string bezierForceFile;

  bool obj;
  bool off;

  bool topForce;
  bool botForce;
  bool rightForce;
  bool leftForce;
  bool frontForce;
  bool backForce;

  bool topFixed;
  bool botFixed;
  bool rightFixed;
  bool leftFixed;
  bool frontFixed;
  bool backFixed;

  bool constantForce;
  bool impulseForce;
  bool bezierForce;

  int corner; // +x+y+z = 0, -x+y+z = 1, -x-y+z = 2, -x+y+z = 3 ...

  double fixedDepth;
  double forceDepth;
  double maxForce;

  //dvec2 boundsX;
  //dvec2 boundsY;
  //dvec2 boundsZ;
  dvec3 forceDirection;

  ProgramSettings() {
    inputMeshObj = "";
    inputMeshOff = "";
    outputMeshObj = "";
    outputMeshOff = "";
    config = "";
    outputForce = "test.txt";
    bezierForceFile = "";

    obj = false;
    off = false;

    topForce = false;
    botForce = false;
    rightForce = false;
    leftForce = false;
    frontForce = false;
    backForce = false;

    topFixed = false;
    botFixed = false;
    rightFixed = false;
    leftFixed = false;
    frontFixed = false;
    backFixed = false;

    constantForce = false;
    impulseForce = false;
    bezierForce = false;

    corner = -1;

    fixedDepth = 0.0;
    forceDepth = 0.0;
    maxForce = 0.0;

    //boundsX = dvec2(0.0, 0.0);
    //boundsY = dvec2(0.0, 0.0);
    //boundsZ = dvec2(0.0, 0.0);
    forceDirection = dvec3(0.0, 0.0, 0.0);
  }

  ~ProgramSettings() { }
};
