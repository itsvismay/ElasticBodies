#pragma once

#include <string>
#include <glm/glm.hpp>

using namespace std;
using namespace glm;

// container for arguements basically

class ProgramSettings {
public:
  string mesh;
  string outputMesh;
  string config;
  string outputForce;
  string bezierForceFile;
  string domainBezierFile;

  // steps
  bool needsAlignment;
  bool cubeAlign;
  bool topCubeAlign;
  bool botCubeAlign;
  bool rightCubeAlign;
  bool leftCubeAlign;
  bool frontCubeAlign;
  bool backCubeAlign;
  bool needsForce;
  bool cubeForce;
  bool topCubeForce;
  bool botCubeForce;
  bool rightCubeForce;
  bool leftCubeForce;
  bool frontCubeForce;
  bool backCubeForce;
  bool domainForce;
  bool domainXTrace;
  bool domainYTrace;
  bool domainZTrace;
  bool domainBezier;
  bool domainPlane;
  bool constantForce;
  bool impulseForce;
  bool bezierForce;

  int corner; // +x+y+z = 0, -x+y+z = 1, -x-y+z = 2, -x+y+z = 3 ...
  float depth;
  float forceDepth;
  float maxForce;
  vec2 domainX;
  vec2 domainY;
  vec2 domainZ;
  vec3 direction;

  ProgramSettings() {
    mesh = "";
    config = "";
    outputForce = "";
    bezierForceFile = "";
    domainBezierFile = "";
    needsAlignment = false;
    cubeAlign = false;
    topCubeAlign = false;
    botCubeAlign = false;
    rightCubeAlign = false;
    leftCubeAlign = false;
    frontCubeAlign = false;
    backCubeAlign = false;
    cubeForce = false;
    topCubeForce = false;
    botCubeForce = false;
    rightCubeForce = false;
    leftCubeForce = false;
    frontCubeForce = false;
    backCubeForce = false;
    domainForce = false;
    domainXTrace = false;
    domainYTrace = false;
    domainZTrace = false;
    domainBezier = false;
    domainPlane = false;
    constantForce = false;
    impulseForce = false;
    bezierForce = false;
    corner = -1;
    depth = 0.0f;
    forceDepth = 0.0f;
    maxForce = 0.0f;
    domainX = vec2(0.0f, 0.0f);
    domainY = vec2(0.0f, 0.0f);
    domainZ = vec2(0.0f, 0.0f);
    direction = vec3(0.0f, 0.0f, 0.0f);
  }

  ~ProgramSettings() { }
};
