#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <glm/glm.hpp>
#include "mesh.h"

#define MAX_LINE_LENGTH 1024

using namespace std;

class Loader {
public:

  static void loadMesh(Mesh* mesh, char* file) {
    ifstream ifs = open(file);
    char line[MAX_LINE_LENGTH];

    while (ifs.peek() != EOF) {
      ifs.getline(line, sizeof(line), '\n');
      const char* token = line + strspn(line, " \t");
      
      if (token[0] == 0) continue;
      if (token[0] == 'v' && isSpace(token[1])) parseVert(mesh, token += 2);
      if (token[0] == 'f' && isSpace(token[1])) parseFace(mesh, token += 2);
    }

    ifs.close();
  }

  static bool isSpace(char token) {
    return token == ' ';
  }

  static void parseVert(Mesh* mesh, char* token) {
    // to be implemented
  }

  static void parseFace(Mesh* mesh, char* token) {
    // to be implemented
  }
};
