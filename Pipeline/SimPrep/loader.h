#pragma once

#include <iostream>
#include <fstream>
#include <string.h>
#include <glm/glm.hpp>
#include "mesh.h"

#define MAX_LINE_LENGTH 1024

using namespace std;

class Loader {
public:

  static void loadMeshOBJ(Mesh* mesh, const char* file) {
    ifstream ifs(file, ifstream::in);
    char line[MAX_LINE_LENGTH];

    while (ifs.peek() != EOF) {
      ifs.getline(line, sizeof(line), '\n');
      char* token = line + strspn(line, " \t");

      if (token[0] == 0) continue;
      if (token[0] == 'v' && isSpace(token[1])) parseVert(mesh, token += 2);
      if (token[0] == 'f' && isSpace(token[1])) parseFace(mesh, token += 2);
    }

    ifs.close();
  }

  static void loadMeshOFF(Mesh* mesh, const char* file) {
    ifstream ifs(file, ifstream::in);
    char line[MAX_LINE_LENGTH];

    int* numberOfVerts = new int();
    int* numberOfFaces = new int();
    int* numberOfNorms = new int(); // ignored because not used

    ifs.getline(line, sizeof(line), '\n');
    char* token = line + strspn(line, " \t");

    if (token[0] == 'O' && token[1] == 'F' && token[2] == 'F' && isSpace(token[3])) parseOffData(mesh, token += 4, numberOfVerts, numberOfFaces, numberOfNorms);

    for (int i = 0; i < (*numberOfVerts); i++) {
      ifs.getline(line, sizeof(line), '\n');
      char* token = line + strspn(line, " \t");
      parseVert(mesh, token);
    }

    for (int i = 0; i < (*numberOfFaces); i++) {
      ifs.getline(line, sizeof(line), '\n');
      char* token = line + strspn(line, " \t");
      parseFace(mesh, token += 2);
    }

    delete numberOfVerts;
    delete numberOfFaces;
    delete numberOfNorms;
  }

  static void parseOffData(Mesh* mesh, char* token, int* verts, int* faces, int* norms) {
    dvec3 vert(0.0,0.0,0.0);
    int ind = 0;
    char* tok = strtok(token, " ");

    *verts = atoi(tok);
    tok = strtok(NULL, " ");
    *faces = atoi(tok);
    tok = strtok(NULL, " ");
    *norms = atoi(tok);

    cout << "Verts: " << *verts << " Faces: " << *faces << " Norms: " << *norms << endl;
  }

  static bool isSpace(char token) {
    return token == ' ';
  }

  static void parseVert(Mesh* mesh, char* token) {
    dvec3 vert(0.0,0.0,0.0);
    int ind = 0;
    char* tok = strtok(token, " ");
    while (tok != NULL && ind != 3)
    {
      vert[ind] = atof(tok);
      ind++;
      tok = strtok(NULL, " ");
    }
    mesh->addVert(vert);
  }

  static void parseFace(Mesh* mesh, char* token) {
    ivec3 face(0, 0, 0);
    int ind = 0;
    char* tok = strtok(token, " ");
    while (tok != NULL && ind != 3)
    {
      face[ind] = atoi(tok);
      ind++;
      tok = strtok(NULL, " ");
    }
    mesh->addFace(face);
  }
};
