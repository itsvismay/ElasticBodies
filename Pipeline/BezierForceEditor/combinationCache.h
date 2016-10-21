#pragma once

#include <vector>

using namespace std;

class CombinationCache {
private:
  static CombinationCache* instance;
  CombinationCache();
  int calculateCombination(int n,int r);
  int findCombination(int n,int r);
  int factorial(int val);
  vector<vector<int>> cache;
  void debug();
public:
  ~CombinationCache();
  static CombinationCache* getInstance();
  static void initialize();
  static void destroy();
  int getCombination(int n,int r);
};
