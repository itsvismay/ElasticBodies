#pragma once

#include <vector>

using namespace std;

class CombinationCache {
private:
  static CombinationCache* instance;
  CombinationCache();
  long calculateCombination(long n,long r);
  long findCombination(long n,long r);
  long factorial(long val);
  vector<vector<long>> cache;
public:
  ~CombinationCache();
  static CombinationCache* getInstance();
  static void initialize();
  static void destroy();
  long getCombination(long n,long r);
  void debug();
  void specialComb(long n, long r);
};
