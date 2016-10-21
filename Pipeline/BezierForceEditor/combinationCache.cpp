#include "combinationCache.h"
#include <iostream>


CombinationCache* CombinationCache::instance = 0x0;

CombinationCache::CombinationCache() {
  // not really needed
}

CombinationCache::~CombinationCache() {
  // not really needed
}

CombinationCache* CombinationCache::getInstance() {
  return instance;
}

void CombinationCache::initialize() {
  if(!instance)
    instance = new CombinationCache();
}

void CombinationCache::destroy() {
  if(instance)
    delete instance;
}

int CombinationCache::getCombination(int n,int r) {
  int val = -1;
  if(n < cache.size())
    val = findCombination(n,r);
  if(val==-1)
    return calculateCombination(n,r);
  return val;
}

int CombinationCache::findCombination(int n,int r) {
  return cache[n][r].val;
}

int CombinationCache::calculateCombination(int n,int r) {
  // first increment the bounds if necessary
  while(cache->getSize() < n+1) {
    vector<int> newAry;
    for(int i=0;i<cache->size()+1;i++)
      newAry.push_back(-1);
    cache.push_back(newAry);
  }
  int result = factorial(n) / (factorial(r) * factorial(n-r));
  cache[n][r] = result;
  return result;
}

int CombinationCache::factorial(int val) {
  if(val == 0)
    return 1;
  int product = 1;
  for(int i=1;i<=val;i++)
    product *= i;
  return product;
}

void CombinationCache::debug() {
  cout << "CombinationCache Debug: " << endl;
  cout << "Cache Size: " << cache->getSize() << endl;
  for(int i=0;i<cache.size();i++)
    cout << "In-Cache Size: " << cache[i].size() << endl;
}
