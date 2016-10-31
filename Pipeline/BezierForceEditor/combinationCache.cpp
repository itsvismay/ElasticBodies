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

long CombinationCache::getCombination(long n,long r) {
  long val = -1;
  if(n < cache.size())
    val = findCombination(n,r);
  if(val==-1)
    return calculateCombination(n,r);
  return val;
}

long CombinationCache::findCombination(long n,long r) {
  return cache[n][r];
}

long CombinationCache::calculateCombination(long n,long r) {
  // first increment the bounds if necessary
  while(cache.size() < n+1) {
    vector<long> newAry;
    for(long i=0;i<cache.size()+1;i++)
      newAry.push_back(-1);
    cache.push_back(newAry);
  }
  long result = factorial(n) / (factorial(r) * factorial(n-r));
  cache[n][r] = result;
  return result;
}

void CombinationCache::specialComb(long n, long r) {
  cout << "factorial n :: " << factorial(n) << endl;
  cout << "factorial r :: " << factorial(r) << endl;
  cout << "factorial n-r :: " << factorial(n-r) << endl;
  cout << "Result :: " << factorial(n) / (factorial(r) * factorial(n-r)) << endl;
}

long CombinationCache::factorial(long val) {
  if(val == 0)
    return 1;
  long product = 1;
  for(long i=1;i<=val;i++)
    product *= i;
  return product;
}

void CombinationCache::debug() {
  cout << "CombinationCache Debug: " << endl;
  cout << "Cache Size: " << cache.size() << endl;
  for(long i=0;i<cache.size();i++)
    cout << "In-Cache Size: " << cache[i].size() << endl;
  cout << "13 choose 2 :: " << calculateCombination(13,2) << endl;
  cout << "20 choose 5 :: " << calculateCombination(20,5) << endl;
  cout << "8 choose 2 :: " << calculateCombination(8, 2) << endl;
  cout << "9 choose 2 :: " << calculateCombination(9, 2) << endl;
  cout << "10 choose 2 :: " << calculateCombination(10, 2) << endl;
  cout << "11 choose 2 :: " << calculateCombination(11, 2) << endl;
  cout << "12 choose 2 :: " << calculateCombination(12, 2) << endl;
  cout << "13 choose 2 :: " << calculateCombination(13, 2) << endl;
  specialComb(13,2);
}
