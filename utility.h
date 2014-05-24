#ifndef UTILITY_H__
#define UTILITY_H__

#include <vector>
#include <string>
#include <cstdlib>
#include <unordered_set>

using namespace std;

inline int StringToInt(string x) {
  return atoi(x.c_str());
}

inline double StringToDouble(string x) {
  return atof(x.c_str());
}

template<class T>
vector<T> USetToVector(const unordered_set<T>& x) {
  vector<T> ret;
  for (auto& y: x) {
    ret.push_back(y);
  }
  return ret;
}



#endif
