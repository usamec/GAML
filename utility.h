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

inline vector<int> InvertPath(const vector<int> x) {
  vector<int> ret;
  if (x.empty()) return ret;
  for (int i = x.size() - 1; i >= 0; i--) {
    if (x[i] >= 0) 
      ret.push_back(x[i] ^ 1);
    else
      ret.push_back(x[i]);
  }
  return ret;
}

inline void ReversePath(vector<int>& path) {
  reverse(path.begin(), path.end());
  for (int i = 0; i < path.size(); i++) {
    if (path[i] >= 0) {
      path[i] ^= 1;
    }
  }
}


#endif
