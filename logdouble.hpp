#ifndef LOG_DOUBLE
#define LOG_DOUBLE

#include <cmath>
#include <algorithm>
#include <iostream>
#include <limits>

using std::min;
using std::max;

class logdouble {
 public:
  double logval;

  logdouble() : logval(-std::numeric_limits<double>::infinity()) {}
  logdouble(double x) : logval(log(x)) {} 

  logdouble operator+=(const logdouble &b) {
    if (isinf(logval) && logval < 0) {
      logval = b.logval;
      return *this;
    }
    if (isinf(b.logval) && b.logval < 0) {
      return *this;
    }
    logval = max(logval, b.logval) + log1p(exp(min(logval, b.logval) - max(logval, b.logval))); 
    return *this;
  }
  logdouble operator*=(const logdouble &b) {
    logval += b.logval;
    return *this;
  }
};

inline logdouble operator+(const logdouble& a, const logdouble& b) {
  if (isinf(a.logval) && a.logval < 0) {
    return b;
  }
  if (isinf(b.logval) && b.logval < 0) {
    return a;
  }
  logdouble ret;
  ret.logval = max(a.logval, b.logval) + log1p(exp(min(a.logval, b.logval) - max(a.logval, b.logval))); 
  return ret;
}

inline logdouble operator*(const logdouble& a, const logdouble& b) {
  logdouble ret;
  ret.logval = a.logval + b.logval;
  return ret;
}

inline logdouble operator^(const logdouble& a, const double& b) {
  logdouble ret;
  ret.logval = a.logval * b;
  return ret;
}

inline logdouble operator/(const logdouble& a, const logdouble& b) {
  logdouble ret;
  ret.logval = a.logval - b.logval;
  return ret;
}

inline std::ostream& operator<<(std::ostream& os, logdouble l) {
  os << l.logval;
  return os;
}

inline bool operator<(const logdouble& a, const logdouble& b) {
  return a.logval < b.logval;
}

inline bool operator>(const logdouble& a, const logdouble& b) {
  return a.logval > b.logval;
}
#endif
