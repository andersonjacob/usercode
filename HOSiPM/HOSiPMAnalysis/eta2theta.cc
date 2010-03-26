#include <cmath>

double eta2theta( double eta ) {
  return 2.*atan(exp(-eta));
}

double eta2sintheta( double eta ) {
  return sin(eta2theta(eta));
}

double theta2eta( double theta) {
  return log(tan(theta/2.));
}
