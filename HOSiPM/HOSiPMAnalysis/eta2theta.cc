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

double phidiff( double phi1, double phi2, double maxphi ) {
  double dif1 = fabs(phi1-phi2);
  double dif2 = fabs(phi1+maxphi-phi2);
  double dif3 = fabs(phi1-maxphi-phi2);
  if ((dif1<dif2) && (dif1<dif3))
    return dif1;
  else if (dif2<dif3)
    return dif2;
  else
    return dif3;
  return 0.;
}
