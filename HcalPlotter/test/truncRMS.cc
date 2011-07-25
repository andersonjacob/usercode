#include <vector>
#include <map>

#include "TString.h"
#include "TMinuit.h"
#include "TFitter.h"
#include "TF1.h"

//double cutoff;
unsigned int Ndata = 1000000;
double * BarrelE = 0;
double * EBE = 0;
double * HOE = 0;
double lastRms = 10.0;
std::map< TString, double * > data;
TFitter const * theFitter = 0;

//double rangeMin = 0;
//double rangeMax = 500;

void truncRMS(int &npar, double */*gin*/, double &retval, double *parArray,
	      int /*iflag*/) {
  //std::cout << npar << '\n';
  assert(data.size() > 0);
  assert(npar>=1);
  double sum2 = 0/*, sumw = 0*/;
  double Esum = 0.;
  sum2 = 0.;
  unsigned int N = 0;
  double MPV = parArray[0];
  std::map< TString, double * >::iterator row;
  unsigned int rowN = 1;
  for (unsigned int datum = 0; datum<Ndata; ++datum) {
    Esum = 0.;
    for (rowN = 0; rowN < data.size() ; ++rowN) {
      row = data.find(theFitter->GetParName(rowN+1));
      if (row != data.end())
	Esum += row->second[datum]*parArray[rowN+1];
    }
    ++N;
    sum2 += ((Esum-MPV)*(Esum-MPV));
  }
  lastRms = sqrt(sum2/N);
  retval = sqrt(sum2/N)/MPV;
  //std::cout << "trunc rms: " << retval << '\n';
}

// void setCutoff(double cut) {
//   cutoff = cut;
// }

// void setRange(double min, double max) {
//   rangeMin = min;
//   rangeMax = max;
// }

void setData(unsigned int Nd, double * datapts, TString dataName) {
  if (Nd < Ndata)
    Ndata = Nd;
  std::map< TString, double * >::iterator pts = data.find(dataName);
  double * currData = 0;
  if (pts != data.end()) {
    delete[] pts->second;
    pts->second = new double[Ndata];
    currData = pts->second;
  } else {
    currData = new double[Ndata];
    pts = data.insert(std::pair< TString, double * >(dataName, currData)).first;
  }
  for (unsigned int i=0; i<Ndata; ++i) {
    currData[i] = datapts[i];
  }
  std::cout << dataName << " data added " << Ndata << " pts.\n"
	    << "data size: " << data.size() << '\n';
}

void hookupMinuit(TFitter &min) {
  min.SetFCN(truncRMS);
  theFitter = &min;
}

double truncFuncHO(double * x, double * par) {
  double rmsPars[3] = {0};
  rmsPars[0] = x[0];
  rmsPars[1] = par[0];
  rmsPars[2] = par[1];
  double retval = -1;
  double gin[2] = {0};
  int iflag = 0;
  int nPars = 3;
  truncRMS(nPars, gin, retval, rmsPars, iflag);
  return retval;
}

TF1 * makeFuncHO() {
  return new TF1("hoFunc", truncFuncHO, 0.0, 2.5, 2);
}
