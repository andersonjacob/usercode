#include "TH1.h"
#include "TMinuit.h"
#include "TFitter.h"
#include "TF1.h"

//double cutoff;
unsigned int Ndata = 0;
double * BarrelE = 0;
double * EBE = 0;
double * HOE = 0;
double lastRms = 10.0;
//double rangeMin = 0;
//double rangeMax = 500;

void truncRMS(int &npar, double */*gin*/, double &retval, double *parArray,
	      int /*iflag*/) {
  //std::cout << npar << '\n';
  assert(BarrelE);
  assert(HOE);
  assert(npar>=1);
  double sum = 0/*, sumw = 0*/;
  double Esum;
  // TH1D dhist("dhist", "dhist", 100, rangeMin, rangeMax);
  // for (unsigned int datum = 0; datum<Ndata; ++datum) {
  //   Esum = BarrelE[datum] + parArray[0]*HOE[datum] + parArray[1]*EBE[datum];
  //   dhist.Fill(Esum);
  // }
  // for (int i=dhist.GetMaximumBin()-2; i <= dhist.GetMaximumBin()+2; ++i) {
  //   sum += dhist.GetBinCenter(i)*dhist.GetBinContent(i);
  //   sumw += dhist.GetBinContent(i);
  // }
  // double mpv = sum/sumw;
  //printf("mpv: %f\n", mpv);
  sum = 0.;
  unsigned int N = 0;
  double MPV = parArray[2];
  for (unsigned int datum = 0; datum<Ndata; ++datum) {
    Esum = BarrelE[datum] + parArray[0]*HOE[datum] + parArray[1]*EBE[datum];
    //    if ( ((Esum - MPV) > -2.*lastRms)/* && 
    //	 ((Esum - MPV) <  5.*cutoff*MPV)*/ ) {
      ++N;
      sum += ((Esum-MPV)*(Esum-MPV));
    //    }
  }
  lastRms = sqrt(sum/N);
  retval = sqrt(sum/N)/MPV;
}

// void setCutoff(double cut) {
//   cutoff = cut;
// }

// void setRange(double min, double max) {
//   rangeMin = min;
//   rangeMax = max;
// }

void setData(unsigned int Nd, double * BE, double * OE, double * EE) {
  Ndata = Nd;
  if (BarrelE) delete[] BarrelE;
  if (HOE) delete[] HOE;
  if (EBE) delete[] EBE;
  BarrelE = new double[Ndata];
  HOE = new double[Ndata];
  EBE = new double[Ndata];
  for (unsigned int i=0; i<Ndata; ++i) {
    BarrelE[i] = BE[i];
    HOE[i] = OE[i];
    EBE[i] = EE[i];
  }
}

void hookupMinuit(TFitter &min) {
  min.SetFCN(truncRMS);
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
