// -*- c++ -*-
#include <vector>
#include "TRandom3.h"

#ifndef SIPMMODEL_H
#define SIPMMODEL_H

class SiPMModel {
public:
  // SiPMModel(unsigned int NP, double rTime1, double rPct1, double rTime2, 
  // 	    double tau = -1.0);
  SiPMModel();
  SiPMModel(unsigned int NP, double tau);
  SiPMModel(SiPMModel const& other);
  virtual ~SiPMModel() { }

  void resetSiPM() { std::fill(_SiPM.begin(), _SiPM.end(), -999.); }
  virtual double hitPixels(unsigned int pes, double tempDiff =0.,
			   double photonTime = 0.);
  // virtual double recoveryModel(double t);
  virtual double totalCharge() const { return totalCharge(theLastHitTime); }
  virtual double totalCharge(double time) const;
  // virtual void expRecover(double dt);
  // virtual void recoverForTime(double time, double dt = 0.);

  unsigned int getNP() const { return _NP; }
  // double getrTime1() const { return _rTime1; }
  // double getrPct1() const { return _recoveryPct1; }
  // double getrTime2() const { return _rTime2; }
  double getTau() const { return _tau; }
  double getCrossTalk() const { return _crossTalk; }
  double getTempDep() const { return _tempDep; }

  void setNP(unsigned int NP);
  void setTau(double tau) { _tau = tau; }
  void setCrossTalk(double xtalk);
  void setTempDep(double dTemp);

  static double hitPixelsQ(unsigned int NP, unsigned int pes, double eff = 1.0,
			   double tempDep = 0., double dT = 0., 
			   double xtalk = 0.);
  static void setSeed(unsigned int seed) { rnd.SetSeed(seed); }
  static double correctSaturation(double output, unsigned int NP, 
				  float xtalk = 0.);

protected:

  double cellCharge(double deltaTime) const;

  unsigned int _NP;
  std::vector<double> _SiPM;
  double _tau;
  double _crossTalk;
  double _tempDep;
  double theLastHitTime;

  static TRandom3 rnd;

};

#endif
