/*
  model of a SiPM device that accurately depicts the response to a
  photon impulse.  The photon impulse is assumed to 100% efficient in
  potentially producing an avalanche in a micro-pixel.  The assumption
  is that all of the pixels in the sensor are uniformly illuminated.
  The device is described by a few parameters.
  
  NP : the number of micro-pixels in the sensor

  crossTalk : probability that a single photo-electron will trigger an 
    additional photo-electron
  tempDep : temperature dependence of the gain as a percent of nominal
    per unit temperature.

  For most SiPMs the micro-pixels recover as an RC circuit.  This can be
    modeled using the tau parameter.  Generally one would use either the
    RC based recovery model or the one described below, but not both
    simultanously
  tau : the micro-pixel characteristic recharge time (RC time constant)


  For Zecotek sensors the recharge time is characterized by two
    periods with a knee between them.  The next 3 parameters capture
    this behaviour.  
  rTime1 : This is the time of the knee.
  rPct1 : this is the fractional charge of a pixel at time rTime1 
  rTime2 : This is the time when a pixel is full recharged.

 */

#include "SiPMModel.h"
#include "TMath.h"
#include <assert.h>

TRandom3 SiPMModel::rnd = TRandom3(0);

SiPMModel::SiPMModel(unsigned int NP, double rTime1, double rPct1, 
		     double rTime2, double tau) : 
  _NP(NP), _SiPM(NP), _rTime1(rTime1), _recoveryPct1(rPct1), _rTime2(rTime2),
  _tau(tau), _crossTalk(0), _tempDep(0)
{
  if (tau <= 0) _tau = rTime1;
  resetSiPM();
}

SiPMModel::SiPMModel() : _NP(0), _rTime1(0), _recoveryPct1(0), _rTime2(0), _tau(0), _crossTalk(0), _tempDep(0)
{ 
}

SiPMModel::SiPMModel(unsigned int NP, double tau) :
  _NP(NP), _SiPM(NP), _rTime1(tau), _recoveryPct1(0.6), _rTime2(10*tau),
  _tau(tau), _crossTalk(0), _tempDep(0)
{
  resetSiPM();
}

SiPMModel::SiPMModel(SiPMModel const& other) :
  _NP(other._NP), _SiPM(other._SiPM), _rTime1(other._rTime1), 
  _recoveryPct1(other._recoveryPct1), _rTime2(other._rTime2), _tau(other._tau),
  _crossTalk(other._crossTalk), _tempDep(other._tempDep)
{
}

double SiPMModel::hitPixels(unsigned int pes, double fraction, 
			    double tempDiff) {
  // response to light impulse with pes input photons.  The return is the number
  // of micro-pixels hit.  If a fraction other than 0. is supplied then the
  // micro-pixel doesn't fully discharge.  The tempDiff is the temperature 
  // difference from nominal and is used to modify the relative strength of a
  // hit pixel.  Pixels which are fractionally charged return a fractional
  // number of hit pixels.

  unsigned int pixel;
  double sum = 0.;
  if ((_crossTalk > 0.) && (_crossTalk < 1.)) {
    double realcross = pes/(1. - _crossTalk) - pes;
    pes += rnd.Poisson(realcross);
  }

  for (unsigned int pe = 0; pe < pes; ++pe) {
    pixel = rnd.Integer(_NP);
    sum += (_SiPM[pixel]*(1 + (tempDiff*_tempDep)));
    _SiPM[pixel] = fraction;
  }

  return sum;
}

double SiPMModel::recoveryModel(double t) {
  // the Zecotek recover model as a function of t.
  if (t < 0) return 1.0;
  if (t == 0) return 0.;
  double lt = TMath::Log10(t);
  double lrt1 = TMath::Log10(_rTime1);
  double tau = _recoveryPct1/lrt1;
  if (t < _rTime1) return lt*tau;
  double lrt2 = TMath::Log10(_rTime2);
  tau = (1.-_recoveryPct1)/(lrt2-lrt1);
  double val = tau*(lt-lrt1)+_recoveryPct1;
  return ((val<1.) ? val : 1.);
}

double SiPMModel::totalCharge() const {
  // sum of the micro-pixels.  NP is a fully charged device.
  // 0 is a fullly depleted device.
  double tot = 0.;
  for(unsigned int i=0; i<_NP; ++i) tot += _SiPM[i];
  return tot;
}

void SiPMModel::expRecover(double dt) {
  // recover each micro-pixel using the RC model.  For this to work well.
  // dt << tau (typically dt = 0.2*tau or less)
  double newval;
  for (unsigned int i=0; i<_NP; ++i) {
    newval = _SiPM[i] + (1 - _SiPM[i])*dt/_tau;
    _SiPM[i] = (newval < 1.0) ? newval : 1.0;
  }
}

void SiPMModel::recoverForTime(double time, double dt) {
  // apply the RC recover model to the pixels for time.  If dt is not
  // positive then tau/5 will be used for dt.
  if (dt <= 0.)
    dt = _tau/5.;
  for (double t = 0; t <= time; t += dt)
    expRecover(dt);
}

double SiPMModel::hitPixelsQ(unsigned int NP, unsigned int pes, double eff,
			     double tempDep, double dT, double xtalk) {

  // This is a quick, analytic, static model of a sipm.  It doesn't keep a 
  // pixel history.  All parameters are passed in and the number of hit
  // cells is returned.  If eff < 1. then the effective number of cells is
  // reduced by that fraction making some dead space.  This implementation
  // can return a value for an arbitrarily large number of pes with only
  // generating 2 random numbers.

  if (pes < 1) return 0;
  assert((eff >= 0.)&&(eff<=1.));
  if ((xtalk > 0.) && (xtalk < 1.))
    pes += rnd.Poisson(pes/(1-xtalk)-pes);
  double alpha = pes/double(NP);
  double interp = 1. - TMath::Exp(-alpha);
  if (interp > 1.) interp = 1.;
  interp *= eff*NP;
  double sig2 = NP*TMath::Exp(-alpha)*(1-(1+alpha)*TMath::Exp(-alpha));
  unsigned int q;
  while (true) {
    q = static_cast<unsigned int>(rnd.Gaus(interp, TMath::Sqrt(sig2+interp*(1-eff)))+0.5);
    if ((q<=NP))
      return q*(1-(dT*tempDep));
  }
}

void SiPMModel::setNP(unsigned int NP) {
  // change the number of micro-pixels in the sensor.

  _NP = NP;
  _SiPM.resize(_NP);
  resetSiPM();
}

void SiPMModel::setCrossTalk(double xTalk) {
  // set the cross-talk probability

  if(xTalk < 0 || xTalk > 1) {
    _crossTalk = 0.;
  } else {
    _crossTalk = xTalk;
  }   

}
void SiPMModel::setTempDep(double dTemp){
  // set the temperature dependence

    _tempDep = dTemp;
}
