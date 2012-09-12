#include <iostream>
#include <vector>
#include "TRandom3.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TMath.h"
#include "TGraph2D.h"

#include "QIE8.cc"
#include "QIE10.cc"
#include "SiPMModel.cc"
/*
12345678911234567892123456789312345678941234567895123456789612345678971234567898
*/
using std::cout;
using std::endl;

static double const tempConstraint = 0.1;

double correctSaturation(double charge, int NP, double /*prehit*/,
			 double xtalk = 0.) {
  if (charge >= NP) return /*1.e6*/7.7*NP;
  if (charge < NP*0.05) return charge;
  double val = TMath::Log(1.0 - charge/NP);
  val *= -NP;
  if ((xtalk > 0.) && (xtalk < 1.))
    val *= (1-xtalk);
  return val;
}

double respSlope(double val, int NP, double prehit = 0.) {
  if (val < 0) return 0.;
  assert(prehit >= 0);
  double slope = TMath::Exp(-val/NP)/NP;
  slope *= (1-prehit)*NP;
  return slope;
}

double analyticError(unsigned int NP, int pes, double val, double eff = 1.0) {
  double alpha = pes/double(NP);
  double sig2 = NP*TMath::Exp(-alpha)*(1-(1+alpha)*TMath::Exp(-alpha));
  return TMath::Sqrt(sig2 + val*(1-eff));
}

double timeDepResponse(SiPMModel& sipm, unsigned int pes, double dT = 0.) {
  static TF1 Y11Shape("Y11Shape", "exp(-0.0635-0.1518*x)*x**2.528", 0., 80.);

  int nbins = int(50./sipm.getTau()*5+0.5);
  double dt = 50./nbins;

  TH1D pedist("pedist", "pedist", nbins, 0., 50.);
  // for (int pe = 0; pe < pes; ++pe)
  //   pedist.Fill(Y11Shape.GetRandom());
  pedist.FillRandom("Y11Shape", pes);

  // pedist.Draw();
  // gPad->Update();
  // gPad->WaitPrimitive();

  double sum = 0.;
  for (int tbin = 1; tbin <= nbins; ++tbin) {
    sum += sipm.hitPixels((unsigned int)pedist.GetBinContent(tbin), 0., dT);
    sipm.expRecover(dt);
  }

  return sum;
}

void quickSim (Int_t NP = 15000, Int_t trials = 500, Double_t pctdamage = 0.,
	       double tempCoef = 0., double crossTalk = 0.,
	       Int_t QIEver = 0, Double_t gain = 1.0, 
	       Bool_t gainMatch = kFALSE, double timeConst = 10.,
	       bool doTimeDepSignal = false,
	       Color_t tcolor = kBlue, int method = 1) {
  Int_t const steps = 200;
  double cutOff = 15.;
  if (doTimeDepSignal)
    cutOff *= 5;
  TF1 * linear = new TF1("linear", "[0]*x", 1./(NP*1.1), cutOff);
  linear->SetLineWidth(2);
  linear->SetLineColor(kBlack);
  linear->SetLineStyle(kDotted);
  linear->SetParameter(0, 1.0);
  //Double_t x[steps], y[steps], ey[steps], ex[steps];
  TGraphErrors *Response = new TGraphErrors(steps);
  TGraphErrors *NormResponse = new TGraphErrors(steps);
  TGraph *NormRMS = new TGraph(steps);
  TGraphErrors * Corrected = new TGraphErrors(steps);
  Int_t pes = 0, step = 0, trial;
  Double_t val, sum, sum2, rms, corrval, corrsum, corrsum2, corrrms, wsum;
  SiPMModel SiPM(NP, timeConst);
  SiPM.setCrossTalk(crossTalk);
  SiPM.setTempDep(tempCoef);
  static TRandom3 rnd(0);
//   static TCanvas *c1 = 0;
//   static TCanvas *c2 = 0;
  Double_t lastStep = 0., lastpes = 0.;
  Double_t platCnt = 0;
//   for(step = 0; step < steps; step++) {
  char buf[200];
  TFile responsePlots("ToyResponsePlots.root", "update");
  TH1F::SetDefaultBufferSize(trials);
  //Int_t minStep = int(NP*10./100. + 0.5);
  Int_t nextStep;
  double expSteps = TMath::Log(NP*cutOff)/steps;
  int maxCharge = int(cutOff*NP*gain+0.5);
  if (QIEver == 8)
    maxCharge = 10000;
  else if (QIEver == 10)
    maxCharge = 300000;
  if (gainMatch) {
    gain = maxCharge/double(NP);
    std::cout << "gain: " << gain << " fC/pe\n";
  }
  do {
    //pes += 5*step;
    //pes += Int_t(step*sqrt(step));
    nextStep = int(pes*expSteps + 0.5); 
    pes += nextStep + 1;
    sum = 0.;
    sum2 = 0.;
    wsum = 0.;
    corrsum = 0.;
    corrsum2 = 0.;
    if (step%10 == 0) {
      cout << "Simulating " << pes << " PE's\n";
    }
    for(trial = 0; trial < trials; trial++) {
      double dT = rnd.Uniform(-tempConstraint,tempConstraint);

      if (method)
	val = SiPMModel::hitPixelsQ(NP, rnd.Poisson(pes), 1.-pctdamage,
				    tempCoef, dT, crossTalk);
      else {
	SiPM.resetSiPM();
	if (doTimeDepSignal)
	  val = timeDepResponse(SiPM, rnd.Poisson(pes), dT);
	else
	  val = SiPM.hitPixels(rnd.Poisson(pes), 0., dT);
      }
      val *= gain;
      if (QIEver == 8) {
	wsum += qie8.binWidth(val);
	val = qie8.findEdge(val);
      } else if (QIEver == 10) {
	wsum += qie10.binWidth(val);
	val = qie10.findEdge(val);
      }
      sum += val;
      sum2 += val*val;
      corrval = correctSaturation(val/gain, NP, pctdamage, crossTalk);

      corrsum += corrval;
      corrsum2 += corrval*corrval;
    }

    val = sum/Double_t(trials)/gain;
    rms = (sum2 - sum*sum/double(trials))/double(trials-1);
    if (rms < 0.) rms = 0.;
    rms = sqrt(rms)/gain;

    corrval = corrsum/trials;
    corrrms = (corrsum2 - corrsum*corrsum/double(trials))/double(trials-1);
    if (corrrms < 0) corrrms = 0.;
    corrrms = sqrt(corrrms);

    // rms = sqrt((sum2 - sum*sum/Double_t(trials))/
    // 	       Double_t(trials-1))/gain;

    // corrval = corrsum/trials;
    // corrrms = sqrt((corrsum2 - corrsum*corrsum/trials)/(trials-1));

    NormResponse->SetPoint(step, pes/Double_t(NP), 
			   val/double(NP));
    double slope = respSlope(pes, NP, pctdamage);
    double bwidth = wsum/double(trials)/gain;
    double error = analyticError(NP, pes, val, 1-pctdamage);
    double rawError = error;
    error = sqrt(error*error + bwidth*bwidth/12);

    // if (pes/double(NP) > 6.) {
    //   std::cout << "pes: " << pes << " NP: " << NP << '\n'
    // 		<< "slope: " << slope << " bin width: " << bwidth << '\n'
    // 		<< "val: " << val << " rms: " << rms << " analytic error: " 
    // 		<< rawError << " bin corrected error: " << error
    // 		<< '\n'
    // 		<< "corrected: " << corrval << " corr rms: " << corrrms 
    // 		<< " analytic error/slope: " << (error/slope) << '\n'
    // 		<< "pct error: " << (rms/slope/double(pes)) << '\n'
    // 		<< '\n';
    // }

    corrrms = sqrt(val + rawError*rawError + bwidth*bwidth/12);
    rms = sqrt(rms*rms + bwidth*bwidth/12);

    //if (rms < error) rms = error;
    // if (corrrms < error/slope) corrrms = error/slope;

    NormResponse->SetPointError(step, 0., rms/double(NP));
    NormRMS->SetPoint(step, pes/Double_t(NP),
		      corrrms/slope/double(pes));

    Corrected->SetPoint(step, pes/double(NP), 
			corrval/NP);
    Corrected->SetPointError(step, 0., corrrms/slope/double(NP));
    if (val <= lastStep) {
      platCnt++;
      val = lastStep;
    } else platCnt = 0;
    lastStep = val;
    lastpes = pes;
    Response->SetPoint(step, pes, val);
    Response->SetPointError(step, 0., rms);
  } while (/*(platCnt < 25) && */(++step < steps) && (pes < cutOff*NP));
  Response->Set(step);
  NormResponse->Set(step);
  NormRMS->Set(step);
  Corrected->Set(step);
  std::cout << "Nsteps: " << step << "\n";

  Response->SetMarkerColor(tcolor);
  Response->SetLineColor(tcolor);
  Response->SetFillColor(tcolor);

  sprintf(buf, "Response vs. PEs (%i pixels, %2.0f%% pre-hit, QIE%i)", 
	  NP, pctdamage*100, QIEver);
  Response->SetTitle(buf);
  sprintf(buf, "RespvPE%ipx%1.0fd", NP, pctdamage*100);
  Response->SetName(buf);
  Response->Write();

  NormResponse->SetMarkerStyle(kFullCircle);
  NormResponse->SetMarkerSize(0.2);
  NormResponse->SetMarkerColor(tcolor);
  NormResponse->SetLineColor(tcolor);
  NormResponse->SetLineWidth(2);
  NormResponse->SetFillColor(tcolor);
  NormResponse->SetFillStyle(3006);

  sprintf(buf, "Normalized Response vs. PEs "
	  "(%i pixels, %2.0f%% pre-hit, QIE%i)", 
	  NP, pctdamage*100, QIEver);
  NormResponse->SetTitle(buf);
  sprintf(buf, "NormRespvPE%ipx%1.0fd", NP, pctdamage*100);
  NormResponse->SetName(buf);
  NormResponse->Write();

  NormRMS->SetMarkerColor(tcolor+1);
  NormRMS->SetLineColor(tcolor+1);
  NormRMS->SetLineWidth(2);
  NormRMS->SetFillColor(tcolor+1);

  sprintf(buf, "Fractional error vs. PEs (%i pixels, %2.0f%% pre-hit, QIE%i)", 
	  NP, pctdamage*100, QIEver);
  NormRMS->SetTitle(buf);
  sprintf(buf, "ErrvPE%ipx%1.0fd", NP, pctdamage*100);
  NormRMS->SetName(buf);
  NormRMS->Write();

  // NormResponse->Draw("APLZ");
  // NormResponse->Draw("APLZ");

  // NormRMS->Draw("APL");

  Corrected->SetTitle(TString::Format("Corrected Response (%i pixels, %2.0f%% "
				      "pre-hit, QIE%i)", 
				      NP, pctdamage*100, QIEver));
  Corrected->SetName(TString::Format("CorrResp%ipx%1.0fd", 
				     NP, pctdamage*100));
  Corrected->SetMarkerColor(tcolor-1);
  Corrected->SetLineColor(tcolor-1);
  Corrected->SetFillColor(tcolor-1);
  Corrected->SetFillStyle(3018);
  Corrected->SetLineWidth(2);
  Corrected->Write();

  TCanvas * c1 = new TCanvas("c1", "c1", 600, 600);
  c1->Update();
  c1->SetLogx();
  c1->SetLogy();
  c1->SetGridx();
  c1->SetGridy();

  linear->Draw();
  NormResponse->Draw("LX");
  NormRMS->Draw("L");
  Corrected->Draw("3L");
  c1->Update();
  Corrected->SetMaximum(cutOff);
  c1->Modified();

  TCanvas * c2 = new TCanvas("c2", "correction", 600, 600);
  // c2->Update();
  c2->SetLogx();
  c2->SetLogy();
  c2->SetGridx();
  c2->SetGridy();
  linear->Draw();
  Corrected->Draw("3L");

  responsePlots.Close();
}
