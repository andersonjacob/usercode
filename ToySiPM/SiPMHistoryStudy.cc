#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include "TList.h"
#include "TCollection.h"
#include "TString.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TRegexp.h"
//#include "QIEBinEdges.h"

#include "QIE10.cc"
#include "SiPMModel.cc"

//static TRandom3 rnd;
//static int GenMode;
ofstream outfs;
static int bxspacing = 25;
static double Y11Times[2] = {15, 35};
static int Y11Time = int(Y11Times[0] + Y11Times[1] + 0.5);
static int dt = 1;
static int freq = 100;

double TriangleSmear(double * x, double * par) {
  double h = 2/(par[0] + par[1]);
  if (x[0] <= par[0])
    return h/par[0]*x[0];
  double val = h - h/par[1]*(x[0]-par[0]);
  return (val > 0) ? val : 0.0;
}


double SiPMHistory (int pes, double t, SiPMModel& SiPM) {
  return SiPM.hitPixels(pes, SiPM.recoveryModel(t));
}

double applyHistory (SiPMModel& SiPM, int past, TList& hists, int intperx,
		     float pescale, int& totalpes, double * last3bx) {
  //static TRandom3 rnd;
  double pes;
  totalpes = 0;
  TH1D * hist;
  double bxresp;
  // double last3bx[3];
  //std::map<int, double> history;
  int pesonsipm;
  for (int t = past; t >= 0; t-=bxspacing) {
    TIter nextHist(&hists);
    bxresp = 0.;
    int tmpInts = gRandom->Poisson(intperx);
    pesonsipm = 0;
    while ((hist = (TH1D *)nextHist())) {
      //double norm = hist->Integral();
      for (int scatter=0; scatter < tmpInts; ++scatter) {
	pes = hist->GetRandom()*pescale;
	pes = gRandom->PoissonD(pes);
	pesonsipm += int(pes+0.5);
      }
    }
    totalpes += pesonsipm;
    if (pesonsipm > 0) {
      bxresp = SiPMHistory (pesonsipm, t, SiPM);
    }

    if (t<3*bxspacing) {
      //std::cout << "t: " << t << '\n';
      if (t>bxspacing) last3bx[2] = bxresp;
      if (t>0) last3bx[1] = bxresp;
      last3bx[0] = bxresp;
    }
    //outfs.flush();
  }
  double totalcharge = SiPM.totalCharge();
  // outfs << totalpes << ' ';
  // outfs << totalcharge << ' '
  // 	<< NP << ' ';
  // for (int i=2; i>=0; --i) outfs << last3bx[i] << ' ';
  return totalcharge;
}

double applyTestPulse(SiPMModel& SiPM, int testp, double peperns, 
		      int * genpebins = 0, 
		      float * resppebins = 0) {
  if (testp < 1) return double(testp);
  static TF1 smearf("smearf", TriangleSmear, 0, Y11Time, 2);
  static bool init = true;
  if (init) {
    smearf.SetParameters(Y11Times);
    init = false;
  }
  double pulseResp = 0;
  TH1I smearHits("smearHits", "smearHits", Y11Time/dt, 0, Y11Time);
  smearHits.FillRandom("smearf", testp);
  for (int bin = 1; bin <= smearHits.GetNbinsX(); ++bin) {
    SiPM.hitPixels(gRandom->Poisson(peperns), 0.0);
    double binResp = 
      SiPM.hitPixels(int(smearHits.GetBinContent(bin)), 0.0);
    if ((genpebins) && (resppebins)) {
      genpebins[bin-1] = int(smearHits.GetBinContent(bin));
      resppebins[bin-1] = binResp;
    }
    pulseResp += binResp;
    SiPM.expRecover(dt);
  }
  return pulseResp;
}

void EDUHistoryStudy (int etaBin = 3, int intperx = 1, float TestPulse[19] = 0,
		      int NP = 15000, int past = 1000000, 
		      double rTime1 = 160, double recoveryPct1 = 0.47,
		      double rTime2 = 1E6, float pescale = 1.0, 
		      TString histfilename = "peSpectrumByLayer7EtaBins.root",
		      float Test2[19] = 0) {
  //GenMode = gm;
  SiPMModel layerSiPM(NP,rTime1,recoveryPct1,rTime2);
  // double * layerSiPM = new double[NP];
  int testOffset = 0;
  //first layer of next depth
  int const depths = 5;
  int depthBounds[depths] = {1,5,11,17,19};
  if (etaBin > 4) {
    depthBounds[0] = 2;
    depthBounds[1] = 6;
    depthBounds[2] = 10;
    depthBounds[3] = 19;
    testOffset = 1;
  }
  int depth = 0;
  TFile histFile(histfilename.Data());
  TH1D * layerHist = 0;
  char name[50];
  double pulseResp;
  //double damage;
  TList hists;
  double pulseResp2;
  int TP, TP2;
  double dpulseResp = 0., dpulseResp2 = 0., dtest = 0., dtest2 = 0., 
    dtotalcharge = 0, dlast3bx[3] = {0.,0.,0.}, last3bx[3];
  int dtotalpes = 0, dNP = 0;
  float * respForm = new float[Y11Time/dt];
  int * genForm = new int[Y11Time/dt];
  for (int l = 0; l < 19; ++l) {
    if (l >= depthBounds[depth]) {
      outfs << dtotalpes << ' ' << dtotalcharge << ' ' << dNP << ' ';
      for (int bx = 0; bx < 3; ++bx) {
	outfs << dlast3bx[bx] << ' ';
	dlast3bx[bx] = 0.;
      }
      outfs << depth << ' '
	    << dpulseResp << ' ' << dtest << ' '
	    << dpulseResp2 << ' ' << dtest2 << ' ';
      for (int tstep = 0; tstep < Y11Time/dt; ++tstep)
	outfs << respForm[tstep] << ' ';
      for (int tstep = 0; tstep < Y11Time/dt; ++tstep)
	outfs << genForm[tstep] << ' ';
      outfs << '\n';
      // std::cout << "Depth: " << depth
      // 		<< " Test1: " << dtest << " Response: " << dpulseResp
      // 		<< " Test2: " << dtest2 << " Response: " << dpulseResp2
      // 		<< std::endl;
      dpulseResp = 0.;
      dpulseResp2 = 0.;
      dtest = 0;
      dtest2 = 0;
      dtotalpes = 0;
      dtotalcharge = 0.;
      dNP = 0;
      depth++;
    }
    layerSiPM.resetSiPM();
    // resetSiPM(NP, layerSiPM);
    sprintf(name, "peSpectrum_etaBin%i_layer%i", etaBin, l);
    layerHist = 0;
    histFile.GetObject(name, layerHist);
    if (!layerHist) continue;
    if (layerHist->GetRMS() <= 0.) continue;
    hists.Add(layerHist);
    // outfs << "applying History with GenMode " << GenMode << std::endl;
    int totalpes = 0;
    double tot = applyHistory(layerSiPM, past, hists, intperx, pescale, 
			      totalpes, last3bx);
    int testl = l;
    double peperns = double(totalpes)/past;
    if (testl >= testOffset) testl -= testOffset;
    TP2 = 0;
    pulseResp2 = 0.;
    if (Test2) {
      TP2 = gRandom->Poisson(int(Test2[testl]+0.5));
      SiPMModel SiPM2(layerSiPM);
      pulseResp2 = applyTestPulse(SiPM2, TP2, peperns);
      // dtest2 += Test2[testl];
      dtest2 += TP2;
      dpulseResp2 += pulseResp2;
    }
    TP = gRandom->Poisson(int(TestPulse[testl]+0.5));
    for (int tstep = 0; tstep < Y11Time/dt; ++tstep) {
      genForm[tstep] = 0;
      respForm[tstep] = 0.;
    }
    pulseResp = applyTestPulse(layerSiPM, TP, peperns, genForm, respForm);
    // std::cout << "layer: " << testl << " Test: " << TestPulse[testl]
    // 	      << " pes: " << TP << " response: " << pulseResp
    // 	      << " peperns: " << peperns
    // 	      << std::endl;
    dpulseResp += pulseResp;
    // dtest += TestPulse[testl];
    dtest += TP;
    dtotalpes += totalpes;
    dtotalcharge += tot;
    dNP += NP;
    for (int bx = 0; bx < 3; ++bx) dlast3bx[bx] += last3bx[bx];
    hists.Clear();
  }
  outfs << dtotalpes << ' ' << dtotalcharge << ' ' << dNP << ' ';
  for (int bx = 0; bx < 3; ++bx) {
    outfs << dlast3bx[bx] << ' ';
    dlast3bx[bx] = 0.;
  }
  outfs << depth << ' '
	<< dpulseResp << ' ' << dtest << ' '
	<< dpulseResp2 << ' ' << dtest2 << ' ';
  for (int tstep = 0; tstep < Y11Time/dt; ++tstep)
    outfs << respForm[tstep] << ' ';
  for (int tstep = 0; tstep < Y11Time/dt; ++tstep)
    outfs << genForm[tstep] << ' ';
  outfs << '\n';
  outfs.flush();
  delete[] respForm;
  delete[] genForm;
  histFile.Close();
}

void ODUHistoryStudy (int etaBin = 3, int intperx = 1, float TestPulse[19] = 0,
		      int NP = 400, int past = 1000, 
		      double rTime1 = 200, double recoveryPct1 = 0.91,
		      double rTime2 = 1000, float pescale = 1.0, 
		      TString histfilename = "peSpectrumByLayer7EtaBins.root",
		      float Test2[19] = 0) {
  // double * depthSiPM = 0;//new double[NP];
  SiPMModel depthSiPM(NP, rTime1, recoveryPct1, rTime2);
  int testOffset = 0;

  //first layer of next depth
  int depthBounds[5] = {1,5,11,17,19};
  if (etaBin > 4) {
    depthBounds[0] = 2;
    depthBounds[1] = 6;
    depthBounds[2] = 12;
    depthBounds[3] = 19;
    testOffset = 1;
  }

  TFile histFile(histfilename.Data());
  TH1D * layerHist = 0;
  char name[50];
  double pulseResp, pulseResp2;
  int depth = 0;
  TList hists;
  hists.Clear();
  int testp = 0, testp2 = 0;
  double truetp = 0., truetp2 = 0.;

  //set number of pixels for the depths
  int LargeNP = int(NP*4.5+0.5);
  int HONP = int(7.*NP+0.5);
  depthSiPM.setNP(LargeNP);
  int currNP = LargeNP;
  double tot;
  double last3bx[3];
  double peperns;
  int totalpes;
  float * respForm = new float[Y11Time/dt];
  int * genForm = new int[Y11Time/dt];
  for (int l = 0; l < 19; ++l) {
    if (l >= depthBounds[depth]) {
      tot = applyHistory(depthSiPM, past, hists, intperx, pescale, 
			 totalpes, last3bx);
      peperns = double(totalpes)/past;
      pulseResp2 = 0.;
      if (Test2) {
	SiPMModel SiPM2(depthSiPM);
	pulseResp2 = applyTestPulse(SiPM2, testp2, peperns);
      }
      for (int tstep = 0; tstep < Y11Time/dt; ++tstep) {
	genForm[tstep] = 0;
	respForm[tstep] = 0.;
      }
      pulseResp = applyTestPulse(depthSiPM, testp, peperns,
				 genForm, respForm);

      outfs << totalpes << ' ' << tot << ' ' << currNP << ' ';
      for (int bx = 0; bx < 3; ++bx) outfs << last3bx[bx] << ' ';
      outfs << depth << ' '
	    << pulseResp << ' ' << testp << ' '
	    << pulseResp2 << ' ' << testp2 << ' ';
      for (int tstep = 0; tstep < Y11Time/dt; ++tstep)
	outfs << respForm[tstep] << ' ';
      for (int tstep = 0; tstep < Y11Time/dt; ++tstep)
	outfs << genForm[tstep] << ' ';
      outfs << '\n';
      ++depth;
      if (depth==4) {
	depthSiPM.setNP(HONP);
	currNP = HONP;
      } else {
	depthSiPM.setNP(LargeNP);
	currNP = LargeNP;
      }
      depthSiPM.resetSiPM();
      hists.Clear();
      testp = 0;
      truetp = 0.;
      testp2 = 0;
      truetp2 = 0.;
    }
    sprintf(name, "peSpectrum_etaBin%i_layer%i", etaBin, l);
    layerHist = 0;
    histFile.GetObject(name, layerHist);
    if (!layerHist) continue;
    if (layerHist->GetRMS() <= 0.) continue;
    hists.Add(layerHist);
    int testl = l;
    if (testl>=testOffset) testl -= testOffset;
    testp += gRandom->Poisson(int(TestPulse[testl]+0.5));
    truetp += TestPulse[testl];
    if (Test2) {
      testp2 += gRandom->Poisson(int(Test2[testl]+0.5));
      truetp2 += Test2[testl];
    }
  }
  tot = applyHistory(depthSiPM, past, hists, intperx, pescale, 
		     totalpes, last3bx);
  peperns = double(totalpes)/past;
  pulseResp2 = 0.;
  if (Test2) {
    SiPMModel SiPM2(depthSiPM);
    pulseResp2 = applyTestPulse(SiPM2, testp2, peperns);
  }
  for (int tstep = 0; tstep < Y11Time/dt; ++tstep) {
    genForm[tstep] = 0;
    respForm[tstep] = 0.;
  }
  pulseResp = applyTestPulse(depthSiPM, testp, peperns,
			     genForm, respForm);

  outfs << totalpes << ' ' << tot << ' ' << ' ' << currNP << ' ';
  for (int bx = 0; bx < 3; ++bx) outfs << last3bx[bx] << ' ';
  outfs << depth << ' '
	<< pulseResp << ' ' << truetp << ' '
	<< pulseResp2 << ' ' << truetp2 << ' ';
  for (int tstep = 0; tstep < Y11Time/dt; ++tstep)
    outfs << respForm[tstep] << ' ';
  for (int tstep = 0; tstep < Y11Time/dt; ++tstep)
    outfs << genForm[tstep] << ' ';
  outfs << '\n';
  outfs.flush();
  hists.Clear();
  histFile.Close();
}

/*
  SiPMHistoryStudy parameter description:
  trials : number of simulated events to do
  etaBin : 1-7 eta bin to use.  4 ietas
  intperx : average number of interactions per bunch crossing (pile up)
  NP : pixels/mm^2 of the SiPM
  SiPMtype : ODU or EDU, 1 is EDU anything else is ODU
  past : number of ns in the past to run the simulation.  
         Should probably be the same as rTime2.
  rTime1 : time in ns of the elbow in the SiPM recovery
  recoveryPct1 : fraction of the recovery at rTime1, from 0 to 1
  rTime2 : time at which SiPM fully recovers
  out_fname : file nane to put the results.
  pescale : use 1.  adjust the pe/gev which is 30
  histfilename : file that has the minBias information for the pileup simulation
                 don't change unless you have a newer file
  test1 : could be 500 or 1000, I have two files
          LayerPes500GeV.txt and LayerPes1000GeV.txt
  secondTest : could be 500 or 1000, this will pick up the second test file to 
               do it at the same time.
  seed : random seed being used for the simulation.  I suggest 0.
 */
void SiPMHistoryStudy (int trials, int etaBin, int intperx, int NP, 
		       int SiPMtype, int past,  
		       double rTime1, double recoveryPct1,
		       double rTime2, 
		       TString out_fname = "ResponseOutput.asc",
		       float pescale = 1., //int gm = 0,
		       TString histfilename = "peSpectrumByLayer7EtaBins.root",
		       int test1 = 500, 
		       int secondTest = 0,
		       unsigned int seed = 0) {
  //static TRandom3 rnd;
  // GenMode = gm;
  if (seed > 0) {
    gRandom->SetSeed(seed);
    SiPMModel::setSeed(seed);
  }

  outfs.open(out_fname.Data(), ios::trunc);
  outfs << "# " << " etaBin " << etaBin 
	<< " intperx " << intperx
	<< " NP " << NP 
	<< " SiPMtype " << ((SiPMtype == 1) ? "EDU" : "ODU")
	<< " rTime1 " << TString::Format("%0f", rTime1)
	<< " rPct1 " << TString::Format("%0.2f", recoveryPct1)
	<< " rTime2 " << TString::Format("%0f", rTime2)
	<< " pescale " << TString::Format("%0f", pescale)
	<< " Test1 " << test1 
	<< " Test2 " << secondTest
	<< '\n';
  // read files with the test pulses.
  TTree pes;
  TString TestFileName = TString::Format("LayerPes%iGeV.txt", test1);
  pes.ReadFile(TestFileName.Data(), 
      "L0:L1:L2:L3:L4:L5:L6:L7:L8:L9:L10:L11:L12:L13:L14:L15:L16:L17:L18");
  float TestP[19];
  char bname[4];
  //pes.Print();
  for (int i=0; i<19; ++i) {
    sprintf(bname, "L%i", i);
    pes.SetBranchAddress(bname, &(TestP[i]));
  }
  TTree pes2;
  TString test2 = TString::Format("LayerPes%iGeV.txt", secondTest);
  float TestP2[19];
  if (secondTest > 0) {
    pes2.ReadFile(test2.Data(), 
	 "L0:L1:L2:L3:L4:L5:L6:L7:L8:L9:L10:L11:L12:L13:L14:L15:L16:L17:L18");
    //pes2.Print();
    for (int i=0; i<19; ++i) {
      sprintf(bname, "L%i", i);
      pes2.SetBranchAddress(bname, &(TestP2[i]));
    }
  }
  for (int trial = 0; trial<trials; ++trial) {
    // pick a random test pulse.
    pes.GetEntry(trial);
    if (secondTest > 0)
      pes2.GetEntry(trial);
    if (trial%freq == 0)
      std::cout << "starting trial " << trial+1 << " of " << trials << '\n';
    //std::cout << " pescale: " << pescale;
    // scale back test for lower pde
    for (int l=0; l<19; l++) {
      TestP[l] *= pescale;
      TestP2[l] *= pescale;
//       std::cout << " L" << l << " pes: " << TestP[l];
    }

    // do one test
    if (SiPMtype == 1)
      EDUHistoryStudy(etaBin, intperx, TestP, NP, past, rTime1, 
		      recoveryPct1, rTime2, pescale, histfilename,
		      ((secondTest>0)? TestP2 : 0));
    else
      ODUHistoryStudy(etaBin, intperx, TestP, NP, past, rTime1,
		      recoveryPct1, rTime2, pescale, histfilename,
		      ((secondTest>0)? TestP2 : 0));
  }
  outfs.close();
}
