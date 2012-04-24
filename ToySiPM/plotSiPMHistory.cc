#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TRegexp.h"
#include "TSystem.h"
#include "TStyle.h"

/*
  plotSiPMHistory parameters
  fname : input file name which is the output from SiPMHistoryStudy
  plotdir : directory to put the plots.  It needs to already exist
  intperx : pileup interactions per bunch crossing
  TestGev : 500 or 1000 that was used in SiPMHistoryStudy
  TestGev2 : 500 or 1000 that was used in SiPMHistorySudy
  pxfname : name of where to put the latex table data that summarized the 
            pixels that were still active on average
  respfname : name of file to put the latex table data that summarized the 
              response
  occfname : name of the file to put the latex table data that will summarize
             the average occupancy from pile up.
 */
void plotSiPMHistory(TString fname, 
		     TString plotdir, 
		     //int intperx, int TestGev, int TestGev2, 
		     TString pxfname = "", 
		     TString respfname = "" , TString occfname = "") {
  gStyle->SetOptStat(0);
  ifstream ifs(fname.Data(), ifstream::in);
  TString firstLine;
  firstLine.ReadLine(ifs);
  ifs.close();
  std::cout << firstLine << '\n';
  TTree Resps;
  if (Resps.ReadFile(fname.Data(), "TotPe/I:LivePxIntegral/F:NP/I:ped[3]/F:Layer/I:Response/F:GenPE/I:Response2/F:GenPE2/I:respForm[50]/F:genForm[50]/F") <= 0) return;
  if (Resps.GetEntries() <= 0) return;
  //Resps.Print();
  TRegexp fre("[0-9\\.]+");
  TRegexp ire("[0-9]+");
  int intperx = 1;
  TString tmps = firstLine(TRegexp("intperx [0-9]*"));
  // std::cout << tmps 
  // 	    << '\n';
  tmps = tmps(ire);
  std::cout << tmps 
	    << ' ';
  intperx = tmps.Atoi();
  int TestGev = 500;
  tmps = firstLine(TRegexp("Test1 [0-9]*"));
  tmps = tmps(ire, 5);
  std::cout << tmps 
	    << ' ';
  TestGev = tmps.Atoi();
  int TestGev2 = 1000;
  tmps = firstLine(TRegexp("Test2 [0-9]*"));
  tmps = tmps(ire, 5);
  std::cout << tmps 
	    << ' ';
  TestGev2 = tmps.Atoi();
  TString SiPMType = firstLine(TRegexp("[EO]DU"));
  tmps = firstLine(TRegexp("rPct1 [0-9\\.]*"));
  tmps = tmps(fre, 5);
  std::cout << tmps 
	    << ' ';
  double rPct1 = tmps.Atof();
  tmps = firstLine(TRegexp("rTime1 [0-9\\.]*"));
  tmps = tmps(fre, 6);
  std::cout << tmps 
	    << ' ';
  double rTime1 = tmps.Atof();
  tmps = firstLine(TRegexp("rTime2 [0-9\\.]*"));
  tmps = tmps(fre, 6);
  std::cout << tmps 
	    << ' ';
  double rTime2 = tmps.Atof();
  int NP = int(Resps.GetMaximum("NP"));
  tmps = firstLine(TRegexp("NP [0-9]*"));
  tmps = tmps(ire);
  std::cout << tmps 
	    << ' ';
  NP = tmps.Atoi();
  tmps = firstLine(TRegexp("etaBin [1-7]"));
  tmps = tmps(ire);
  std::cout << tmps 
	    << ' ';
  int etaBin = tmps.Atoi();
  tmps = firstLine(TRegexp("pescale [0-9\\.]*"));
  tmps = tmps(fre);
  std::cout << tmps 
	    << ' ';
  double pescale = tmps.Atof();

  TString devDir = TString::Format("NP_%i_%s_rTime1_%0.0f_rPct1_%0.2f_rTime2_%0.0f_pescale_%0.2f",
				   NP, SiPMType.Data(), rTime1, rPct1, rTime2, pescale);
  TString PUDir = TString::Format("etaBin_%i_PU_%i", etaBin, intperx);
  std::cout << '\n' 
	    << plotdir << '/' << devDir << '/' << PUDir
	    << '\n';
  gSystem->mkdir(plotdir + "/" + devDir + "/" + PUDir, true);
  TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
  std::cout << std::endl;
  int depths = int(Resps.GetMaximum("Layer"));
  char name[100];
  char title[200];
  char drawcmd[200];
  char cut[100];
  TH1F * pixels;
  TH1F * resp;
  TH1F * resp2;
  TH1F * occ;
  TH1F * peds;
  TString outfname;
  bool doPx = true;
  if (pxfname.Length() == 0) {
    pxfname = "/dev/null";
    // doPx = false;
  }
  bool doPed = doPx;
  bool doOcc = true;
  if (occfname.Length() == 0) {
    occfname = "/dev/null";
    // doOcc = false;
  }

  if (respfname.Length() == 0) {
    respfname = "/dev/null";
  }

  ofstream pxout(pxfname.Data(), ios::out | ios::app);
  ofstream respout(respfname.Data(), ios::out | ios::app);
  ofstream occout(occfname.Data(), ios::out | ios::app);
  pxout << intperx << " PU & tot";
  respout << intperx << " PU & " << TestGev << " GeV";
  stringstream pxline2(ios::in | ios::out);
  pxline2 << " & frac";
  stringstream respline2(ios::in | ios::out);
  respline2 << " & " << TestGev2 << " GeV";
//   stringstream respline2(ios::in | ios::out);
//   respline2 << " & frac";
  occout << intperx << " PU";
  // //first layer of next depth
  // int const maxDepth = 5;
  // int depthBounds[maxDepth] = {1,5,9,17,19};
  // if (etaBin > 4) {
  //   depthBounds[0] = 2;
  //   depthBounds[1] = 6;
  //   depthBounds[2] = 10;
  //   depthBounds[3] = 19;
  //   testOffset = 1;
  // }
  // double dpulseResp = 0., dpulseResp2 = 0., dtest = 0., dtest2 = 0.;

  for (int d=0; d<depths; ++d) {
    //pedestal
    if (doPed) {
      sprintf(name, "PedsDepth%i", d);
      sprintf(title, "pedestal: %s depth %i,%i pixel/mm^{2},%i PU,"
	      "%0.0f rTime1,%0.2f rPct1,%0.0f rTime2",
	      SiPMType.Data(), d, NP, intperx, rTime1, rPct1, rTime2);
      //double maxPed = Resps.GetMaximum("ped0");
      sprintf(drawcmd, "ped>>Peds");
      sprintf(cut, "Layer==%i", d);
      //std::cout << drawcmd << ", " << cut << "\n";
      Resps.Draw(drawcmd, cut);
      gDirectory->GetObject("Peds", peds);
      //std::cout << name << "\n";
      if (peds) {
	//std::cout << "got peds" << "\n";
	if (peds->GetEntries() < 1.) continue;
	peds->SetTitle(title);
	peds->Draw();
	c1->Update();

	outfname = plotdir + '/' + devDir + '/' + PUDir;
	outfname += "/";
	outfname += name;
	outfname += ".eps";
	c1->Print(outfname.Data());
      }
      delete peds;
    }

    //std::cout << pxdensity << " " << NP << std::endl;
    if (doPx) {
      sprintf(name, "LivePixelsDepth%i", d);
      sprintf(title, "Live pixels: %s depth %i,%i pixel/mm^{2},%i PU,"
	      "%0.0f rTime1,%0.2f rPct1,%0.0f rTime2",
	      SiPMType.Data(), d, NP, intperx, rTime1, rPct1, rTime2);
      sprintf(drawcmd, "LivePxIntegral>>LivePixels(20)");
      sprintf(cut, "Layer==%i", d);
      Resps.Draw(drawcmd, cut);
      sprintf(drawcmd, "LivePxIntegral/NP>>AvgPx(20)");
      Resps.Draw(drawcmd, cut);
      //gDirectory->ls();
      gDirectory->GetObject("LivePixels", pixels);
      TH1F * avgpx;
      gDirectory->GetObject("AvgPx", avgpx);
      //std::cout << name << "\n";
      if (pixels) {
	//std::cout << "got pixels" << "\n";
	if (pixels->GetEntries() < 1.) continue;
	if (avgpx->GetEntries() < 1.) continue;
	pixels->SetTitle(title);
	pixels->Draw();
	c1->Update();

	pxout << " & "  << std::setprecision(0) << std::fixed << pixels->GetMean();
	pxout.unsetf(ios_base::floatfield);
	pxline2  << " & " << std::setprecision(3) << avgpx->GetMean();

	outfname = plotdir + '/' + devDir + '/' + PUDir;
	outfname += '/';
	outfname += name;
	outfname += ".eps";
	c1->Print(outfname.Data());
      }
      delete pixels;
      delete avgpx;
    }

    sprintf(name, "RespDepth%iTest%iGeV", d, TestGev);
    sprintf(title, "Response %i GeV: %s depth %i,%i pixel/mm^{2},%i PU,"
	    "%0.0f rTime1,%0.2f rPct1,%0.0f rTime2",
	    TestGev, SiPMType.Data(), d, NP, intperx, rTime1, rPct1, rTime2);
    sprintf(drawcmd, "Response/GenPE>>RespHist(20,0.,1.01)");
    sprintf(cut, "(Layer==%i)&&(GenPE>0)", d);
    //std::cout << drawcmd << ", " << cut << "\n";
    Resps.Draw(drawcmd, cut);
    //gDirectory->ls();
    gDirectory->GetObject("RespHist", resp);
    //std::cout << name << "\n";
    if (resp) {
      //std::cout << "got pixels" << "\n";
      if (resp->GetEntries() < 1.) continue;
      resp->SetTitle(title);
      resp->Draw();
      c1->Update();

      respout << " & " << std::setprecision(4) << resp->GetMean();

      outfname = plotdir + '/' + devDir + '/' + PUDir;
      outfname += '/';
      outfname += name;
      //outfname.ReplaceAll('.', "-");
      outfname += ".eps";
      c1->Print(outfname.Data());
    }
    delete resp;

    if (TestGev2 > 0) {
      sprintf(name, "RespDepth%iTest%iGeV", d, TestGev2);
      sprintf(title, "Response %i GeV: %s depth %i,%i pixel device,%i PU,"
	      "%0.0f rTime1,%0.2f rPct1,%0.0f rTime2",
	      TestGev2, SiPMType.Data(), d, NP, intperx, rTime1, rPct1, rTime2);
      sprintf(drawcmd, "Response2/GenPE2>>RespHist2(20,0.,1.01)");
      sprintf(cut, "(Layer==%i)&&(GenPE2>0)", d);
      //std::cout << drawcmd << ", " << cut << "\n";
      Resps.Draw(drawcmd, cut);
      //gDirectory->ls();
      gDirectory->GetObject("RespHist2", resp2);
      //std::cout << name << "\n";
      if (resp2) {
	//std::cout << "got pixels" << "\n";
	if (resp2->GetEntries() < 1.) continue;
	resp2->SetTitle(title);
	resp2->Draw();
	c1->Update();

	respline2 << " & " << std::setprecision(4) << resp2->GetMean();

	outfname = plotdir + '/' + devDir + '/' + PUDir;
	outfname += '/';
	outfname += name;
	//outfname.ReplaceAll('.', "-");
	outfname += ".eps";
	c1->Print(outfname.Data());
      }
      delete resp2;
    }

    if (doOcc) {
      sprintf(name, "OccDepth%i", d);
      sprintf(title, "Occupancy: %s depth %i,%i pixel device,%i PU,"
	      "%0.0f rTime1,%0.2f rPct1,%0.0f rTime2",
	      SiPMType.Data(), d, NP, intperx, rTime1, rPct1, rTime2);
      sprintf(drawcmd, "TotPe>>peHist(10)");
      sprintf(cut, "Layer==%i", d);
      //std::cout << drawcmd << ", " << cut << "\n";
      Resps.Draw(drawcmd, cut);
      //gDirectory->ls();
      gDirectory->GetObject("peHist", occ);
      //std::cout << name << "\n";
      if (occ) {
	//std::cout << "got pixels" << "\n";
	occ->SetTitle(title);
	occ->Draw();
	c1->Update();
	
	occout << " & " << int(occ->GetMean()+0.5);
	
	outfname = plotdir + '/' + devDir + '/' + PUDir;
	outfname += '/';
	outfname += name;
	//outfname.ReplaceAll('.', "-");
	outfname += ".eps";
	c1->Print(outfname.Data());
      }
      delete occ;
    }
  }

  pxout << " \\\\ \n"
	<< pxline2.str() << "\\\\ \n\\hline \n";
  respout << " \\\\ \n";
  if (TestGev2 > 0)
    respout << respline2.str() << "\\\\ \n";
  respout << "\\hline \n";
  occout << " \\\\ \n\\hline \n";

  occout.close();
  respout.close();
  pxout.close();
  delete c1;
}
