#include "TROOT.h"
#include "TTree.h"
#include "TCut.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include "TMath.h"

#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooArgSet.h"
#include "RooPlot.h"
#include "RooArgList.h"
#include "RooAddPdf.h"
#include "RooExponential.h"
#include "RooBifurGauss.h"
#include "RooCBShape.h"
#include "RooHist.h"
#include "RooFitResult.h"
#include "RooChebyshevPDF.h"
#include "RooWorkspace.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"


char const * jetCut = "abs(JetPFCor_Eta[0]-JetPFCor_Eta[1])<2.5 && "
  "sqrt(JetPFCor_Pt[0]**2+JetPFCor_Pt[1]**2+2*JetPFCor_Pt[0]*JetPFCor_Pt[1]*cos(JetPFCor_Phi[0]-JetPFCor_Phi[1]))>40. && "
"abs(abs(abs(JetPFCor_Phi[0]-event_met_pfmetPhi)-TMath::Pi())-TMath::Pi())>0.5 && "
  "event_met_pfmet>25. && W_mt>40. && JetPFCor_Pt[0]>30."
  " && JetPFCor_Pt[1]>30. && JetPFCor_Pt[2]<30."
  " && JetPFCor_bDiscriminator[0]<1.19"
  " && JetPFCor_bDiscriminator[1]<1.19";

//   " && abs(abs(abs(JetPFCor_Phi[0]-event_met_pfmetPhi)-TMath::Pi())-TMath::Pi()) >0.5"
//    " && abs(abs(abs(W_phi-atan2((JetPFCor_Py[0]+JetPFCor_Py[1]),(JetPFCor_Px[0]+JetPFCor_Px[1])))-TMath::Pi())-TMath::Pi())>2.74";
      
// "sqrt(JetPFCor_Pt[0]**2+JetPFCor_Pt[1]**2+2*cos(abs(abs(abs(JetPFCor_Phi[0]-JetPFCor_Phi[1])-TMath::Pi())-TMath::Pi())) ) >50 && "

//   "sqrt(JetPFCor_Pt[0]**2+JetPFCor_Pt[1]**2+2*JetPFCor_Pt[0]*JetPFCor_Pt[1]*cos(abs(abs(abs(JetPFCor_Phi[0]-JetPFCor_Phi[1])-TMath::Pi())-TMath::Pi()))) >30. && "
//"abs(JetPFCor_Eta[0]-JetPFCor_Eta[1])<1.5 &&"

char const * muCut = "W_muon_pt>25."
   " && (W_muon_trackiso+W_muon_hcaliso+W_muon_ecaliso-event_RhoForLeptonIsolation*3.141592653589*0.09)/W_muon_pt<0.1"
   " && abs(W_muon_d0bsp)<0.02";

char const * eleCut = "W_electron_et>30."
   " && (W_electron_trackiso+W_electron_hcaliso+W_electron_ecaliso-event_RhoForLeptonIsolation*3.141592653589*0.09)/W_electron_pt<0.05"
   " && W_electron_isWP80"
   " && ((abs(W_electron_eta)<1.5 &&"
   " abs(W_electron_deltaphi_in)<0.03"
   " && abs(W_electron_deltaeta_in)<0.004)"
   " || (abs(W_electron_eta)>1.5"
   " && abs(W_electron_deltaphi_in)<0.02"
   " && abs(W_electron_deltaeta_in)<0.005))"
   " && sqrt((W_electron_vx-event_BeamSpot_x)**2+(W_electron_vy-event_BeamSpot_y)**2)<0.02";

/*
muCut = "W_muon_pt>25."
   " && (W_muon_trackiso+W_muon_hcaliso+W_muon_ecaliso-event_RhoForLeptonIsolation*3.141592653589*0.09)/W_muon_pt<0.1"
   " && abs(abs(abs(JetPFCor_Phi[0]-(W_muon_phi)-TMath::Pi())-TMath::Pi()))>2.1"
   " && abs(abs(abs(JetPFCor_Phi[1]-(W_muon_phi)-TMath::Pi())-TMath::Pi()))>2.3"
   " && abs(W_muon_d0bsp)<0.02";

eleCut = "W_electron_et>30."
   " && (W_electron_trackiso+W_electron_hcaliso+W_electron_ecaliso-event_RhoForLeptonIsolation*3.141592653589*0.09)/W_electron_pt<0.05"
   " && W_electron_isWP80 && ((abs(W_electron_eta)<1.5"
   " && abs(W_electron_deltaphi_in)<0.03"
   " && abs(W_electron_deltaeta_in)<0.004)"
   " || (abs(W_electron_eta)>1.5 && abs(W_electron_deltaphi_in)<0.02"
   " && abs(W_electron_deltaeta_in)<0.005))"
   " && abs(abs(abs(JetPFCor_Phi[0]-(W_electron_phi)-TMath::Pi())-TMath::Pi()))>1.8"
   " && abs(abs(abs(JetPFCor_Phi[1]-(W_electron_phi)-TMath::Pi())-TMath::Pi()))>2."
   " && sqrt((W_electron_vx-event_BeamSpot_x)**2+(W_electron_vy-event_BeamSpot_y)**2)<0.02";
*/

void RooWjjSidebandFitter(bool fitAll = false, bool useHists = true) {

  TString const sideband1("(Mass2j_PFCor<65.)&&(Mass2j_PFCor>30.)");
  TString const sideband2("(Mass2j_PFCor>95.)&&(Mass2j_PFCor<125.)");
  TCut massCut("(Mass2j_PFCor>30.)&&(Mass2j_PFCor<125.)");

  TCut* jetsSel;
  TCut* mu_selection;
  TCut* ele_selection;

  jetsSel = new TCut (jetCut);
  mu_selection = new TCut ( (*jetsSel) && TCut(muCut) );
  ele_selection = new TCut ( (*jetsSel) && TCut(eleCut) );
  
  *mu_selection = TCut( *mu_selection && massCut );
  *ele_selection = TCut( *ele_selection && massCut );

  TTree * treeTemp;


  // TFile wjetsShape_mu_file("/uscms_data/d2/kalanand/WjjTrees/WmunuJets_CMSSW415-Spring11MC_WJets.root", "READ");
  // TFile wjetsShape_mu_file("/uscms_data/d2/kalanand/WjjTrees/WmunuJets_CMSSW423-Summer11MC_WJets.root", "READ");
  TFile wjetsShape_mu_file("/uscms_data/d2/kalanand/WjjTrees/ReducedTree/RD_WmunuJets_CMSSW423-Summer11MC_WJets.root", "READ");
  treeTemp = (TTree*) wjetsShape_mu_file.Get("WJet");
  //ActivateTreeBranches(*treeTemp);
  gROOT->cd();
  TTree* tree1 = treeTemp->CopyTree(*mu_selection);


  // TFile wjetsShape_ele_file("/uscms_data/d2/kalanand/WjjTrees/WenuJets_CMSSW415-Spring11MC_WJets.root", "READ");
  // TFile wjetsShape_ele_file("/uscms_data/d2/kalanand/WjjTrees/WenuJets_CMSSW423-Summer11MC_WJets.root", "READ");
  TFile wjetsShape_ele_file("/uscms_data/d2/kalanand/WjjTrees/ReducedTree/RD_WenuJets_CMSSW423-Summer11MC_WJets.root", "READ");
  treeTemp = (TTree*) wjetsShape_ele_file.Get("WJet");
  //ActivateTreeBranches(*treeTemp, true);
  gROOT->cd();
  TTree* tree2 = treeTemp->CopyTree(*ele_selection);
  RooRealVar mVjj("fit_mlvjj", "m_{l#nu jj}", 160., 700.0, "GeV/c^{2}");
  mVjj.setBins(int((mVjj.getMax()-mVjj.getMin())/20.+0.5));
  RooRealVar mjj("Mass2j_PFCor", "m_{jj}", 20., 600., "GeV/c^{2}");
  mjj.setBins(int((mjj.getMax()-mjj.getMin())/10.+0.5));
  RooCategory channel("channel", "channel");
  channel.defineType("muon");
  channel.defineType("electron");

  RooDataSet * edata = new RooDataSet("edata", "edata", RooArgSet(mVjj,mjj),
				      RooFit::Import(*tree2));
  RooDataSet * mdata = new RooDataSet("mdata", "mdata", RooArgSet(mVjj,mjj),
				      RooFit::Import(*tree1));
  wjetsShape_mu_file.Close();
  wjetsShape_ele_file.Close();
  delete tree1;
  delete tree2;

  edata->Print();
  mdata->Print();
  RooDataSet wjetsData("wjetsData", "wjetsData", RooArgSet(mVjj,mjj),
		       RooFit::Index(channel),
		       RooFit::Import("muon", *mdata), 
		       RooFit::Import("electron", *edata));
  delete edata;
  delete mdata;
  wjetsData.Print("v");

  //Build pdf
  RooRealVar sigmaL1("sigmaL1", "#sigma_{L,1}", 22., 0., 1000., "GeV/c^{2}");
  sigmaL1.setPlotLabel(sigmaL1.GetTitle());
  RooRealVar sigmaR1("sigmaR1", "#sigma_{R,1}", 90., 0., 1000., "GeV/c^{2}");
  sigmaR1.setPlotLabel(sigmaR1.GetTitle());
  RooRealVar mean1("mean1", "<m>_{1}", 350., 0., 1000., "GeV/c^{2}");
  mean1.setPlotLabel(mean1.GetTitle());

  RooRealVar sigmaL2("sigmaL2", "#sigma_{L,2}", 22., 0., 2000., "GeV/c^{2}");
  sigmaL2.setPlotLabel(sigmaL2.GetTitle());
  RooRealVar sigmaR2("sigmaR2", "#sigma_{R,2}", 150., 0., 2000., "GeV/c^{2}");
  sigmaR2.setPlotLabel(sigmaR2.GetTitle());
  RooRealVar mean2("mean2", "<m>_{2}", 375, 0., 2000., "GeV/c^{2}");
  mean2.setPlotLabel(mean2.GetTitle());

  RooRealVar sigmaL3("sigmaL3", "#sigma_{L,3}", 20., 0., 3000., "GeV/c^{2}");
  sigmaL3.setPlotLabel(sigmaL3.GetTitle());
  RooRealVar sigmaR3("sigmaR3", "#sigma_{R,3}", 300., 0., 3000., "GeV/c^{2}");
  sigmaR3.setPlotLabel(sigmaR3.GetTitle());
  RooRealVar mean3("mean3", "<m>_{3}", 250., 0., 3000., "GeV/c^{2}");
  mean3.setPlotLabel(mean3.GetTitle());

  // RooArgList coefs("coefs");
  // for(int c = 1; c <= 14; ++c) {
  //   RooRealVar cn(TString::Format("c%i", c), TString::Format("c_{%i}", c),
  // 		  0., -10., 10.);
  //   cn.setPlotLabel(cn.GetTitle());
  //   cn.setError(0.01);
  //   coefs.addClone(cn);
  // }

  // coefs.Print("v");

  RooBifurGauss g1("g1", "g1", mVjj, mean1, sigmaL1, sigmaR1);
  RooBifurGauss g2("g2", "g2", mVjj, mean2, sigmaL2, sigmaR2);
  RooBifurGauss g3("g3", "g3", mVjj, mean2, sigmaL2, sigmaR3);

  RooRealVar fg1("fg1", "f_{g1}", 0.4, 0., 1.);
  fg1.setPlotLabel(fg1.GetTitle());
  RooRealVar fg2("fg2", "f_{g2}", 0.65, 0., 1.);
  fg2.setPlotLabel(fg2.GetTitle());
  
  RooAddPdf sum("sum", "sum", RooArgList(g1, g2,g3), 
  		RooArgList(fg1,fg2), true);
  // RooChebyshevPDF sum("sum", "sum", mVjj, coefs);
  RooWorkspace ws("ws", "ws");
  ws.import(sum, RooFit::RenameAllNodes("low"), RooFit::RenameAllVariablesExcept("low", "MassV2j_PFCor"));
  ws.import(sum, RooFit::RenameAllNodes("high"), RooFit::RenameAllVariablesExcept("high", "MassV2j_PFCor"), RooFit::RecycleConflictNodes());
  ws.Print();

  RooAbsPdf * lowPdf = ws.pdf("sum_low");
  RooAbsPdf * highPdf = ws.pdf("sum_high");

  RooArgSet * lowparams = lowPdf->getParameters(RooArgSet(mVjj,channel));
  lowparams->readFromFile("initialLowParams.txt");
  
  TH1 * lowH = wjetsData.createHistogram("lowH", mVjj, RooFit::Cut(sideband1));
  lowH->Print();
  RooDataSet lowData("lowData", "lowData", RooArgSet(mVjj),
  		     RooFit::Import(wjetsData), RooFit::Cut(sideband1));
  RooDataHist lowDhist("lowDhist", "lowDhist", RooArgList(mVjj), lowH);
  RooHistPdf lowHist("lowHist", "lowHist", RooArgSet(mVjj), lowDhist);

  RooFitResult * fr = 0;
  if (useHists)
    lowPdf = &lowHist;
  else
    fr = lowPdf->fitTo(lowData, RooFit::Minos(false),
		       RooFit::Save(true));
  lowPdf->Print();

  RooArgSet * highparams = highPdf->getParameters(RooArgSet(mVjj,channel));
  highparams->readFromFile("initialHighParams.txt");
  
  TH1 * highH = wjetsData.createHistogram("highH", mVjj, 
					  RooFit::Cut(sideband2));
  highH->Print();
  RooDataSet highData("highData", "highData", RooArgSet(mVjj),
		      RooFit::Import(wjetsData), RooFit::Cut(sideband2));
  RooDataHist highDhist("highDhist", "highDhist", RooArgList(mVjj), highH);
  RooHistPdf highHist("highHist", "highHist", RooArgSet(mVjj), highDhist);
  
  RooFitResult * fr2 = 0;
  if (useHists)
    highPdf = &highHist;
  else
    fr2 = highPdf->fitTo(highData, RooFit::Minos(false),
			 RooFit::Save(true));
  highPdf->Print();

  double chi2 = -1.0;
  if (fitAll) {
    lowparams->setAttribAll("Constant", true);
    highparams->setAttribAll("Constant", true);
    RooRealVar ghigh("ghigh", "ghigh", 0.5, 0., 1.0);
    RooAddPdf all("all", "all", RooArgList(*highPdf,*lowPdf),
		  RooArgList(ghigh));
    RooFitResult * fr_all = all.fitTo(wjetsData, RooFit::Minos(false),
				      RooFit::Save(true));

    RooPlot * mf = mVjj.frame(RooFit::Name("mVjj_plot"));
    wjetsData.plotOn(mf, RooFit::Name("data_pts"));
    all.plotOn(mf, RooFit::Name("wjets_shape"));
    chi2 = mf->chiSquare(1);

    RooHist * resHist = mf->residHist("data_pts", "wjets_shape");

    TCanvas * c1 = new TCanvas("shape", "shape");
    c1->cd();
    mf->Draw();
    TCanvas * c2 = new TCanvas("residuals", "residuals");
    c2->cd();
    resHist->Draw("ap");
    c2->SetGridy();
    resHist->GetXaxis()->SetTitle(mVjj.GetTitle());
    resHist->GetYaxis()->SetTitle("Fit residuals");
    c2->Modified();

    c1->Print("mVjj-MC-template-fit.png");
    c1->Print("mVjj-MC-template-fit.pdf");
    c1->Print("mVjj-MC-template-fit.eps");
    c2->Print("mVjj-MC-template-fit-residuals.png");
    c2->Print("mVjj-MC-template-fit-residuals.pdf");
    c2->Print("mVjj-MC-template-fit-residuals.eps");

    fr_all->Print("v");
    std::cout << "full chi^2/ndf: " << chi2;
    double prob = TMath::Prob(chi2*mVjj.getBins()-1, mVjj.getBins()-1);
    std::cout << " prob: " << prob << '\n';
    delete fr_all;
  }
  
  RooPlot * mf_low = mVjj.frame(RooFit::Name("mVjj_low_plot"));
  lowData.plotOn(mf_low, RooFit::Name("data_pts"));
  lowPdf->plotOn(mf_low, RooFit::Name("wjets_shape"));
  //sum.paramOn(mf_low);
  //g2.plotOn(mf_low);

  RooHist * resHist_low = mf_low->residHist("data_pts", "wjets_shape");
  

  TCanvas * c3 = new TCanvas("shape_low", "shape_low");
  c3->cd();
  mf_low->Draw();
  TCanvas * c4 = new TCanvas("residuals_low", "residuals_low");
  c4->cd();
  resHist_low->Draw("ap");
  c4->SetGridy();
  resHist_low->GetXaxis()->SetTitle(mVjj.GetTitle());
  resHist_low->GetYaxis()->SetTitle("Fit residuals");
  c4->Modified();

  RooPlot * mf_high = mVjj.frame(RooFit::Name("mVjj_high_plot"));
  highData.plotOn(mf_high, RooFit::Name("data_pts"));
  highPdf->plotOn(mf_high, RooFit::Name("wjets_shape"));
  //sum.paramOn(mf_high);
  //g2.plotOn(mf_high);

  RooHist * resHist_high = mf_high->residHist("data_pts", "wjets_shape");
  

  TCanvas * c5 = new TCanvas("shape_high", "shape_high");
  c5->cd();
  mf_high->Draw();
  TCanvas * c6 = new TCanvas("residuals_high", "residuals_high");
  c6->cd();
  resHist_high->Draw("ap");
  c6->SetGridy();
  resHist_high->GetXaxis()->SetTitle(mVjj.GetTitle());
  resHist_high->GetYaxis()->SetTitle("Fit residuals");
  c6->Modified();

  // c1->Print("MassV2j_PFCor-MC-template-fit.png");
  // c1->Print("MassV2j_PFCor-MC-template-fit.pdf");
  // c1->Print("MassV2j_PFCor-MC-template-fit.eps");
  // c2->Print("MassV2j_PFCor-MC-template-fit-residuals.png");
  // c2->Print("MassV2j_PFCor-MC-template-fit-residuals.pdf");
  // c2->Print("MassV2j_PFCor-MC-template-fit-residuals.eps");


  lowparams->writeToFile("lastWjjLowTemplateFit.txt");
  highparams->writeToFile("lastWjjHighTemplateFit.txt");

  if (fr) {
    fr->Print();
    delete fr;
  }
  if (fr2) {
    fr2->Print();
    delete fr2;
  }
  delete lowparams;
  delete highparams;
  delete lowH;
  delete highH;
}
