// -*- C++ -*-
//
// Package:    HOAllLayers
// Class:      HOAllLayers
// 
/**\class HOAllLayers HOAllLayers.cc HOSiPM/HOAllLayers/src/HOAllLayers.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  "Jacob Anderson"
//         Created:  Thu Sep  3 09:02:21 CDT 2009
// $Id: HOAllLayers.cc,v 1.3 2010/08/24 18:03:03 andersj Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimDataFormats/CaloTest/interface/HcalTestNumbering.h"

#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HBHEDataFrame.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
// #include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"

#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"

#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

// #include "TH1.h"
// #include "TFile.h"
#include "TTree.h"

// #include "RooGlobalFunc.h"
// #include "RooDataSet.h"
// #include "RooRealVar.h"
// #include "RooCategory.h"
// #include "RooArgSet.h"
// #include "RooLandau.h"
// #include "RooGaussian.h"
// #include "RooAddPdf.h"
// #include "RooPlot.h"
// #include "RooFitResult.h"

//
// class decleration
//

int const depths = 5;

class HOAllLayers : public edm::EDAnalyzer {
public:
  explicit HOAllLayers(const edm::ParameterSet&);
  ~HOAllLayers();
  
  // static bool isChannelDead(HcalDetId const& id);

  static inline double eta2theta( double eta ) {
    return 2.*atan(exp(-eta));
  }

  static inline double eta2sintheta( double eta ) {
    return sin(eta2theta(eta));
  }


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  // ----------member data ---------------------------
  TTree * hohits;
  TTree * jetTree;
  int maxEta;
  int maxPhi;
  double simHitL[20];
  double sumHBsimhits;
  double sumHEsimhits;
  double sumHOsimhits;
  double genp;
  double geneta;
  double genphi;
  double genpt;
  int isSig;
  double mipE;
  double HOThreshold;
  bool doFit;
  bool doJets;
  int centralEta;
  int centralPhi;
  double nominal_fC[depths];
  double nominal_fC4[depths];
  double nominal_fC9[depths];
  double nominal_fC16[depths];
  double nominal_fC25[depths];
  double nominal_fC144[depths];
  double ped_nominal_fC[depths];
  double ped_nominal_fC4[depths];
  double ped_nominal_fC9[depths];
  double ped_nominal_fC16[depths];
  double ped_nominal_fC25[depths];
  double ped_nominal_fC144[depths];
  double ho_fC;
  double ho_fC4;
  double ho_fC9;
  double ho_fC16;
  double ho_fC25;
  double ho_fC144;
  double ped_ho_fC;
  double ped_ho_fC4;
  double ped_ho_fC9;
  double ped_ho_fC16;
  double ped_ho_fC25;
  double ped_ho_fC144;
  double HB_E;
  double HB_E4;
  double HB_E9;
  double HB_E25;
  double HB_E49;
  double HB_E81;
  double EB_E;
  double EB_E4;
  double EB_E9;
  double EB_E25;
  double EB_E49;
  double EB_E81;
  double HO_E;
  double HO_E_ped;
  double HO_E4;
  double HO_E9;
  double HO_E25;
  double HO_E49;
  double HO_E81;

  double genJetEta, genJetPhi, genJetmass, genJetEt, genJetenergy, 
    genJetpt, genJetemEnergy, genJethadEnergy;
  double caloJetEta, caloJetPhi, caloJetmass, caloJetEt, caloJetenergy, 
    caloJetpt, caloJetemEnergy, caloJethadEnergy, caloJetouterEnergy;
  double caloJetWithHOEta, caloJetWithHOPhi, caloJetWithHOmass, 
    caloJetWithHOEt, caloJetWithHOenergy, caloJetWithHOpt, 
    caloJetWithHOemEnergy, caloJetWithHOhadEnergy, caloJetWithHOouterEnergy;

  bool doMuons;
  bool findCenter;
  bool testNumbering;

  // TrackAssociatorParameters assocParams;
  // TrackDetectorAssociator assoc;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

// bool HOAllLayers::isChannelDead(HcalDetId const& id) {
//   static uint32_t deadChannels[26] = { 0x4601028Bu,
// 				       0x4601028Cu,
// 				       0x4601028Du,
// 				       0x4601028Eu,
// 				       0x46010225u,
// 				       0x46010226u,
// 				       0x460104A6u,
// 				       0x46010526u,
// 				       0x460104A7u,
// 				       0x46010527u,
// 				       0x460104A8u,
// 				       0x46010528u,
// 				       0x46010529u,
// 				       0x4601052Au,
// 				       0x4601052Bu,
// 				       0x460100BBu,
// 				       0x4601013Bu,
// 				       0x460101BBu,
// 				       0x46012292u,
// 				       0x46012293u,
// 				       0x460122A3u,
// 				       0x46012323u,
// 				       0x460122A4u,
// 				       0x46012324u,
// 				       0x460122A5u,
// 				       0x46012325u};
//   static std::vector<uint32_t> deadIds(26);
//   static bool inited = false;
//   if (!inited) {
//     for (int i=0; i<26; ++i) deadIds.push_back(deadChannels[i]);
//     std::sort(deadIds.begin(), deadIds.end());
//     inited = true;
//   }
//   std::vector<uint32_t>::const_iterator found = 
//     std::find(deadIds.begin(), deadIds.end(), id.rawId());
//   if ((found != deadIds.end()) && (*found == id.rawId())) return true;
//   else return false;
// }

//
// constructors and destructor
//
HOAllLayers::HOAllLayers(const edm::ParameterSet& iConfig) :
  isSig(1),
  mipE(iConfig.getUntrackedParameter<double>("mipE", 1.)),
  HOThreshold(iConfig.getUntrackedParameter<double>("HOThreshold", -100.)),
  doFit(iConfig.getUntrackedParameter<bool>("doFit", false)),
  doJets(iConfig.getUntrackedParameter<bool>("doJets", false)),
  centralEta(iConfig.getUntrackedParameter<int>("centralEta", 0)),
  centralPhi(iConfig.getUntrackedParameter<int>("centralPhi", -1)),
  doMuons(iConfig.getUntrackedParameter<bool>("doMuons", false)),
  findCenter(iConfig.getUntrackedParameter<bool>("findCenter", true)),
  testNumbering(iConfig.getUntrackedParameter<bool>("testNumbering", true))
  // assocParams(iConfig.getParameter<edm::ParameterSet>("TrackAssociatorParameters"))
{
   //now do what ever initialization is needed
  // assoc.useDefaultPropagator();
  // outf = new TFile(outfname, "recreate");
  // TH1::AddDirectory(true);
  // outf->cd();
  edm::Service<TFileService> fs;
  hohits = fs->make<TTree>("hcalhits", "hcalhits");
  // hohits = new TTree("hcalhits", "hcalhits");
  hohits->Branch("simHitL", simHitL, "simHitL[20]/D");
  hohits->Branch("sumHBsimhits", &sumHBsimhits, "sumHBsimhits/D");
  hohits->Branch("sumHEsimhits", &sumHEsimhits, "sumHEsimhits/D");
  hohits->Branch("sumHOsimhits", &sumHOsimhits, "sumHOsimhits/D");
  hohits->Branch("maxEta", &maxEta, "maxEta/I");
  hohits->Branch("maxPhi", &maxPhi, "maxPhi/I");
  hohits->Branch("genp", &genp, "genp/D");
  hohits->Branch("genpt", &genpt, "genpt/D");
  hohits->Branch("geneta", &geneta, "geneta/D");
  hohits->Branch("genphi", &genphi, "genphi/D");
  hohits->Branch("isSig", &isSig, "isSig/I");
  hohits->Branch("nominal_fC", nominal_fC, "nominal_fC[5]/D");
  hohits->Branch("nominal_fC4", nominal_fC4, "nominal_fC4[5]/D");
  hohits->Branch("nominal_fC9", nominal_fC9, "nominal_fC9[5]/D");
  hohits->Branch("nominal_fC16", nominal_fC16, "nominal_fC16[5]/D");
  hohits->Branch("nominal_fC25", nominal_fC25, "nominal_fC25[5]/D");
  hohits->Branch("nominal_fC144", nominal_fC144, "nominal_fC144[5]/D");
  hohits->Branch("ped_nominal_fC", ped_nominal_fC, "nominal_fC[5]/D");
  hohits->Branch("ped_nominal_fC4", ped_nominal_fC4, "nominal_fC4[5]/D");
  hohits->Branch("ped_nominal_fC9", ped_nominal_fC9, "nominal_fC9[5]/D");
  hohits->Branch("ped_nominal_fC16", ped_nominal_fC16, "nominal_fC16[5]/D");
  hohits->Branch("ped_nominal_fC25", ped_nominal_fC25, "nominal_fC25[5]/D");
  hohits->Branch("ped_nominal_fC144", ped_nominal_fC144, "nominal_fC144[5]/D");
  hohits->Branch("ho_fC", &ho_fC, "ho_fC/D");
  hohits->Branch("ho_fC4", &ho_fC4, "ho_fC4/D");
  hohits->Branch("ho_fC9", &ho_fC9, "ho_fC9/D");
  hohits->Branch("ho_fC16", &ho_fC16, "ho_fC16/D");
  hohits->Branch("ho_fC25", &ho_fC25, "ho_fC25/D");
  hohits->Branch("ho_fC144", &ho_fC144, "ho_fC144/D");
  hohits->Branch("ped_ho_fC", &ped_ho_fC, "ho_fC/D");
  hohits->Branch("ped_ho_fC4", &ped_ho_fC4, "ho_fC4/D");
  hohits->Branch("ped_ho_fC9", &ped_ho_fC9, "ho_fC9/D");
  hohits->Branch("ped_ho_fC16", &ped_ho_fC16, "ho_fC16/D");
  hohits->Branch("ped_ho_fC25", &ped_ho_fC25, "ho_fC25/D");
  hohits->Branch("ped_ho_fC144", &ped_ho_fC144, "ho_fC144/D");
  hohits->Branch("HB_E", &HB_E, "HB_E/D");
  hohits->Branch("HB_E4", &HB_E4, "HB_E4/D");
  hohits->Branch("HB_E9", &HB_E9, "HB_E9/D");
  hohits->Branch("HB_E25", &HB_E25, "HB_E25/D");
  hohits->Branch("HB_E49", &HB_E49, "HB_E49/D");
  hohits->Branch("HB_E81", &HB_E81, "HB_E81/D");
  hohits->Branch("EB_E", &EB_E, "EB_E/D");
  hohits->Branch("EB_E4", &EB_E4, "EB_E4/D");
  hohits->Branch("EB_E9", &EB_E9, "EB_E9/D");
  hohits->Branch("EB_E25", &EB_E25, "EB_E25/D");
  hohits->Branch("EB_E49", &EB_E49, "EB_E49/D");
  hohits->Branch("EB_E81", &EB_E81, "EB_E81/D");
  hohits->Branch("HO_E", &HO_E, "HO_E/D");
  hohits->Branch("HO_E_ped", &HO_E_ped, "HO_E_ped/D");
  hohits->Branch("HO_E4", &HO_E4, "HO_E4/D");
  hohits->Branch("HO_E9", &HO_E9, "HO_E9/D");
  hohits->Branch("HO_E25", &HO_E25, "HO_E25/D");
  hohits->Branch("HO_E49", &HO_E49, "HO_E49/D");
  hohits->Branch("HO_E81", &HO_E81, "HO_E81/D");

  jetTree = 0;
  if (doJets) {
    jetTree = fs->make<TTree>("jetTree", "jetTree");
    jetTree->Branch("genJetEta", &genJetEta, "genJetEta/D");
    jetTree->Branch("genJetPhi", &genJetPhi, "genJetPhi/D");
    jetTree->Branch("genJetmass", &genJetmass, "genJetmass/D");
    jetTree->Branch("genJetEt", &genJetEt, "genJetEt/D");
    jetTree->Branch("genJetenergy", &genJetenergy, "genJetenergy/D");
    jetTree->Branch("genJetpt", &genJetpt, "genJetpt/D");
    jetTree->Branch("genJetemEnergy", &genJetemEnergy, "genJetemEnergy/D");
    jetTree->Branch("genJethadEnergy", &genJethadEnergy, "genJethadEnergy/D");
    jetTree->Branch("caloJetEta", &caloJetEta, "caloJetEta/D");
    jetTree->Branch("caloJetPhi", &caloJetPhi, "caloJetPhi/D");
    jetTree->Branch("caloJetmass", &caloJetmass, "caloJetmass/D");
    jetTree->Branch("caloJetEt", &caloJetEt, "caloJetEt/D");
    jetTree->Branch("caloJetenergy", &caloJetenergy, "caloJetenergy/D");
    jetTree->Branch("caloJetpt", &caloJetpt, "caloJetpt/D");
    jetTree->Branch("caloJetemEnergy", &caloJetemEnergy, "caloJetemEnergy/D");
    jetTree->Branch("caloJethadEnergy", &caloJethadEnergy, 
		    "caloJethadEnergy/D");
    jetTree->Branch("caloJetouterEnergy", &caloJetouterEnergy, 
		    "caloJetouterEnergy/D");
    jetTree->Branch("caloJetWithHOEta", &caloJetWithHOEta, 
		    "caloJetWithHOEta/D");
    jetTree->Branch("caloJetWithHOPhi", &caloJetWithHOPhi, 
		    "caloJetWithHOPhi/D");
    jetTree->Branch("caloJetWithHOmass", &caloJetWithHOmass, 
		    "caloJetWithHOmass/D");
    jetTree->Branch("caloJetWithHOEt", &caloJetWithHOEt, "caloJetWithHOEt/D");
    jetTree->Branch("caloJetWithHOenergy", &caloJetWithHOenergy, 
		    "caloJetWithHOenergy/D");
    jetTree->Branch("caloJetWithHOpt", &caloJetWithHOpt, "caloJetWithHOpt/D");
    jetTree->Branch("caloJetWithHOemEnergy", &caloJetWithHOemEnergy, 
		    "caloJetWithHOemEnergy/D");
    jetTree->Branch("caloJetWithHOhadEnergy", &caloJetWithHOhadEnergy, 
		    "caloJetWithHOhadEnergy/D");
    jetTree->Branch("caloJetWithHOouterEnergy", &caloJetWithHOouterEnergy, 
		    "caloJetWithHOouterEnergy/D");
    
  }
}

HOAllLayers::~HOAllLayers()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  // outf->Close();
  // delete outf;
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
HOAllLayers::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  Handle<reco::GenParticleCollection> genparts;
  iEvent.getByType(genparts);
  reco::GenParticleCollection::const_iterator genPart = genparts->begin();

  genp = genPart->p();
  geneta = genPart->eta();
  genphi = genPart->phi();
  genpt = genPart->pt();

  Handle<PCaloHitContainer> HcalHits;
  iEvent.getByLabel("g4SimHits", "HcalHits", HcalHits);
  if (!HcalHits.isValid()) {
    std::cout << "no SimHits" << std::endl;
    return;
  }

  std::map<uint32_t,double> homap;
  std::map<uint32_t,double>::iterator hitsum;
  std::map<int, double> towermap;
  homap.clear();
  int det, z, depth, eta, phi, layer;
  PCaloHitContainer::const_iterator simhit;
  sumHBsimhits = 0.;
  sumHEsimhits = 0.;
  sumHOsimhits = 0.;
  for (simhit = HcalHits->begin(); simhit != HcalHits->end(); ++simhit) {
    if (testNumbering) {
      HcalTestNumbering::unpackHcalIndex(simhit->id(), det, z, depth, eta,
					 phi, layer);
    } else {
      HcalDetId tmpId(simhit->id());
      det = tmpId.subdetId();
      z = (tmpId.zside()>0) ? 1 : 0;
      depth = tmpId.depth();
      eta = tmpId.ietaAbs();
      phi = tmpId.iphi();
      layer = depth;
    }
    int tmpIndex = eta*1000 + phi*10 + z;
    if (det == HcalBarrel)
      sumHBsimhits += simhit->energy();
    else if (det == HcalEndcap)
      sumHEsimhits += simhit->energy();
    else if (det == HcalOuter)
      sumHOsimhits += simhit->energy();
    std::map<uint32_t,double>::iterator lb = homap.lower_bound(simhit->id());
    if ((lb != homap.end()) && (simhit->id() == lb->first)) {
      homap[simhit->id()] += simhit->energy();
    } else {
      homap.insert(lb, std::map<uint32_t,double>::value_type(simhit->id(),
							     simhit->energy()));
    }
    std::map<int, double>::iterator lbt = towermap.lower_bound(tmpIndex);
    if ( (lbt != towermap.end()) && (tmpIndex == lbt->first) ) {
      towermap[tmpIndex] += simhit->energy();
    } else {
      towermap.insert(lbt, std::map<int,double>::value_type(tmpIndex,
							    simhit->energy()));
    }
  }

  maxEta = centralEta;
  maxPhi = centralPhi;
  std::map<int,double>::iterator tower;
  double maxE = 0.;
  for (tower = towermap.begin(); tower != towermap.end(); ++tower) {
    // std::cout << " index: " << tower->first << '\n';
    z = tower->first%10;
    phi = (tower->first/10)%100;
    eta = (tower->first/1000)%100;
    // std::cout << " z: " << z
    // 	      << " eta: " << eta
    // 	      << " phi: " << phi
    // 	      << '\n';
    int tmpEta = ((z>0) ? eta : -eta);
    if ((tmpEta == maxEta) && (phi==maxPhi)) {
      maxE += tower->second;
    }
    if (tower->second > maxE) {
      maxE = tower->second;
      maxEta = tmpEta;
      maxPhi = phi;
    }
  }

  std::cout << " maxEta: " << maxEta
  	    << " maxPhi: " << maxPhi
  	    << " maxE: " << maxE
	    << " HB simhits: " << sumHBsimhits
	    << " HO simhits: " << sumHOsimhits
	    << " towers: " << towermap.size()
  	    << '\n';
  if (findCenter) {
    centralEta = maxEta;
    centralPhi = maxPhi;
  }

  for (int l = 0; l<20; ++l)
    simHitL[l] = 0.;
  for (hitsum = homap.begin(); hitsum != homap.end(); ++hitsum) {
    if (testNumbering) {
      HcalTestNumbering::unpackHcalIndex(hitsum->first, det, z, depth, eta,
					 phi, layer);
    } else {
      HcalDetId tmpId(hitsum->first);
      det = tmpId.subdetId();
      z = (tmpId.zside()>0) ? 1 : 0;
      depth = tmpId.depth();
      eta = tmpId.ietaAbs();
      phi = tmpId.iphi();
      layer = depth;
    }
    int tmpEta = (z) ? eta : -eta;
    if ((tmpEta == centralEta) && (phi == centralPhi)) {
      if ( (layer < 20) && (layer >= 0) )
  	simHitL[layer] += hitsum->second;
    }
  }

  for (int i = 0; i<depths; ++i) {
    nominal_fC[i] = 0.;
    nominal_fC4[i] = 0.;
    nominal_fC9[i] = 0.;
    nominal_fC16[i] = 0.;
    nominal_fC25[i] = 0.;
    nominal_fC144[i] = 0.;
    ped_nominal_fC[i] = 0.;
    ped_nominal_fC4[i] = 0.;
    ped_nominal_fC9[i] = 0.;
    ped_nominal_fC16[i] = 0.;
    ped_nominal_fC25[i] = 0.;
    ped_nominal_fC144[i] = 0.;
  }

  ho_fC = 0.;
  ho_fC4 = 0.;
  ho_fC9 = 0.;
  ho_fC16 = 0.;
  ho_fC25 = 0.;
  ho_fC144 = 0.;
  ped_ho_fC = 0.;
  ped_ho_fC4 = 0.;
  ped_ho_fC9 = 0.;
  ped_ho_fC16 = 0.;
  ped_ho_fC25 = 0.;
  ped_ho_fC144 = 0.;

  HB_E = 0.;
  HO_E = 0.;
  EB_E = 0.;
  HB_E4 = 0.;
  HO_E4 = 0.;
  EB_E4 = 0.;
  HB_E9 = 0.;
  HO_E9 = 0.;
  EB_E9 = 0.;
  HB_E25 = 0.;
  HO_E25 = 0.;
  EB_E25 = 0.;
  HB_E49 = 0.;
  HO_E49 = 0.;
  EB_E49 = 0.;
  HB_E81 = 0.;
  HO_E81 = 0.;
  EB_E81 = 0.;
  HO_E_ped = 0.;

  Handle<HBHEDigiCollection> hbhedigis;
  iEvent.getByLabel("simHcalUnsuppressedDigis", hbhedigis);
  if (!hbhedigis.isValid()) {
    hohits->Fill();
    std::cout << "no digis." << std::endl;
    return;
  }
  HBHEDigiCollection::const_iterator digi;

  double sum;
  int linearEta, deta, dphi;
  int const samples = 4;
  if (centralEta < 0) ++centralEta;
  for (digi = hbhedigis->begin(); digi != hbhedigis->end(); ++digi) {
    sum = 0.;
    for (int i = digi->presamples(); 
	 (i < digi->size()) && (i < digi->presamples() + samples); ++i) {
      sum += digi->sample(i).nominal_fC();
    }
    linearEta = (digi->id().ieta() < 0) ? 
      digi->id().ieta()+1 : digi->id().ieta();
    deta = centralEta - linearEta;
    dphi = centralPhi - digi->id().iphi();
    if (dphi > 35) {
      dphi -= 72;
    } else if (dphi < -36) {
      dphi += 72;
    }
    if ( (deta < 6) && (deta > -7) && (dphi < 6) && (dphi > -7) &&
	 (digi->id().depth() < depths) ) {
      if ( (deta < 3) && (deta > -3) && (dphi < 3) && (dphi > -3) ) {
	if ( (deta < 2) && (deta > -3) && (dphi < 2) && (dphi > -3) ) {
	  if ( (deta < 2) && (deta > -2) && (dphi < 2) && (dphi > -2) ) {
	    if ( (deta < 1) && (deta > -2) && (dphi < 1) && (dphi > -2) ) {
	      if ( (deta < 1) && (deta > -1) && (dphi < 1) && (dphi > -1) ) {
		nominal_fC[digi->id().depth()] += sum;
	      } // 1x1 group
	      nominal_fC4[digi->id().depth()] += sum;
	    } // 2x2 group
	    nominal_fC9[digi->id().depth()] += sum;
	  } // 3x3 group
	  nominal_fC16[digi->id().depth()] += sum;
	} // 4x4 group
	nominal_fC25[digi->id().depth()] += sum;
      } // 5x5 group
      nominal_fC144[digi->id().depth()] += sum;
    } // 12x12 group
    deta = -1*centralEta - linearEta;
    dphi = 72 - centralPhi - digi->id().iphi();
    if (dphi > 35) {
      dphi -= 72;
    } else if (dphi < -36) {
      dphi += 72;
    }
    if ( (deta < 6) && (deta > -7) && (dphi < 6) && (dphi > -7) &&
	 (digi->id().depth() < depths) ) {
      if ( (deta < 3) && (deta > -3) && (dphi < 3) && (dphi > -3) ) {
	if ( (deta < 2) && (deta > -3) && (dphi < 2) && (dphi > -3) ) {
	  if ( (deta < 2) && (deta > -2) && (dphi < 2) && (dphi > -2) ) {
	    if ( (deta < 1) && (deta > -2) && (dphi < 1) && (dphi > -2) ) {
	      if ( (deta < 1) && (deta > -1) && (dphi < 1) && (dphi > -1) ) {
		ped_nominal_fC[digi->id().depth()] += sum;
	      } // 1x1 group
	      ped_nominal_fC4[digi->id().depth()] += sum;
	    } // 2x2 group
	    ped_nominal_fC9[digi->id().depth()] += sum;
	  } // 3x3 group
	  ped_nominal_fC16[digi->id().depth()] += sum;
	} // 4x4 group
	ped_nominal_fC25[digi->id().depth()] += sum;
      } // 5x5 group
      ped_nominal_fC144[digi->id().depth()] += sum;
    } // 12x12 group
  }

  Handle<HODigiCollection> hodigis;
  iEvent.getByLabel("simHcalUnsuppressedDigis", hodigis);
  if (!hodigis.isValid()) {
    hohits->Fill();
    std::cout << "no digis." << std::endl;
    return;
  }
  HODigiCollection::const_iterator hodigi;

  for (hodigi = hodigis->begin(); hodigi != hodigis->end(); ++hodigi) {
    sum = 0.;
    for (int i = hodigi->presamples(); 
	 (i < hodigi->size()) && (i < hodigi->presamples() + samples); ++i) {
      sum += hodigi->sample(i).nominal_fC();
    }
    linearEta = (hodigi->id().ieta() < 0) ? 
      hodigi->id().ieta()+1 : hodigi->id().ieta();
    deta = centralEta - linearEta;
    dphi = centralPhi - hodigi->id().iphi();
    if (dphi > 35) {
      dphi -= 72;
    } else if (dphi < -36) {
      dphi += 72;
    }
    if ( (deta < 6) && (deta > -7) && (dphi < 6) && (dphi > -7) ) {
      if ( (deta < 3) && (deta > -3) && (dphi < 3) && (dphi > -3) ) {
	if ( (deta < 2) && (deta > -3) && (dphi < 2) && (dphi > -3) ) {
	  if ( (deta < 2) && (deta > -2) && (dphi < 2) && (dphi > -2) ) {
	    if ( (deta < 1) && (deta > -2) && (dphi < 1) && (dphi > -2) ) {
	      if ( (deta < 1) && (deta > -1) && (dphi < 1) && (dphi > -1) ) {
		ho_fC += sum;
	      } // 1x1 group
	      ho_fC4 += sum;
	    } // 2x2 group
	    ho_fC9 += sum;
	  } // 3x3 group
	  ho_fC16 += sum;
	} // 4x4 group
	ho_fC25 += sum;
      } // 5x5 group
      ho_fC144 += sum;
    } // 12x12 group
    deta = -1*centralEta - linearEta;
    dphi = 72 - centralPhi - hodigi->id().iphi();
    if (dphi > 35) {
      dphi -= 72;
    } else if (dphi < -36) {
      dphi += 72;
    }
    if ( (deta < 6) && (deta > -7) && (dphi < 6) && (dphi > -7) ) {
      if ( (deta < 3) && (deta > -3) && (dphi < 3) && (dphi > -3) ) {
	if ( (deta < 2) && (deta > -3) && (dphi < 2) && (dphi > -3) ) {
	  if ( (deta < 2) && (deta > -2) && (dphi < 2) && (dphi > -2) ) {
	    if ( (deta < 1) && (deta > -2) && (dphi < 1) && (dphi > -2) ) {
	      if ( (deta < 1) && (deta > -1) && (dphi < 1) && (dphi > -1) ) {
		ped_ho_fC += sum;
	      } // 1x1 group
	      ped_ho_fC4 += sum;
	    } // 2x2 group
	    ped_ho_fC9 += sum;
	  } // 3x3 group
	  ped_ho_fC16 += sum;
	} // 4x4 group
	ped_ho_fC25 += sum;
      } // 5x5 group
      ped_ho_fC144 += sum;
    } // 12x12 group
  }
  
//   SumSimHits->Fill(sumsimhits);

  // edm::Handle<HBHERecHitCollection> hbherh;  
  // iEvent.getByLabel("hbhereco",hbherh);
  // if (!hbherh.isValid()) {
  //   std::cout << "no hbhe rechits" << std::endl;
  //   hohits->Fill();
  //   return;
  // }

  edm::Handle<HORecHitCollection> horh;
  iEvent.getByLabel("horeco",horh);
  if (!horh.isValid()) {
    std::cout << "no ho rechits" << std::endl;
    hohits->Fill();
    return;
  }

  HORecHitCollection::const_iterator hit;
  for (hit = horh->begin(); hit != horh->end(); ++hit) {
    HcalDetId tmpId(hit->id().rawId());
    linearEta = (tmpId.ieta() < 0) ? 
      tmpId.ieta()+1 : tmpId.ieta();
    deta = centralEta - linearEta;
    dphi = centralPhi - tmpId.iphi();
    if (dphi > 35) {
      dphi -= 72;
    } else if (dphi < -36) {
      dphi += 72;
    }
    if ( (deta < 5) && (deta > -5) && (dphi < 5) && (dphi > -5) ) {   
      if ( (deta < 4) && (deta > -4) && (dphi < 4) && (dphi > -4) ) {   
	if ( (deta < 3) && (deta > -3) && (dphi < 3) && (dphi > -3) ) {   
	  if ( (deta < 2) && (deta > -2) && (dphi < 2) && (dphi > -2) ) {
	    if ( (deta < 1) && (deta > -2) && (dphi < 1) && (dphi > -2) ) {
	      if ( (deta < 1) && (dphi > -1) && (dphi < 1) && (dphi > -1) ) {
		HO_E += hit->energy();
	      } // 1x1 gourp
	      if (hit->energy() > HOThreshold)
		HO_E4 += hit->energy();
	    } // 2x2 group
	    if (hit->energy() > HOThreshold)
	      HO_E9 += hit->energy();
	  } // 3x3 group
	  if (hit->energy() > HOThreshold)
	    HO_E25 += hit->energy();
	} // 5x5 group
	if (hit->energy() > HOThreshold)
	  HO_E49 += hit->energy();
      } // 7x7 group
      if (hit->energy() > HOThreshold)
	HO_E81 += hit->energy();
    } // 9x9 group
    if ((tmpId.ieta()==-1*centralEta) && (tmpId.iphi()==72-centralPhi)) {
      HO_E_ped += hit->energy();
    }
  }

  Handle<CaloTowerCollection> towers;
  iEvent.getByLabel("towerMaker", towers);
  if (!towers.isValid()) {
    std::cout << "no caloTowers" << std::endl;
    hohits->Fill();
    return;
  }
  CaloTowerCollection::const_iterator ctower;
  for (ctower = towers->begin(); ctower != towers->end(); ++ctower) {
    linearEta = (ctower->ieta() < 0) ? 
      ctower->ieta()+1 : ctower->ieta();
    deta = centralEta - linearEta;
    dphi = centralPhi - ctower->iphi();
    if (dphi > 35) {
      dphi -= 72;
    } else if (dphi < -36) {
      dphi += 72;
    }
    if ( (deta < 5) && (deta > -5) && (dphi < 5) && (dphi > -5) ) {
      if ( (deta < 4) && (deta > -4) && (dphi < 4) && (dphi > -4) ) {
	if ( (deta < 3) && (deta > -3) && (dphi < 3) && (dphi > -3) ) {
	  if ( (deta < 2) && (deta > -2) && (dphi < 2) && (dphi > -2) ) {
	    if ( (deta < 1) && (deta > -2) && (dphi < 1) && (dphi > -2) ) {
	      if ( (deta < 1) && (dphi > -1) && (dphi < 1) && (dphi > -1) ) {
		HB_E += ctower->hadEnergy();
		EB_E += ctower->emEnergy();
	      } // 1x1 group
	      HB_E4 += ctower->hadEnergy();
	      EB_E4 += ctower->emEnergy();
	    } // 2x2 group
	    HB_E9 += ctower->hadEnergy();
	    EB_E9 += ctower->emEnergy();
	  } // 3x3 group
	  HB_E25 += ctower->hadEnergy();
	  EB_E25 += ctower->emEnergy();
	} // 5x5 group
	HB_E49 += ctower->hadEnergy();
	EB_E49 += ctower->emEnergy();
      } // 7x7 group
      HB_E81 += ctower->hadEnergy();
      EB_E81 += ctower->emEnergy();
    } // 9x9 group
  }

  hohits->Fill();
  std::cout.flush();

  if (!doJets) return;

  genJetpt = -10.;
  caloJetpt = -10.;
  caloJetWithHOpt = -10.;

  Handle<reco::GenJetCollection> ak7GenJets;
  iEvent.getByLabel("ak7GenJets", ak7GenJets);
  if (!ak7GenJets.isValid()) {
    std::cout << "no gen jets" << std::endl;
    return;
  }

  reco::GenJetCollection::const_iterator genJet;
  for (genJet = ak7GenJets->begin(); genJet != ak7GenJets->end(); ++genJet) {
    if (genJet->pt() > genJetpt) {
      genJetEta = genJet->eta();
      genJetPhi = genJet->phi();
      genJetmass = genJet->mass();
      genJetEt = genJet->et();
      genJetenergy = genJet->energy();
      genJetpt = genJet->pt();
      genJetemEnergy = genJet->emEnergy();
      genJethadEnergy = genJet->hadEnergy();
    }
  }

  Handle<reco::CaloJetCollection> ak7CaloJets;
  iEvent.getByLabel("ak7CaloJets", ak7CaloJets);
  if (!ak7CaloJets.isValid()) {
    std::cout << "no calo jets" << std::endl;
    return;
  }

  reco::CaloJetCollection::const_iterator caloJet;
  for (caloJet = ak7CaloJets->begin(); caloJet != ak7CaloJets->end(); 
       ++caloJet) {
    if (caloJet->pt() > caloJetpt) {
      caloJetEta = caloJet->eta();
      caloJetPhi = caloJet->phi();
      caloJetmass = caloJet->mass();
      caloJetEt = caloJet->et();
      caloJetenergy = caloJet->energy();
      caloJetpt = caloJet->pt();
      caloJetemEnergy = caloJet->emEnergyInEB();
      caloJethadEnergy = caloJet->hadEnergyInHB();
      caloJetouterEnergy = caloJet->hadEnergyInHO();
    }
  }

  Handle<reco::CaloJetCollection> ak7CaloJetsWithHO;
  iEvent.getByLabel("ak7CaloJetsWithHO", ak7CaloJetsWithHO);
  if (!ak7CaloJetsWithHO.isValid()) {
    std::cout << "no calo jets with HO" << std::endl;
    return;
  }
  
  for (caloJet = ak7CaloJetsWithHO->begin(); 
       caloJet != ak7CaloJetsWithHO->end(); ++caloJet) {
    if (caloJet->pt() > caloJetWithHOpt) {
      caloJetWithHOEta = caloJet->eta();
      caloJetWithHOPhi = caloJet->phi();
      caloJetWithHOmass = caloJet->mass();
      caloJetWithHOEt = caloJet->et();
      caloJetWithHOenergy = caloJet->energy();
      caloJetWithHOpt = caloJet->pt();
      caloJetWithHOemEnergy = caloJet->emEnergyInEB();
      caloJetWithHOhadEnergy = caloJet->hadEnergyInHB();
      caloJetWithHOouterEnergy = caloJet->hadEnergyInHO();
    }
  }

  jetTree->Fill();
  std::cout.flush();
}


// ------------ method called once each job just before starting event loop  ------------
void 
HOAllLayers::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HOAllLayers::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HOAllLayers);
