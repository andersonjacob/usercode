// -*- C++ -*-
//
// Package:    HOWithReco
// Class:      HOWithReco
// 
/**\class HOWithReco HOWithReco.cc HOSiPM/HOWithReco/src/HOWithReco.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  "Jacob Anderson"
//         Created:  Thu Sep  3 09:02:21 CDT 2009
// $Id$
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

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
//#include "SimDataFormats/CaloTest/interface/HcalTestNumbering.h"

#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HORecHit.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
// #include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"

#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"

#include "TH1.h"
#include "TFile.h"
#include "TTree.h"

#include "RooGlobalFunc.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooArgSet.h"
#include "RooLandau.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"

//
// class decleration
//

class HOWithReco : public edm::EDAnalyzer {
public:
  explicit HOWithReco(const edm::ParameterSet&);
  ~HOWithReco();
  
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
  TString outfname;
  TTree * hohits;
  TFile * outf;
  int maxEta;
  int maxPhi;
  int HOmaxEta;
  int HOmaxPhi;
  double maxSimHit;
  double sumsimhits;
  double maxbarrelsimhit;
  double sumbarrelsimhits;
  double HO3x3E;
  double sumHOE;
  double centralHORecHit;
  double HB3x3E;
  double centralHBRecHit;
  double Barrel3x3E;
  double centralBarrelTower;
  double genp;
  double geneta;
  double genphi;
  double genpt;
  int isSig;
  double mipE;
  bool doFit;
  int centralEta;
  int centralPhi;

  bool doMuons;
  bool findCenter;
  bool findWithBarrel;

  // TrackAssociatorParameters assocParams;
  // TrackDetectorAssociator assoc;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

// bool HOWithReco::isChannelDead(HcalDetId const& id) {
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
HOWithReco::HOWithReco(const edm::ParameterSet& iConfig) :
  outfname(iConfig.getUntrackedParameter<std::string>("outfname", 
						      "HOSiPMMuons.root")),
  outf(0), isSig(1),
  mipE(iConfig.getUntrackedParameter<double>("mipE", 1.)),
  doFit(iConfig.getUntrackedParameter<bool>("doFit", false)),
  centralEta(iConfig.getUntrackedParameter<int>("centralEta", 0)),
  centralPhi(iConfig.getUntrackedParameter<int>("centralPhi", -1)),
  doMuons(iConfig.getUntrackedParameter<bool>("doMuons", false)),
  findCenter(iConfig.getUntrackedParameter<bool>("findCenter", true)),
  findWithBarrel(iConfig.getUntrackedParameter<bool>("findWithBarrel", true))
  // assocParams(iConfig.getParameter<edm::ParameterSet>("TrackAssociatorParameters"))
{
   //now do what ever initialization is needed
  // assoc.useDefaultPropagator();
  outf = new TFile(outfname, "recreate");
  TH1::AddDirectory(true);
  outf->cd();
  hohits = new TTree("hohits", "hohits");
  hohits->Branch("HOmaxEta", &HOmaxEta, "HOmaxEta/I");
  hohits->Branch("HOmaxPhi", &HOmaxPhi, "HOmaxPhi/I");
  hohits->Branch("maxSimHit", &maxSimHit, "maxSimHit/D");
  hohits->Branch("sumsimhits", &sumsimhits, "sumsimhits/D");
  hohits->Branch("maxEta", &maxEta, "maxEta/I");
  hohits->Branch("maxPhi", &maxPhi, "maxPhi/I");
  hohits->Branch("maxbarrelsimhit", &maxbarrelsimhit, "maxbarrelsimhit/D");
  hohits->Branch("sumbarrelsimhits", &sumbarrelsimhits, "sumbarrelsimhits/D");
  hohits->Branch("centralHORecHit", &centralHORecHit, "centralHORecHit/D");
  hohits->Branch("HO3x3E", &HO3x3E, "HO3x3E/D");
  hohits->Branch("centralHBRecHit", &centralHBRecHit, "centralHBRecHit/D");
  hohits->Branch("HB3x3E", &HB3x3E, "HB3x3E/D");
  hohits->Branch("sumHOE", &sumHOE, "sumHOE/D");
  hohits->Branch("Barrel3x3E", &Barrel3x3E, "Barrel3x3E/D");
  hohits->Branch("centralBarrelTower", &centralBarrelTower, 
		 "centralBarrelTower/D");
  hohits->Branch("genp", &genp, "genp/D");
  hohits->Branch("genpt", &genpt, "genpt/D");
  hohits->Branch("geneta", &geneta, "geneta/D");
  hohits->Branch("genphi", &genphi, "genphi/D");
  hohits->Branch("isSig", &isSig, "isSig/I");
}


HOWithReco::~HOWithReco()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  outf->Close();
  delete outf;
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
HOWithReco::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  Handle<reco::GenParticleCollection> genparts;
  iEvent.getByType(genparts);
  reco::GenParticleCollection::const_iterator genPart = genparts->begin();

  genp = genPart->p();
  geneta = genPart->eta();
  genphi = genPart->phi();
  genpt = genPart->pt();

  // GlobalPoint genvertex(genPart->vx(), genPart->vy(), genPart->vz());
  // GlobalVector genmomentum(genPart->px(), genPart->py(), genPart->pz());

  // TrackDetMatchInfo genMatch = assoc.associate(iEvent, iSetup, genmomentum,
  // 					       genvertex, genPart->charge(),
  // 					       assocParams);

  Handle<PCaloHitContainer> HcalHits;
  iEvent.getByLabel("g4SimHits", "HcalHits", HcalHits);
  if (!HcalHits.isValid()) {
    std::cout << "no SimHits" << std::endl;
    return;
  }

  std::map<uint32_t,double> homap;
  std::map<uint32_t,double>::iterator hitsum;
  int hiEta = centralEta+1, loEta = centralEta-1;
  int hiPhi = centralPhi+1, loPhi = centralPhi-1;
  homap.clear();
  PCaloHitContainer::const_iterator simhit;
  //int det,z,depth,eta,phi,layer;
  for (simhit = HcalHits->begin(); simhit != HcalHits->end(); ++simhit) {
    HcalDetId det(simhit->id());
    // std::cout << " det: " << det.subdet()
    // 		<< " eta: " << det.ieta()
    // 		<< " phi: " << det.iphi()
    // 		<< " depth: " << det.depth()
    // 		<< " energy: " << simhit->energy()*eta2sintheta(geneta)
    // 		<< "\n";
    std::map<uint32_t,double>::iterator lb = homap.lower_bound(simhit->id());
    if ((lb != homap.end()) && (simhit->id() == lb->first)) {
      homap[simhit->id()] += simhit->energy()*eta2sintheta(geneta);
    } else {
      lb = homap.insert(lb, std::map<uint32_t,double>::value_type(simhit->id(),
				simhit->energy()*eta2sintheta(geneta)));
    }
    // std::cout << " tower: " << lb->first
    // 		<< " sum: " << homap[simhit->id()]
    // 		<< "\n";
  }

  HOmaxEta = 0;
  HOmaxPhi = -1;
  double HOmaxE = 0.;
  sumsimhits = 0.;
  for (hitsum = homap.begin(); hitsum != homap.end(); ++hitsum) {
    HcalDetId det(hitsum->first);
    if (det.subdet() == HcalOuter) {
      sumsimhits += hitsum->second;
      if ((det.ieta() == HOmaxEta) && (det.iphi() == HOmaxPhi)) {
	HOmaxE += hitsum->second;
      }
      if ((hitsum->second > HOmaxE) && (det.ietaAbs()<10)) {
	HOmaxE += hitsum->second;
	HOmaxEta = det.ieta();
	HOmaxPhi = det.iphi();
      }
    }
  }

  std::cout << " outer simHit maxEta: " << HOmaxEta
	    << " maxPhi: " << HOmaxPhi
	    << " maxE: " << HOmaxE
	    << " sumsimhits HO: " << sumsimhits
	    << "\n";
  maxEta = HOmaxEta;
  maxPhi = HOmaxPhi;
  if (findCenter) {
    centralEta = maxEta;
    centralPhi = maxPhi;
    hiEta = centralEta+1;
    loEta = centralEta-1;
    if (hiEta==0) ++hiEta;
    if (loEta==0) --loEta;
    hiPhi = centralPhi+1;
    loPhi = centralPhi-1;
    if (hiPhi > 72) hiPhi = hiPhi-72;
    if (loPhi < 1) loPhi = 72+loPhi;
  }

  maxSimHit = 0.;
  for (hitsum = homap.begin(); hitsum != homap.end(); ++hitsum) {
    HcalDetId det(hitsum->first);
    if (det.subdet() == HcalOuter) {
      if ((det.ieta() == maxEta) && (det.iphi() == maxPhi)) {
	maxSimHit += hitsum->second;
      }
    }
  }

//   SumSimHits->Fill(sumsimhits);

  sumbarrelsimhits = 0.;
  maxbarrelsimhit = 0.;
  if (findWithBarrel) {
    maxEta = 0;
    maxPhi = -1;
    double maxE = 0.;
    for (hitsum = homap.begin(); hitsum != homap.end(); ++hitsum) {
      HcalDetId det(hitsum->first);
      if (det.subdet() == HcalBarrel) {
	sumbarrelsimhits += hitsum->second;
	if ((det.ieta() == maxEta) && (det.iphi() == maxPhi)) {
	  maxE += hitsum->second;
	}
	if ((hitsum->second > maxE) && (det.ietaAbs()<10)) {
	  maxE += hitsum->second;
	  maxEta = det.ieta();
	  maxPhi = det.iphi();
	}
      }
    }

    maxbarrelsimhit = maxE;
    centralEta = maxEta;
    centralPhi = maxPhi;
    hiEta = centralEta+1;
    loEta = centralEta-1;
    if (hiEta==0) ++hiEta;
    if (loEta==0) --loEta;
    hiPhi = centralPhi+1;
    loPhi = centralPhi-1;
    if (hiPhi > 72) hiPhi = hiPhi-72;
    if (loPhi < 1) loPhi = 72+loPhi;
    std::cout << " barrel tower maxEta: " << maxEta
	      << " maxPhi: " << maxPhi
	      << " maxE: " << maxE
	      << " sumsimhits HB: " << sumbarrelsimhits
	      << "\n";
  }

  HO3x3E = 0;
  sumHOE = 0.;
  centralHORecHit = 0.;

  Handle<HORecHitCollection> HORecHits;
  iEvent.getByLabel("horeco",HORecHits);
  if (!HORecHits.isValid()) {
    hohits->Fill();
    std::cout.flush();
    return;
  }
  HORecHitCollection::const_iterator hit;
  homap.clear();
  for (hit = HORecHits->begin(); hit != HORecHits->end(); ++hit) {
    HcalDetId tmpId(hit->id().rawId());
    double hitEnergy = hit->energy()*eta2sintheta(geneta)/mipE;
    if ((tmpId.ieta() == maxEta) && (tmpId.iphi() == maxPhi)) {
      centralHORecHit += hitEnergy;
    }
    if ((tmpId.ieta() == centralEta) && (tmpId.iphi() == centralPhi)) {
      isSig = 1;
      HO3x3E += hitEnergy;
    } else if ( (centralEta != 0) && (tmpId.ieta()>=loEta) && 
		(tmpId.ieta()<=hiEta) ) {
      if ((loPhi <= hiPhi) && 
	  (tmpId.iphi() <= hiPhi) && (tmpId.iphi() >= loPhi)) {
	isSig = 2;
	HO3x3E += hitEnergy;
      } else if ((loPhi > hiPhi) &&
		 ((tmpId.iphi() <= hiPhi) || (tmpId.iphi() >= loPhi))) {
	isSig = 2;
	HO3x3E += hitEnergy;
      } else {
	isSig = 0;
      }
    } else {
      isSig = 0;
    }
    sumHOE += hitEnergy;
  }

  HB3x3E = 0.;
  centralHBRecHit = 0.;

  Handle<HBHERecHitCollection> hbheRecHits;
  iEvent.getByLabel("hbhereco",hbheRecHits);
  if (!hbheRecHits.isValid()) {
    hohits->Fill();
    std::cout.flush();
    return;
  }
  HBHERecHitCollection::const_iterator rech;

  for (rech = hbheRecHits->begin(); rech != hbheRecHits->end(); ++rech) {
    HcalDetId tmpId(rech->id().rawId());
    double rechEnergy = rech->energy()*eta2sintheta(geneta);
    if ((tmpId.ieta() == maxEta) && (tmpId.iphi() == maxPhi)) {
      centralHBRecHit += rechEnergy;
    }
    if ( (centralEta != 0) && (tmpId.ieta() <= hiEta) && 
	 (tmpId.ieta() >= loEta) ) {
      if ((loPhi <= hiPhi) &&
	  (tmpId.iphi() >= loPhi) && (tmpId.iphi() <= hiPhi)) {
	HB3x3E += rechEnergy;
      } else if ((loPhi > hiPhi) &&
		 ((tmpId.iphi() <= hiPhi) || (tmpId.iphi() >= loPhi))) {
	HB3x3E += rechEnergy;
      }
    }
  }

  Barrel3x3E = 0.;
  centralBarrelTower = 0.;

  Handle<CaloTowerCollection> towers;
  iEvent.getByLabel("towerMaker", towers);
  if (!towers.isValid()) {
    hohits->Fill();
    return;
  }
  CaloTowerCollection::const_iterator tower;
  for (tower = towers->begin(); tower != towers->end(); ++tower) {
    double towerEnergy = 
      ( tower->emEnergy()+tower->hadEnergy() )*eta2sintheta(geneta);
    if ( (tower->ieta() == maxEta) && (tower->iphi() == maxPhi) ) {
      centralBarrelTower += towerEnergy;
    }
    if ( (centralEta != 0) && (tower->ieta() <= hiEta) && 
	 (tower->ieta() >= loEta) ) {
      if ((loPhi <= hiPhi) &&
	  (tower->iphi() >= loPhi) && (tower->iphi() <= hiPhi)) {
	Barrel3x3E += towerEnergy;
      } else if ((loPhi > hiPhi) &&
		 ((tower->iphi() <= hiPhi) || (tower->iphi() >= loPhi))) {
	Barrel3x3E += towerEnergy;
      }
    }
  }

  

  hohits->Fill();
  std::cout.flush();
}


// ------------ method called once each job just before starting event loop  ------------
void 
HOWithReco::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HOWithReco::endJob() {
  outf->Write();

  if (doFit) {
    RooRealVar en("sumsimhits", "SimHit Energy", 0., 0.05, "GeV");
    RooRealVar max("maxEta", "max eta", -72, 72);
    RooCategory sig("isSig", "isSig");
    sig.defineType("signal", 1);
    sig.defineType("noise", 0);
    RooRealVar MPV("MPV", "MPV", 0.004, 0.001, 0.05, "GeV");
    RooRealVar width("width", "width", 0.0006, 0., 0.01, "GeV");
    RooDataSet data("data", "data", hohits, RooArgSet(en,max), "maxEta!=0");
    data.Print("v");
    RooRealVar mean("mean", "mean", MPV.getVal(), 0.001, 0.05, "GeV");
    RooRealVar sigma("sigma", "sigma", width.getVal()/2., 0., 0.01, "GeV");
    RooLandau land("land", "land", en, MPV, width);
    RooGaussian g("g", "g", en, mean, sigma);
    
    RooRealVar fg("fg", "fg", 0.2, 0., 1.);
    RooAddPdf sum("sum", "sum", RooArgList(g,land), RooArgList(fg));

    RooFitResult * fr = sum.fitTo(data, RooFit::Minos(false), 
				  RooFit::Range(0.0025, 0.01),
				  RooFit::Extended(false), RooFit::Save(true));
    outf->cd();
    fr->Print("v");
    //double maximum  = sum.maxVal(sum.getMaxVal(RooArgSet(en)));
    //std::cout << "Maximum of fit function: " << maximum << std::endl;
    fr->Write("fitRes");
    RooPlot * plot = en.frame(0.00, 0.015, 30);
    data.plotOn(plot);
    sum.plotOn(plot);
    sum.paramOn(plot);
    plot->SetName("simhitfit");
    plot->Write();
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(HOWithReco);
