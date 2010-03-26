// -*- C++ -*-
//
// Package:    HOSiPMAnalysis
// Class:      HOSiPMAnalysis
// 
/**\class HOSiPMAnalysis HOSiPMAnalysis.cc HOSiPM/HOSiPMAnalysis/src/HOSiPMAnalysis.cc

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

#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalRecHit/interface/HORecHit.h"
#include "DataFormats/HcalDigi/interface/HODataFrame.h"

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
// #include "SimDataFormats/CaloTest/interface/HcalTestNumbering.h"

#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"

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

class HOSiPMAnalysis : public edm::EDAnalyzer {
public:
  explicit HOSiPMAnalysis(const edm::ParameterSet&);
  ~HOSiPMAnalysis();
  
  static bool isChannelDead(HcalDetId const& id);

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  // ----------member data ---------------------------
  TString outfname;
  TH1D * NoiseHits;
  TH1D * SignalHits;
  TH1D * outer3byHits;
  TH1D * SigSimHits;
  TH1D * SumSimHits;
  TH1D * OtherHOSimHits;
  TH1D * SigDigisQ;
  TH1D * OtherHODigisQ;
  TTree * hohits;
  TFile * outf;
  double Barrel3x3E;
  double HOSimHitE;
  double HO3x3E;
  double sumHOE;
  double centralRecHit;
  double centralE;
  int maxEta;
  int maxPhi;
  int genAEta;
  int genAPhi;
  double assocSimE;
  double sumsimhits;
  double muonHOE;
  double muonHOES9;
  int Nmu;
  double assocMuSimE;
  int muonEta;
  int muonPhi;
  double assocMuRecE;
  double genp;
  double geneta;
  double genphi;
  double genpt;
  double centralAdc;
  double sumAdc;
  int isSig;
  double mipE;
  bool doFit;
  int centralEta;
  int centralPhi;

  bool doMuons;
  bool findCenter;

  TrackAssociatorParameters assocParams;
  TrackDetectorAssociator assoc;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

bool HOSiPMAnalysis::isChannelDead(HcalDetId const& id) {
  static uint32_t deadChannels[26] = { 0x4601028Bu,
				       0x4601028Cu,
				       0x4601028Du,
				       0x4601028Eu,
				       0x46010225u,
				       0x46010226u,
				       0x460104A6u,
				       0x46010526u,
				       0x460104A7u,
				       0x46010527u,
				       0x460104A8u,
				       0x46010528u,
				       0x46010529u,
				       0x4601052Au,
				       0x4601052Bu,
				       0x460100BBu,
				       0x4601013Bu,
				       0x460101BBu,
				       0x46012292u,
				       0x46012293u,
				       0x460122A3u,
				       0x46012323u,
				       0x460122A4u,
				       0x46012324u,
				       0x460122A5u,
				       0x46012325u};
  static std::vector<uint32_t> deadIds(26);
  static bool inited = false;
  if (!inited) {
    for (int i=0; i<26; ++i) deadIds.push_back(deadChannels[i]);
    std::sort(deadIds.begin(), deadIds.end());
    inited = true;
  }
  std::vector<uint32_t>::const_iterator found = 
    std::find(deadIds.begin(), deadIds.end(), id.rawId());
  if ((found != deadIds.end()) && (*found == id.rawId())) return true;
  else return false;
}

//
// constructors and destructor
//
HOSiPMAnalysis::HOSiPMAnalysis(const edm::ParameterSet& iConfig) :
  outfname(iConfig.getUntrackedParameter<std::string>("outfname", 
						      "HOSiPMMuons.root")),
  NoiseHits(0), SignalHits(0), hohits(0), outf(0), Barrel3x3E(0), 
  HOSimHitE(0), HO3x3E(0), isSig(0),
  mipE(iConfig.getUntrackedParameter<double>("mipE", 1.)),
  doFit(iConfig.getUntrackedParameter<bool>("doFit", false)),
  centralEta(iConfig.getUntrackedParameter<int>("centralEta", 8)),
  centralPhi(iConfig.getUntrackedParameter<int>("centralPhi", 1)),
  doMuons(iConfig.getUntrackedParameter<bool>("doMuons", false)),
  findCenter(iConfig.getUntrackedParameter<bool>("findCenter", doMuons)),
  assocParams(iConfig.getParameter<edm::ParameterSet>("TrackAssociatorParameters"))
{
   //now do what ever initialization is needed
  assoc.useDefaultPropagator();
  outf = new TFile(outfname, "recreate");
  TH1::AddDirectory(true);
  outf->cd();
  NoiseHits = new TH1D("NoiseHits", "HO NoiseHits", 200, -10., 200.);
  SignalHits = new TH1D("SignalHits", "HO SignalHits", 200, -2., 200.);
  outer3byHits = new TH1D("outer3byHits", "HO outer 3x3 hits", 200, -2., 50.);
  SigSimHits = new TH1D("SigSimHits", "Signal Sim Hit Energy", 100, 0., 0.01);
  SumSimHits = new TH1D("SumSimHits", "Sum Sim Hit Energy", 100, 0., 0.1);
  OtherHOSimHits = new TH1D("OtherHOSimHits", "Other HO Sim Hit Energy", 
			    100, 0., 0.01);
  SigDigisQ = new TH1D("SigDigisQ", "Signal Digis Q (fC)", 101, -0.25, 200.25);
  OtherHODigisQ = new TH1D("OtherHODigisQ", "Other HO Digis Q (fC)", 101, 
			   -0.5, 200.25);
  hohits = new TTree("hohits", "hohits");
  hohits->Branch("Barrel3x3E", &Barrel3x3E, "Barrel3x3E/D");
  hohits->Branch("HO3x3E", &HO3x3E, "HO3x3E/D");
  hohits->Branch("HOSimHitE", &HOSimHitE, "HOSimHitE/D");
  hohits->Branch("centralRecHit", &centralRecHit, "centralRecHit/D");
  hohits->Branch("centralE", &centralE, "centralE/D");
  hohits->Branch("sumHOE", &sumHOE, "sumHOE/D");
  hohits->Branch("sumsimhits", &sumsimhits, "sumsimhits/D");
  hohits->Branch("maxEta", &maxEta, "maxEta/I");
  hohits->Branch("maxPhi", &maxPhi, "maxPhi/I");
  hohits->Branch("genAEta", &genAEta, "genAEta/I");
  hohits->Branch("genAPhi", &genAPhi, "genAPhi/I");
  hohits->Branch("assocSimE", &assocSimE, "assocSimE/D");
  hohits->Branch("muonEta", &muonEta, "muonEta/D");
  hohits->Branch("muonPhi", &muonPhi, "muonPhi/D");
  hohits->Branch("muonHOE", &muonHOE, "muonHOE/D");
  hohits->Branch("muonHOES9", &muonHOES9, "muonHOES9/D");
  hohits->Branch("Nmu", &Nmu, "Nmu/I");
  hohits->Branch("assocMuSimE", &assocMuSimE, "assocMuSimE/D");
  hohits->Branch("assocMuRecE", &assocMuRecE, "assocMuRecE/D");
  hohits->Branch("genp", &genp, "genp/D");
  hohits->Branch("genpt", &genpt, "genpt/D");
  hohits->Branch("geneta", &geneta, "geneta/D");
  hohits->Branch("genphi", &genphi, "genphi/D");
  hohits->Branch("centralAdc", &centralAdc, "centralAdc/D");
  hohits->Branch("sumAdc", &sumAdc, "sumAdc/D");
  hohits->Branch("isSig", &isSig, "isSig/I");
}


HOSiPMAnalysis::~HOSiPMAnalysis()
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
HOSiPMAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  Handle<reco::GenParticleCollection> genparts;
  iEvent.getByType(genparts);
  reco::GenParticleCollection::const_iterator genPart = genparts->begin();

  genp = genPart->p();
  geneta = genPart->eta();
  genphi = genPart->phi();
  genpt = genPart->pt();

  GlobalPoint genvertex(genPart->vx(), genPart->vy(), genPart->vz());
  GlobalVector genmomentum(genPart->px(), genPart->py(), genPart->pz());

  TrackDetMatchInfo genMatch = assoc.associate(iEvent, iSetup, genmomentum,
					       genvertex, genPart->charge(),
					       assocParams);
  if (fabs(geneta)>1.305) return;

  //std::cout << "eta: " << geneta << " pt: " << genpt << std::endl;

  Handle<PCaloHitContainer> HcalHits;
  iEvent.getByLabel("g4SimHits", "HcalHits", HcalHits);
  if (!HcalHits.isValid()) {
    std::cout << "no SimHits" << std::endl;
    return;
  }

  std::map<uint32_t,double> homap;
  std::map<uint32_t,double>::iterator hitsum;
  int hiEta = centralEta+1, loEta = centralEta-1;
  if (hiEta==0) ++hiEta;
  if (loEta==0) ++loEta;
  int hiPhi = centralPhi+1, loPhi = centralPhi-1;
  if (hiPhi > 72) hiPhi = hiPhi-72;
  if (loPhi < 1) loPhi = 72+loPhi;
  homap.clear();
//   int det, z, depth, eta, phi, layer, etaBin;
  //std::cout << hiEta << " " << loEta << "," << hiPhi << " " << loPhi << std::endl;
  PCaloHitContainer::const_iterator simhit;
  for (simhit = HcalHits->begin(); simhit != HcalHits->end(); ++simhit) {
    std::map<uint32_t,double>::iterator lb = homap.lower_bound(simhit->id());
    if ((lb != homap.end()) && (simhit->id() == lb->first)) {
      lb->second += simhit->energy();
    } else {
      homap.insert(lb, std::map<uint32_t,double>::value_type(simhit->id(),
							     simhit->energy()));
    }
  }

  HOSimHitE = 0;
  centralE = 0;
  sumsimhits = 0.;
  maxEta = 0;
  maxPhi = -1;
  double maxE = 0.;
  for (hitsum = homap.begin(); hitsum != homap.end(); ++hitsum) {
    HcalDetId tmpId(hitsum->first);
    if ((tmpId.depth()==4) && (!isChannelDead(tmpId))){
      if (hitsum->second > maxE) {
	maxE = hitsum->second;
	maxEta = tmpId.ieta();
	maxPhi = tmpId.iphi();
      }
    }
  }

  centralE = maxE;
  if ((findCenter) && (maxEta != 0)) {
    centralEta = maxEta;
    centralPhi = maxPhi;
    hiEta = centralEta+1;
    loEta = centralEta-1;
    if (hiEta==0) ++hiEta;
    if (loEta==0) ++loEta;
    hiPhi = centralPhi+1;
    loPhi = centralPhi-1;
    if (hiPhi > 72) hiPhi = hiPhi-72;
    if (loPhi < 1) loPhi = 72+loPhi;
  }

  muonHOE = 0.;
  muonHOES9 = 0.;
  Nmu = 0;
  TrackDetMatchInfo * muMatch = 0;
  muonEta = 0;
  muonPhi = -1;
  assocMuSimE = 0.;
  assocMuRecE = 0.;
  if (doMuons) {
    edm::Handle<reco::MuonCollection> muons;
    //  iEvent.getByLabel("trackerMuons", muons);
    iEvent.getByLabel("muons", muons);

    for ( reco::MuonCollection::const_iterator muon = muons->begin(); 
	  muon != muons->end(); ++muon ) {
      // barrel Only
      // if (fabs(muon->eta()) > 1.305) continue;
      if (!muon->isTrackerMuon()) continue;
      ++Nmu;
      if (muonHOE == 0) {
	muonHOE = muon->calEnergy().ho;
	muonHOES9 = muon->calEnergy().hoS9;
	if (! muon->track().isNull() ) {
	  reco::Track const * track = muon->track().get();
	  muMatch = new TrackDetMatchInfo(assoc.associate(iEvent,iSetup,
							  *track,assocParams));
	}
      }
    }
  }

  assocSimE = 0.;
  genAEta = 0;
  genAPhi = -1;
  double maxAE = 0.;
  double maxmuE = 0.;
  for (std::vector<DetId>::const_iterator aid = genMatch.crossedHOIds.begin();
       aid != genMatch.crossedHOIds.end(); ++aid) {
    hitsum = homap.find(aid->rawId());
    if ((hitsum != homap.end()) && (aid->rawId()==hitsum->first) &&
	(!isChannelDead(HcalDetId(hitsum->first)))) {
      HcalDetId mId(hitsum->first);
      assocSimE += hitsum->second;
      if ((genAPhi==-1) || (hitsum->second > maxAE)) {
	genAEta = mId.ieta();
	genAPhi = mId.iphi();
	maxAE = hitsum->second;
      }
    }
  }

  if (muMatch) {
    for (std::vector<DetId>::const_iterator aid = muMatch->crossedHOIds.begin();
	 aid != muMatch->crossedHOIds.end(); ++aid) {
      hitsum = homap.find(aid->rawId());
      if ((hitsum != homap.end()) && (aid->rawId()==hitsum->first) &&
	  (!isChannelDead(HcalDetId(hitsum->first)))) {
	HcalDetId mId(hitsum->first);
	assocMuSimE += hitsum->second;
	if ((muonEta==-1) || (hitsum->second > maxmuE)) {
	  muonEta = mId.ieta();
	  muonPhi = mId.iphi();
	  maxmuE = hitsum->second;
	}
      }
    }
  }

  for (hitsum = homap.begin(); hitsum != homap.end(); ++hitsum) {
    HcalDetId tmpId(hitsum->first);
//     if (tmpId.depth()==4)
//       std::cout << "eta: " << tmpId.ieta() 
// 		<< " phi: " << tmpId.iphi()
// 		<< " simhitE: " << hitsum->second
// 		<< " id: " << tmpId.rawId()
// 		<< std::endl;
//     HcalTestNumbering::unpackHcalIndex(simhit->id(), det, z, depth, eta,
//                                        phi, layer);
    if ((tmpId.depth()==4) && (!isChannelDead(tmpId))) {

      if ((tmpId.ieta() == centralEta) && (tmpId.iphi() == centralPhi)) {
	//       std::cout << " energy: " << simhit->energy() << "\n";
	SigSimHits->Fill(hitsum->second);
	sumsimhits += hitsum->second;
	if (!findCenter) {
	  centralE = hitsum->second;
	}
      } else {
	OtherHOSimHits->Fill(hitsum->second, 1./(16.*72-1.));
	sumsimhits += hitsum->second;
      }

      if ((tmpId.ieta() <= hiEta) && (tmpId.ieta() >= loEta)) {
	if ((hiPhi >= loPhi) && 
	    (tmpId.iphi() <= hiPhi) && (tmpId.iphi() >= loPhi)) {
	  HOSimHitE += hitsum->second;
	} else if ((loPhi > hiPhi) && 
		   ((tmpId.iphi() <= hiPhi) || (tmpId.iphi() >= loPhi))) {
	  HOSimHitE += hitsum->second;
	}
      }
    }
  }

  SumSimHits->Fill(sumsimhits);

  Handle<HORecHitCollection> HORecHits;
  iEvent.getByLabel("horeco",HORecHits);
  if (!HORecHits.isValid()) {
    hohits->Fill();
    if (muMatch) delete muMatch;
    return;
  }
  HORecHitCollection::const_iterator hit;
  homap.clear();
  HO3x3E = 0;
  sumHOE = 0.;
  centralRecHit = 0.;
  maxmuE = 0.;
  for (hit = HORecHits->begin(); hit != HORecHits->end(); ++hit) {
    HcalDetId tmpId(hit->id().rawId());
    double hitEnergy = hit->energy()/mipE;

    if (muMatch) {
      for (std::vector<DetId>::const_iterator aid = muMatch->crossedHOIds.begin();
	   aid != muMatch->crossedHOIds.end(); ++aid) {
	if (aid->rawId() == tmpId.rawId()) {
	  assocMuRecE += hitEnergy;
	}
      }
    }

    if ((tmpId.ieta() == maxEta) && (tmpId.iphi() == maxPhi)) {
      centralRecHit += hitEnergy;
    }

    if ((tmpId.ieta() == centralEta) && (tmpId.iphi() == centralPhi)) {
      SignalHits->Fill(hitEnergy);
      isSig = 1;
      HO3x3E += hitEnergy;
//       std::cout << " depth: " << tmpId.depth() << "\n";
    } else if ((tmpId.ieta()>=loEta) && (tmpId.ieta()<=hiEta)) {
      if ((loPhi <= hiPhi) && 
	  (tmpId.iphi() <= hiPhi) && (tmpId.iphi() >= loPhi)) {
	outer3byHits->Fill(hitEnergy, 1./8.);
	isSig = 2;
	HO3x3E += hitEnergy;
      } else if ((loPhi > hiPhi) &&
		 ((tmpId.iphi() <= hiPhi) || (tmpId.iphi() >= loPhi))) {
	outer3byHits->Fill(hitEnergy, 1./8.);
	isSig = 2;
	HO3x3E += hitEnergy;
      } else {
	NoiseHits->Fill(hitEnergy, 1./(16.*72.-9.));
	isSig = 0;
      }
    } else {
      NoiseHits->Fill(hitEnergy,1./(16.*72.-25.));
      isSig = 0;
    }
    sumHOE += hitEnergy;
  }

  if (muMatch) delete muMatch;

  if ((centralE>0.015) && (centralRecHit<2.)) {
    std::cout << "large simhit: " << centralE 
	      << " no correspondingly large rechit: " << centralRecHit
	      << std::endl;
  }

  if ((centralE>0.0015) && (centralE<0.0025) && (centralRecHit>3.)) {
    std::cout << "MIP-like simhit: " << centralE 
	      << " with large rechit: " << centralRecHit
	      << std::endl;
  }

  if ((centralE>0.0015) && (centralE<0.003) && (muonHOE < 0.01)) {
    std::cout << "MIP-like simhit: " << centralE
	      << " no muon associated deposit: " << muonHOE
	      << " S9 energy: " << muonHOES9
	      << std::endl;
  }

  Barrel3x3E = 0.;

  Handle<CaloTowerCollection> towers;
  iEvent.getByType(towers);
  if (!towers.isValid()) {
    hohits->Fill();
    return;
  }
  CaloTowerCollection::const_iterator tower;
  for (tower = towers->begin(); tower != towers->end(); ++tower) {
    if ((tower->ieta() <= hiEta) && (tower->ieta() >= loEta)) {
      if ((loPhi <= hiPhi) &&
	  (tower->iphi() >= loPhi) && (tower->iphi() <= hiPhi)) {
	Barrel3x3E += tower->energy();
// 	if (tower->outerEnergy() > 0)
// 	  std::cout << "HO energy in tower: " << tower->outerEnergy() << std::endl;
      } else if ((loPhi > hiPhi) &&
		 ((tower->iphi() <= hiPhi) || (tower->iphi() >= loPhi))) {
	Barrel3x3E += tower->energy();
// 	if (tower->outerEnergy() > 0)
// 	  std::cout << "HO energy in tower: " << tower->outerEnergy() << std::endl;
      }
    }
  }

  centralAdc = 0.;
  sumAdc = 0.;

  Handle<HODigiCollection> hodigis;
  iEvent.getByType(hodigis);
  if (!hodigis.isValid()) {
    hohits->Fill();
    return;
  }
  HODigiCollection::const_iterator digi;
  homap.clear();
  double nomCharge, ped;
  for (digi = hodigis->begin(); digi != hodigis->end(); ++digi) {
    nomCharge = 0.;
    ped = 0;
    if (digi->presamples() > 0) {
      for (int s = 0; s < digi->presamples(); ++s) {
	nomCharge += digi->sample(s).nominal_fC();
      }
      ped = nomCharge/digi->presamples();
    }
    nomCharge = 0.;
    for (int s = digi->presamples(); (s<digi->size()) && (s<digi->presamples()+4); 
	 ++s) {
      nomCharge += digi->sample(s).nominal_fC()-ped;
    }
    if ((digi->id().ieta() == centralEta) && (digi->id().iphi() == centralPhi)) {
      SigDigisQ->Fill(nomCharge);
      centralAdc += nomCharge;
    } else {
      OtherHODigisQ->Fill(nomCharge, 1/(16.*72.-1.));
      sumAdc += nomCharge;
    }
  }
  hohits->Fill();
//   std::cout << std::endl;
}


// ------------ method called once each job just before starting event loop  ------------
void 
HOSiPMAnalysis::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HOSiPMAnalysis::endJob() {
  if (doFit) {
    RooRealVar en("hitEnergy", "RecHit Energy", -2., 25., "GeV");
    RooCategory sig("isSig", "isSig");
    sig.defineType("signal", 1);
    sig.defineType("noise", 0);
    RooRealVar MPV("MPV", "MPV", SignalHits->GetMean(), -2., 50., "GeV");
    RooRealVar width("width", "width", SignalHits->GetRMS(), 0., 50., "GeV");
    RooDataSet data("data", "data", hohits, RooArgSet(en,sig), "(isSig==1)&&(hitEnergy<25.)");
    data.Print("v");
    RooRealVar mean("mean", "mean", MPV.getVal(), -2., 25., "GeV");
    RooRealVar sigma("sigma", "sigma", width.getVal(), 0., 50., "GeV");
    RooLandau land("land", "land", en, MPV, width);
    RooGaussian g("g", "g", en, mean, sigma);
    
    RooRealVar fg("fg", "fg", 0., 1.);
    RooAddPdf sum("sum", "sum", RooArgList(g,land), RooArgList(fg));

    RooFitResult * fr = sum.fitTo(data, RooFit::Minos(false), 
				  RooFit::Extended(false), RooFit::Save(true));
    fr->Write();
    RooPlot * plot = en.frame(-2., 25., 27);
    data.plotOn(plot);
    sum.plotOn(plot);
    sum.paramOn(plot);
    plot->SetName("datafit");
    //plot->Write();
  }

  outf->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(HOSiPMAnalysis);
