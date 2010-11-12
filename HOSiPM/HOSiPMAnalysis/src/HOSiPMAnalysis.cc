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
// $Id: HOSiPMAnalysis.cc,v 1.5 2010/04/15 15:13:51 andersj Exp $
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

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalRecHit/interface/HORecHit.h"
#include "DataFormats/HcalDigi/interface/HODataFrame.h"

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"

#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"

// #include "TH1.h"
#include "TTree.h"

#include "RecoMuon/MuonIdentification/interface/MuonHOAcceptance.h"

//
// class decleration
//

class HOSiPMAnalysis : public edm::EDAnalyzer {
public:
  explicit HOSiPMAnalysis(const edm::ParameterSet&);
  ~HOSiPMAnalysis();
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void resetBranches() ;
  virtual void fillMC(edm::Event const& iEvent, edm::EventSetup const& iSetup) ;
  virtual void fillMuons(edm::Event const& iEvent, 
			 edm::EventSetup const& iSetup) ;

  // ----------member data ---------------------------
  TTree * hohits;
  double HOSimHitE;
  double HO3x3E;
  double centralRecHit;
  double centralE;
  int estEta;
  int maxEta;
  int maxPhi;
  int genAEta;
  int genAPhi;
  double genHOeta;
  double genHOphi;
  int genAccept;
  double assocSimE;
  double sumsimhits;
  double muonHOE;
  double muonHOES9;
  double muonHOmax;
  int Nmu;
  int NHOChan;
  double assocMuSimE;
  int muonEta;
  int muonPhi;
  double hoeta;
  double hophi;
  int muAccept;
  double muEta;
  double muPhi;
  double mup;
  double mupt;
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

  double deta;
  double dphi;

  bool doMuons;
  bool findCenter;
  bool doMC;

  TrackAssociatorParameters assocParams;
  TrackDetectorAssociator assoc;

  //MuonHOAcceptance accept;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HOSiPMAnalysis::HOSiPMAnalysis(const edm::ParameterSet& iConfig) :
  // outfname(iConfig.getUntrackedParameter<std::string>("outfname", 
  // 						      "HOSiPMMuons.root")),
  hohits(0),/* outf(0),*/ HOSimHitE(0), HO3x3E(0),
  mipE(iConfig.getUntrackedParameter<double>("mipE", 1.)),
  doFit(iConfig.getUntrackedParameter<bool>("doFit", false)),
  centralEta(iConfig.getUntrackedParameter<int>("centralEta", 8)),
  centralPhi(iConfig.getUntrackedParameter<int>("centralPhi", 1)),
  deta(iConfig.getUntrackedParameter<double>("delta_eta", 0.)),
  dphi(iConfig.getUntrackedParameter<double>("delta_phi", 0.)),
  doMuons(iConfig.getUntrackedParameter<bool>("doMuons", false)),
  findCenter(iConfig.getUntrackedParameter<bool>("findCenter", doMuons)),
  doMC(iConfig.getUntrackedParameter<bool>("doMC", false)),
  assocParams(iConfig.getParameter<edm::ParameterSet>("TrackAssociatorParameters"))
{
  //now do what ever initialization is needed
  // std::cout << "This is the contructor." << std::endl;
  assoc.useDefaultPropagator();
  // outf = new TFile(outfname, "recreate");
  edm::Service<TFileService> fs;
  hohits = fs->make<TTree>("hohits", "hohits");
  hohits->Branch("HO3x3E", &HO3x3E, "HO3x3E/D");
  hohits->Branch("HOSimHitE", &HOSimHitE, "HOSimHitE/D");
  hohits->Branch("centralRecHit", &centralRecHit, "centralRecHit/D");
  hohits->Branch("centralE", &centralE, "centralE/D");
  hohits->Branch("sumsimhits", &sumsimhits, "sumsimhits/D");
  hohits->Branch("estEta", &estEta, "estEta/I");
  hohits->Branch("maxEta", &maxEta, "maxEta/I");
  hohits->Branch("maxPhi", &maxPhi, "maxPhi/I");
  hohits->Branch("genAEta", &genAEta, "genAEta/I");
  hohits->Branch("genAPhi", &genAPhi, "genAPhi/I");
  hohits->Branch("genHOeta", &genHOeta, "genHOeta/D");
  hohits->Branch("genHOphi", &genHOphi, "genHOphi/D");
  hohits->Branch("genAccept", &genAccept, "genAccept/I");
  hohits->Branch("assocSimE", &assocSimE, "assocSimE/D");
  hohits->Branch("hoeta", &hoeta, "hoeta/D");
  hohits->Branch("hophi", &hophi, "hophi/D");
  hohits->Branch("muAccept", &muAccept, "muAccept/I");
  hohits->Branch("muonEta", &muonEta, "muonEta/I");
  hohits->Branch("muonPhi", &muonPhi, "muonPhi/I");
  hohits->Branch("muonHOE", &muonHOE, "muonHOE/D");
  hohits->Branch("muonHOES9", &muonHOES9, "muonHOES9/D");
  hohits->Branch("muonHOmax", &muonHOmax, "muonHOmax/D");
  hohits->Branch("muEta", &muEta, "muEta/D");
  hohits->Branch("muPhi", &muPhi, "muPhi/D");
  hohits->Branch("mup", &mup, "mup/D");
  hohits->Branch("mupt", &mupt, "mupt/D");
  hohits->Branch("Nmu", &Nmu, "Nmu/I");
  hohits->Branch("NHOChan", &NHOChan, "NHOChan/I");
  hohits->Branch("assocMuSimE", &assocMuSimE, "assocMuSimE/D");
  hohits->Branch("assocMuRecE", &assocMuRecE, "assocMuRecE/D");
  hohits->Branch("genp", &genp, "genp/D");
  hohits->Branch("genpt", &genpt, "genpt/D");
  hohits->Branch("geneta", &geneta, "geneta/D");
  hohits->Branch("genphi", &genphi, "genphi/D");
}


HOSiPMAnalysis::~HOSiPMAnalysis()
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
HOSiPMAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  if ( !MuonHOAcceptance::Inited() ) MuonHOAcceptance::initIds(iSetup);
  resetBranches();

  if (doMC) {
    fillMC(iEvent, iSetup);
  }

  // std::cout << "max HO sim hits (eta,phi): (" << centralEta << ','
  // 	    << centralPhi << ")\n";

  if (doMuons) {
    fillMuons(iEvent, iSetup);
  }

  // std::cout << "gen associator (eta,phi): (" << genAEta << ','
  // 	    << genAPhi << ")\n"
  // 	    << "muon track associator (eta,phi): (" << muonEta << ','
  // 	    << muonPhi << ")\n";
  int hiEta = centralEta+1, loEta = centralEta-1;
  if (hiEta==0) ++hiEta;
  if (loEta==0) ++loEta;
  int hiPhi = centralPhi+1, loPhi = centralPhi-1;
  if (hiPhi > 72) hiPhi = hiPhi-72;
  if (loPhi < 1) loPhi = 72+loPhi;

  Handle<HORecHitCollection> HORecHits;
  iEvent.getByLabel("horeco",HORecHits);
  if (!HORecHits.isValid()) {
    hohits->Fill();
    //if (muMatch) delete muMatch;
    return;
  }
  HORecHitCollection::const_iterator hit;
  HO3x3E = 0;
  centralRecHit = 0.;
  for (hit = HORecHits->begin(); hit != HORecHits->end(); ++hit) {
    HcalDetId tmpId(hit->id().rawId());
    double hitEnergy = hit->energy()/mipE;
    if ((tmpId.ieta() == maxEta) && (tmpId.iphi() == maxPhi)) {
      centralRecHit += hitEnergy;
    }
  }

  std::cout.flush();
  hohits->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
HOSiPMAnalysis::beginJob()
{
}
  
// ------------ method called once each job just after ending the event loop  ------------
void 
HOSiPMAnalysis::endJob() {
}

void HOSiPMAnalysis::resetBranches() {
  HO3x3E = 0.;
  HOSimHitE = 0.;
  centralRecHit = 0.;
  centralE = 0.;
  sumsimhits = 0.;
  estEta = 0;
  maxEta = 0;
  maxPhi = 0;
  genAEta = 0;
  genAPhi = 0;
  genHOeta = -9999;
  genHOphi = -9999;
  genAccept = 0;
  assocSimE = 0.;
  hoeta = -9999.;
  hophi = -9999.;
  muAccept = 0;
  muonEta = 0;
  muonPhi = 0;
  muonHOE = -9999.;
  muonHOES9 = -9999.;
  muonHOmax = -9999.;
  muEta = -9999.;
  muPhi = -9999.;
  mup = -9999.;
  mupt = 9999.;
  Nmu = 0;
  NHOChan = 0;
  assocMuSimE = 0.;
  assocMuRecE = 0.;
  genp = 0.;
  genpt = 0.;
  geneta = -9999.;
  genphi = -9999.;
}

void HOSiPMAnalysis::fillMC(edm::Event const& iEvent, 
			    edm::EventSetup const& iSetup) {

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

  Handle<PCaloHitContainer> HcalHits;
  iEvent.getByLabel("g4SimHits", "HcalHits", HcalHits);
  if (!HcalHits.isValid()) {
    std::cout << "no SimHits" << std::endl;
    return;
  }

  std::map<uint32_t,double> homap;
  std::map<uint32_t,double>::iterator hitsum;
  homap.clear();
//   int det, z, depth, eta, phi, layer, etaBin;
  //std::cout << hiEta << " " << loEta << "," << hiPhi << " " << loPhi << std::endl;
  PCaloHitContainer::const_iterator simhit;
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
      homap[simhit->id()] += simhit->energy();
    } else {
      lb = homap.insert(lb, std::map<uint32_t,double>::value_type(simhit->id(),
				simhit->energy()));
    }
  }

  HOSimHitE = 0;
  centralE = 0;
  sumsimhits = 0.;
  estEta = int(geneta/0.087) + ((geneta>0) ? 1 : -1);
  maxEta = 0;
  //std::cout << "central Eta estimate: " << estEta << '\n';
  maxPhi = -1;
  double maxE = 0.;
  genHOeta = genMatch.trkGlobPosAtHO.Eta();
  genHOphi = genMatch.trkGlobPosAtHO.Phi();
  genAccept = 0;
  if (MuonHOAcceptance::inGeomAccept(genHOeta,genHOphi,deta,dphi)) 
    genAccept += 1;
  if (MuonHOAcceptance::inNotDeadGeom(genHOeta,genHOphi,deta,dphi))
    genAccept += 10;
  if (MuonHOAcceptance::inSiPMGeom(genHOeta,genHOphi,deta,dphi))
    genAccept += 100;
  for (hitsum = homap.begin(); hitsum != homap.end(); ++hitsum) {
    HcalDetId tmpId(hitsum->first);
    if (tmpId.subdet() == HcalOuter) {
      sumsimhits += hitsum->second;
      if ((tmpId.ieta() == maxEta) && (tmpId.iphi() == maxPhi)) {
	maxE += hitsum->second;
      }
      if (hitsum->second > maxE) {
	maxE = hitsum->second;
	maxEta = tmpId.ieta();
	maxPhi = tmpId.iphi();
      }
    }
  }
  
  centralE = maxE;
  if (findCenter) {
    centralEta = maxEta;
    centralPhi = maxPhi;
  }

  assocSimE = 0.;
  genAEta = 0;
  genAPhi = -1;
  double maxAE = 0.;
  for (std::vector<DetId>::const_iterator aid = genMatch.crossedHOIds.begin();
       aid != genMatch.crossedHOIds.end(); ++aid) {
    HcalDetId mId(aid->rawId());
    if ( (!MuonHOAcceptance::isChannelDead(mId)) && (genAEta==0) ) {
	genAEta = mId.ieta();
	genAPhi = mId.iphi();
    }
    hitsum = homap.find(aid->rawId());
    if ((hitsum != homap.end()) && (aid->rawId()==hitsum->first)) {
      assocSimE += hitsum->second;
      if ((genAPhi==-1) || (hitsum->second > maxAE)) {
	genAEta = mId.ieta();
	genAPhi = mId.iphi();
	maxAE = hitsum->second;
      }
    }
  }

}

void HOSiPMAnalysis::fillMuons(edm::Event const& iEvent, 
			       edm::EventSetup const& iSetup) {
  using namespace edm;

  edm::Handle<reco::MuonCollection> muons;
  //  iEvent.getByLabel("trackerMuons", muons);
  iEvent.getByLabel("muons", muons);

  TrackDetMatchInfo * muMatch = 0;
  for ( reco::MuonCollection::const_iterator muon = muons->begin(); 
	muon != muons->end(); ++muon ) {
    // barrel Only
    if (fabs(muon->eta()) > 1.305) continue;
    if (!muon->isTrackerMuon()) continue;
    ++Nmu;
    if (muonHOmax <= 0) {
      if (! muon->track().isNull() ) {
	reco::Track const * track = muon->track().get();
	muMatch = new TrackDetMatchInfo(assoc.associate(iEvent,iSetup,
							*track,assocParams));

	muonHOE = muMatch->crossedEnergy(TrackDetMatchInfo::HORecHits);
	muonHOES9 = muMatch->nXnEnergy(TrackDetMatchInfo::HORecHits,1);
	DetId hoMaxId = 
	  muMatch->findMaxDeposition(TrackDetMatchInfo::HORecHits,1);
	for(std::vector<const HORecHit*>::const_iterator muhit = 
	      muMatch->hoRecHits.begin(); 	 
	    muhit != muMatch->hoRecHits.end(); ++muhit) { 	 
	  if ((*muhit)->id() != hoMaxId) continue; 	 
	  muonHOmax = (*muhit)->energy();
	}

	hophi = muMatch->trkGlobPosAtHO.Phi();
	hoeta = muMatch->trkGlobPosAtHO.Eta();
	muEta = track->eta();
	muPhi = track->phi();
	mup = track->p();
	mupt = track->pt();
	muAccept = 0;
	if (MuonHOAcceptance::inGeomAccept(hoeta, hophi, deta, dphi)) 
	  muAccept += 1;
	if (MuonHOAcceptance::inNotDeadGeom(hoeta, hophi, deta, dphi)) 
	  muAccept += 10;
	if (MuonHOAcceptance::inSiPMGeom(hoeta, hophi, deta, dphi)) 
	  muAccept += 100;
	NHOChan = 0;
	muonEta = 0;
	muonPhi = -1;
	for (std::vector<DetId>::const_iterator aid = 
	       muMatch->crossedHOIds.begin();
	     aid != muMatch->crossedHOIds.end(); ++aid) {
	  HcalDetId mId(aid->rawId());
	  ++NHOChan;
	  if ( (!MuonHOAcceptance::isChannelDead(mId)) && (muonEta==0) ) {
	    muonEta = mId.ieta();
	    muonPhi = mId.iphi();
	  }
	}
	delete muMatch;
      }
    }
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(HOSiPMAnalysis);
