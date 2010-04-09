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
// $Id: HOSiPMAnalysis.cc,v 1.2 2010/04/08 15:59:16 andersj Exp $
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

#include "HOSiPM/HOSiPMAnalysis/interface/HcalAcceptanceId.h"

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

  // ----------member data ---------------------------
  TString outfname;
  TTree * hohits;
  TFile * outf;
  double Barrel3x3E;
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
  int genHitDeadChannel;
  int genSiPM;
  double assocSimE;
  double sumsimhits;
  double muonHOE;
  double muonHOES9;
  int Nmu;
  int NHOChan;
  double assocMuSimE;
  int muonEta;
  int muonPhi;
  double hoeta;
  double hophi;
  int muAccept;
  int muHitDeadChannel;
  int muSiPM;
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

  TrackAssociatorParameters assocParams;
  TrackDetectorAssociator assoc;

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
  outfname(iConfig.getUntrackedParameter<std::string>("outfname", 
						      "HOSiPMMuons.root")),
  hohits(0), outf(0), Barrel3x3E(0), HOSimHitE(0), HO3x3E(0), isSig(0),
  mipE(iConfig.getUntrackedParameter<double>("mipE", 1.)),
  doFit(iConfig.getUntrackedParameter<bool>("doFit", false)),
  centralEta(iConfig.getUntrackedParameter<int>("centralEta", 8)),
  centralPhi(iConfig.getUntrackedParameter<int>("centralPhi", 1)),
  deta(iConfig.getUntrackedParameter<double>("delta_eta", 0.)),
  dphi(iConfig.getUntrackedParameter<double>("delta_phi", 0.)),
  doMuons(iConfig.getUntrackedParameter<bool>("doMuons", false)),
  findCenter(iConfig.getUntrackedParameter<bool>("findCenter", doMuons)),
  assocParams(iConfig.getParameter<edm::ParameterSet>("TrackAssociatorParameters"))
{
   //now do what ever initialization is needed
  assoc.useDefaultPropagator();
  outf = new TFile(outfname, "recreate");
  TH1::AddDirectory(true);
  outf->cd();
  hohits = new TTree("hohits", "hohits");
  hohits->Branch("Barrel3x3E", &Barrel3x3E, "Barrel3x3E/D");
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
  hohits->Branch("genHitDeadChannel", &genHitDeadChannel, 
		 "genHitDeadChannel/I");
  hohits->Branch("genSiPM", &genSiPM, "genSiPM/I");
  hohits->Branch("assocSimE", &assocSimE, "assocSimE/D");
  hohits->Branch("hoeta", &hoeta, "hoeta/D");
  hohits->Branch("hophi", &hophi, "hophi/D");
  hohits->Branch("muAccept", &muAccept, "muAccept/I");
  hohits->Branch("muHitDeadChannel", &muHitDeadChannel, 
		 "muHitDeadChannel/I");
  hohits->Branch("muSiPM", &muSiPM, "muSiPM/I");
  hohits->Branch("muonEta", &muonEta, "muonEta/I");
  hohits->Branch("muonPhi", &muonPhi, "muonPhi/I");
  hohits->Branch("muonHOE", &muonHOE, "muonHOE/D");
  hohits->Branch("muonHOES9", &muonHOES9, "muonHOES9/D");
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

  if (!HcalAcceptanceId::Inited()) HcalAcceptanceId::initIds(iSetup);

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
  if (HcalAcceptanceId::inGeomAccept(genHOeta,genHOphi,deta,dphi)) 
    genAccept += 1;
  if (HcalAcceptanceId::inNotDeadGeom(genHOeta,genHOphi,deta,dphi))
    genAccept += 10;
  if (HcalAcceptanceId::inSiPMGeom(genHOeta,genHOphi,deta,dphi))
    genAccept += 100;
  genHitDeadChannel = 0;
  genSiPM = 0;
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
      if (HcalAcceptanceId::isChannelDead(tmpId)) ++genHitDeadChannel;
      if (HcalAcceptanceId::isChannelSiPM(tmpId)) ++genSiPM;
    }
  }
  
  centralE = maxE;
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

  std::cout << "max HO sim hits (eta,phi): (" << centralEta << ','
	    << centralPhi << ")\n";

  muonHOE = 0.;
  muonHOES9 = 0.;
  Nmu = 0;
  TrackDetMatchInfo * muMatch = 0;
  muonEta = 0;
  muonPhi = -1;
  assocMuSimE = 0.;
  assocMuRecE = 0.;
  hoeta = -30.;
  hophi = -30.;
  double maxmuE = 0.;
  muAccept = 0;
  muHitDeadChannel = 0;
  NHOChan = 0;
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
	  hophi = muMatch->trkGlobPosAtHO.Phi();
	  hoeta = muMatch->trkGlobPosAtHO.Eta();
	  muEta = track->eta();
	  muPhi = track->phi();
	  mup = track->p();
	  mupt = track->pt();
	  muAccept = 0;
	  if (HcalAcceptanceId::inGeomAccept(hoeta, hophi, deta, dphi)) 
	    muAccept += 1;
	  if (HcalAcceptanceId::inNotDeadGeom(hoeta, hophi, deta, dphi)) 
	    muAccept += 10;
	  if (HcalAcceptanceId::inSiPMGeom(hoeta, hophi, deta, dphi)) 
	    muAccept += 100;
	  muHitDeadChannel = 0;
	  muSiPM = 0;
	  NHOChan = 0;
	  maxmuE = 0.;
	  muonEta = 0;
	  muonPhi = -1;
	  assocMuSimE = 0.;
	  for (std::vector<DetId>::const_iterator aid = 
		 muMatch->crossedHOIds.begin();
	       aid != muMatch->crossedHOIds.end(); ++aid) {
	    HcalDetId mId(aid->rawId());
	    ++NHOChan;
	    if (HcalAcceptanceId::isChannelDead(mId)) ++muHitDeadChannel;
	    if (HcalAcceptanceId::isChannelSiPM(mId)) ++muSiPM;
	    if ( (!HcalAcceptanceId::isChannelDead(mId)) && (muonEta==0) ) {
	      muonEta = mId.ieta();
	      muonPhi = mId.iphi();
	    }
	    hitsum = homap.find(mId.rawId());
	    if ((hitsum != homap.end()) && (mId.rawId() == hitsum->first)) {
	      assocMuSimE += hitsum->second;
	      if ((muonEta==0) || (hitsum->second > maxmuE)) {
		muonEta = mId.ieta();
		muonPhi = mId.iphi();
		maxmuE = hitsum->second;
	      }
	    }
	  }
	  delete muMatch;
	}
      }
    }
  }

  assocSimE = 0.;
  genAEta = 0;
  genAPhi = -1;
  double maxAE = 0.;
  for (std::vector<DetId>::const_iterator aid = genMatch.crossedHOIds.begin();
       aid != genMatch.crossedHOIds.end(); ++aid) {
    HcalDetId mId(aid->rawId());
    if ( (!HcalAcceptanceId::isChannelDead(mId)) && (genAEta==0) ) {
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

  std::cout << "gen associator (eta,phi): (" << genAEta << ','
	    << genAPhi << ")\n"
	    << "muon track associator (eta,phi): (" << muonEta << ','
	    << muonPhi << ")\n";
  Handle<HORecHitCollection> HORecHits;
  iEvent.getByLabel("horeco",HORecHits);
  if (!HORecHits.isValid()) {
    hohits->Fill();
    //if (muMatch) delete muMatch;
    return;
  }
  HORecHitCollection::const_iterator hit;
  homap.clear();
  HO3x3E = 0;
  centralRecHit = 0.;
  maxmuE = 0.;
  for (hit = HORecHits->begin(); hit != HORecHits->end(); ++hit) {
    HcalDetId tmpId(hit->id().rawId());
    double hitEnergy = hit->energy()/mipE;
    if ((tmpId.ieta() == maxEta) && (tmpId.iphi() == maxPhi)) {
      centralRecHit += hitEnergy;
    }
    if ((tmpId.ieta() == centralEta) && (tmpId.iphi() == centralPhi)) {
      isSig = 1;
      HO3x3E += hitEnergy;
    } else if ( (centralPhi > 0) && (tmpId.ieta()>=loEta) && 
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
  }

  //if (muMatch) delete muMatch;

  // if ((centralE>0.015) && (centralRecHit<2.)) {
  //   std::cout << "large simhit: " << centralE 
  // 	      << " no correspondingly large rechit: " << centralRecHit
  // 	      << std::endl;
  // }

  // if ((centralE>0.0015) && (centralE<0.0025) && (centralRecHit>3.)) {
  //   std::cout << "MIP-like simhit: " << centralE 
  // 	      << " with large rechit: " << centralRecHit
  // 	      << std::endl;
  // }

  // if ((centralE>0.0015) && (centralE<0.003) && (muonHOE < 0.01)) {
  //   std::cout << "MIP-like simhit: " << centralE
  // 	      << " no muon associated deposit: " << muonHOE
  // 	      << " S9 energy: " << muonHOES9
  // 	      << std::endl;
  // }

  Barrel3x3E = 0.;

//   Handle<CaloTowerCollection> towers;
//   iEvent.getByType(towers);
//   if (!towers.isValid()) {
//     hohits->Fill();
//     return;
//   }
//   CaloTowerCollection::const_iterator tower;
//   for (tower = towers->begin(); tower != towers->end(); ++tower) {
//     if ((tower->ieta() <= hiEta) && (tower->ieta() >= loEta)) {
//       if ((loPhi <= hiPhi) &&
// 	  (tower->iphi() >= loPhi) && (tower->iphi() <= hiPhi)) {
// 	Barrel3x3E += tower->energy();
// // 	if (tower->outerEnergy() > 0)
// // 	  std::cout << "HO energy in tower: " << tower->outerEnergy() << std::endl;
//       } else if ((loPhi > hiPhi) &&
// 		 ((tower->iphi() <= hiPhi) || (tower->iphi() >= loPhi))) {
// 	Barrel3x3E += tower->energy();
// // 	if (tower->outerEnergy() > 0)
// // 	  std::cout << "HO energy in tower: " << tower->outerEnergy() << std::endl;
//       }
//     }
//   }

  centralAdc = 0.;
  sumAdc = 0.;

  // Handle<HODigiCollection> hodigis;
  // iEvent.getByType(hodigis);
  // if (!hodigis.isValid()) {
  //   hohits->Fill();
  //   return;
  // }
  // HODigiCollection::const_iterator digi;
  // homap.clear();
  // double nomCharge, ped;
  // for (digi = hodigis->begin(); digi != hodigis->end(); ++digi) {
  //   nomCharge = 0.;
  //   ped = 0;
  //   if (digi->presamples() > 0) {
  //     for (int s = 0; s < digi->presamples(); ++s) {
  // 	nomCharge += digi->sample(s).nominal_fC();
  //     }
  //     ped = nomCharge/digi->presamples();
  //   }
  //   nomCharge = 0.;
  //   for (int s = digi->presamples(); (s<digi->size()) && (s<digi->presamples()+4); 
  // 	 ++s) {
  //     nomCharge += digi->sample(s).nominal_fC()-ped;
  //   }
  //   if ((digi->id().ieta() == centralEta) && (digi->id().iphi() == centralPhi)) {
  //     centralAdc += nomCharge;
  //   } else {
  //     sumAdc += nomCharge;
  //   }
  // }
  std::cout.flush();
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
    RooRealVar MPV("MPV", "MPV", 1.0, -2., 50., "GeV");
    RooRealVar width("width", "width", 0.05, 0., 50., "GeV");
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
