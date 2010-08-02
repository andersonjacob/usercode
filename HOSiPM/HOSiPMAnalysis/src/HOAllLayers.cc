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
// $Id: HOAllLayers.cc,v 1.1 2010/03/26 16:02:15 andersj Exp $
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
#include "SimDataFormats/CaloTest/interface/HcalTestNumbering.h"

#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HBHEDataFrame.h"

#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
// #include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"

#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"

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
  TString outfname;
  TTree * hohits;
  TFile * outf;
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
  bool doFit;
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
  outfname(iConfig.getUntrackedParameter<std::string>("outfname", 
						      "HOSiPMMuons.root")),
  outf(0), isSig(1),
  mipE(iConfig.getUntrackedParameter<double>("mipE", 1.)),
  doFit(iConfig.getUntrackedParameter<bool>("doFit", false)),
  centralEta(iConfig.getUntrackedParameter<int>("centralEta", 0)),
  centralPhi(iConfig.getUntrackedParameter<int>("centralPhi", -1)),
  doMuons(iConfig.getUntrackedParameter<bool>("doMuons", false)),
  findCenter(iConfig.getUntrackedParameter<bool>("findCenter", true)),
  testNumbering(iConfig.getUntrackedParameter<bool>("testNumbering", true))
  // assocParams(iConfig.getParameter<edm::ParameterSet>("TrackAssociatorParameters"))
{
   //now do what ever initialization is needed
  // assoc.useDefaultPropagator();
  outf = new TFile(outfname, "recreate");
  TH1::AddDirectory(true);
  outf->cd();
  hohits = new TTree("hcalhits", "hcalhits");
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
}


HOAllLayers::~HOAllLayers()
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
  std::map<int, double> towermap;
  // int hiEta = centralEta+1, loEta = centralEta-1;
  // int hiPhi = centralPhi+1, loPhi = centralPhi-1;
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
      det = tmpId.subdet();
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
    // hiEta = centralEta+1;
    // loEta = centralEta-1;
    // if (hiEta==0) ++hiEta;
    // if (loEta==0) --loEta;
    // hiPhi = centralPhi+1;
    // loPhi = centralPhi-1;
    // if (hiPhi > 72) hiPhi = hiPhi-72;
    // if (loPhi < 1) loPhi = 72+loPhi;
  }

  for (int l = 0; l<20; ++l)
    simHitL[l] = 0.;
  for (hitsum = homap.begin(); hitsum != homap.end(); ++hitsum) {
    if (testNumbering) {
      HcalTestNumbering::unpackHcalIndex(hitsum->first, det, z, depth, eta,
					 phi, layer);
    } else {
      HcalDetId tmpId(hitsum->first);
      det = tmpId.subdet();
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
    deta = -3 - linearEta;
    dphi = 36 - digi->id().iphi();
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
    deta = -3 - linearEta;
    dphi = 36 - hodigi->id().iphi();
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

  hohits->Fill();
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
  outf->Write();

  // if (doFit) {
  //   RooRealVar en("sumsimhits", "SimHit Energy", 0., 0.05, "GeV");
  //   RooRealVar max("maxEta", "max eta", -72, 72);
  //   RooCategory sig("isSig", "isSig");
  //   sig.defineType("signal", 1);
  //   sig.defineType("noise", 0);
  //   RooRealVar MPV("MPV", "MPV", 0.004, 0.001, 0.05, "GeV");
  //   RooRealVar width("width", "width", 0.0006, 0., 0.01, "GeV");
  //   RooDataSet data("data", "data", hohits, RooArgSet(en,max), "maxEta!=0");
  //   data.Print("v");
  //   RooRealVar mean("mean", "mean", MPV.getVal(), 0.001, 0.05, "GeV");
  //   RooRealVar sigma("sigma", "sigma", width.getVal()/2., 0., 0.01, "GeV");
  //   RooLandau land("land", "land", en, MPV, width);
  //   RooGaussian g("g", "g", en, MPV, sigma);
    
  //   RooRealVar fg("fg", "fg", 0.1, 0., 0.5);
  //   RooAddPdf sum("sum", "sum", RooArgList(g,land), RooArgList(fg));

  //   RooFitResult * fr = sum.fitTo(data, RooFit::Minos(false), 
  // 				  RooFit::Range(0.0025, 0.01),
  // 				  RooFit::Extended(false), RooFit::Save(true));
  //   outf->cd();
  //   fr->Print("v");
  //   fr->Write("fitRes");
  //   RooPlot * plot = en.frame(0.00, 0.015, 30);
  //   data.plotOn(plot);
  //   sum.plotOn(plot);
  //   sum.paramOn(plot);
  //   plot->SetName("simhitfit");
  //   plot->Write();
  // }
}

//define this as a plug-in
DEFINE_FWK_MODULE(HOAllLayers);
