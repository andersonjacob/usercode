// -*- C++ -*-
//
// Package:    HcalHOTBPlotAnal
// Class:      HcalHOTBPlotAnal
// 
/**\class HcalHOTBPlotAnal HcalHOTBPlotAnal.cc RecoTBCalo/HcalHOTBPlotAnal/src/HcalHOTBPlotAnal.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Phillip R. Dudero
//         Created:  Tue Jan 16 21:11:37 CST 2007
// $Id: HcalQLPlotAnal.cc,v 1.2 2011/07/21 12:09:22 andersj Exp $
//
//


// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "TBDataFormats/HcalTBObjects/interface/HcalTBTriggerData.h"
#include "TBDataFormats/HcalTBObjects/interface/HcalTBEventPosition.h"
#include "TBDataFormats/HcalTBObjects/interface/HcalTBBeamCounters.h"
#include "TBDataFormats/HcalTBObjects/interface/HcalTBTiming.h"

#include "andersj/HcalPlotter/src/HcalQLPlotAnalAlgos.h"
//
// class declaration
//

class HcalHOTBPlotAnal : public edm::EDAnalyzer {
   public:
      explicit HcalHOTBPlotAnal(const edm::ParameterSet&);
      ~HcalHOTBPlotAnal();


   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
  edm::InputTag hbheRHLabel_,hoRHLabel_,hfRHLabel_;
  edm::InputTag hcalDigiLabel_, hcalTrigLabel_;
  edm::InputTag ebRHLabel_;
  bool doCalib_;
  bool doBeamCounters_;
  double calibFC2GeV_;
  int hbDigiCnt_;
  int hoDigiCnt_;
  int hfDigiCnt_;
  int ebCnt_;
  int evtCnt_;
  HcalQLPlotAnalAlgos * algo_;

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
HcalHOTBPlotAnal::HcalHOTBPlotAnal(const edm::ParameterSet& iConfig) :
  hbheRHLabel_(iConfig.getUntrackedParameter<edm::InputTag>("hbheRHtag")),
  hoRHLabel_(iConfig.getUntrackedParameter<edm::InputTag>("hoRHtag")),
  hfRHLabel_(iConfig.getUntrackedParameter<edm::InputTag>("hfRHtag")),
  hcalDigiLabel_(iConfig.getUntrackedParameter<edm::InputTag>("hcalDigiTag")),
  hcalTrigLabel_(iConfig.getUntrackedParameter<edm::InputTag>("hcalTrigTag")),
  ebRHLabel_(iConfig.getUntrackedParameter<edm::InputTag>("ebRHtag")),
  doCalib_(iConfig.getUntrackedParameter<bool>("doCalib",false)),
  doBeamCounters_(iConfig.getUntrackedParameter<bool>("doBeamCounters",true)),
  calibFC2GeV_(iConfig.getUntrackedParameter<double>("calibFC2GeV",0.01))
{
  hbDigiCnt_ = 0;
  hoDigiCnt_ = 0;
  hfDigiCnt_  = 0;
  ebCnt_ = 0;
  evtCnt_ = 0;
  algo_ = new
    HcalQLPlotAnalAlgos(iConfig.getParameter<edm::ParameterSet>("HistoParameters"));
}


HcalHOTBPlotAnal::~HcalHOTBPlotAnal()
{
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
HcalHOTBPlotAnal::analyze(const edm::Event& iEvent, 
			  const edm::EventSetup& iSetup)
{
  evtCnt_++;
  // Step A/C: Get Inputs and process (repeatedly)
  edm::Handle<HcalTBTriggerData> trig;
  iEvent.getByLabel(hcalTrigLabel_,trig);
  if (!trig.isValid()) {
    edm::LogError("HcalHOTBPlotAnal::analyze") << "No Trigger Data found, "
      "skip event";
    return;
  } else {
    algo_->SetEventType(*trig);
  }
  edm::Handle<HcalTBBeamCounters> qadc;
  iEvent.getByLabel(hcalTrigLabel_,qadc);
  if (!qadc.isValid()) {
    edm::LogError("HcalHOTBPlotAnal::analyze") << 
      "No TB beam counter data found.  Skipping the event.";
    return;
  } else {
    algo_->setBeamCounters(*qadc);
    if (doBeamCounters_)
      algo_->processBeamCounters();
  }
  edm::Handle<HcalTBEventPosition> pos;
  iEvent.getByLabel(hcalTrigLabel_,pos);
  if (!pos.isValid()) {
    edm::LogWarning("HcalHOTBPlotAnal::analyse") << "no TB position data "
      "found.  Skipping for this event.";
  } else {
    algo_->setHBTableEtaPhi(pos->hbheTableEta(), pos->hbheTablePhi());
  }
  edm::Handle<HBHEDigiCollection> hbhedg;
  iEvent.getByLabel(hcalDigiLabel_,hbhedg);
  if (!hbhedg.isValid()) {
    edm::LogWarning("HcalHOTBPlotAnal::analyze") << "One of HBHE "
      "Digis/RecHits not found";
  } else {
    algo_->processDigi(*hbhedg);
    hbDigiCnt_ += hbhedg->size();
  }
  edm::Handle<HBHERecHitCollection> hbherh;  
  iEvent.getByLabel(hbheRHLabel_,hbherh);
  if (!hbherh.isValid()) {
    edm::LogWarning("HcalHOTBPlotAnal::analyze") << "One of HBHE "
      "Digis/RecHits not found";
  } else {
    algo_->processRH(*hbherh,*hbhedg);
  }

  edm::Handle<HODigiCollection> hodg;
  iEvent.getByLabel(hcalDigiLabel_,hodg);
  if (!hodg.isValid()) {
    // can't find it!
    edm::LogWarning("HcalHOTBPlotAnal::analyze") << "One of HO Digis/RecHits "
      "not found";
  } else {
    algo_->processDigi(*hodg);
    hoDigiCnt_ += hodg->size();
  }
  edm::Handle<HORecHitCollection> horh;
  iEvent.getByLabel(hoRHLabel_,horh);
  if (!horh.isValid()) {
    // can't find it!
    edm::LogWarning("HcalHOTBPlotAnal::analyze") << "One of HO Digis/RecHits "
      "not found";
  } else {
    algo_->processRH(*horh,*hodg);
  }
  
  edm::Handle<EcalRecHitCollection> ebrh;
  iEvent.getByLabel(ebRHLabel_, ebrh);
  if (!ebrh.isValid()) {
    edm::LogWarning("HcalHOTBPlotAnal::analyze") << "EB rec hits not found";
  } else {
    algo_->processRH(*ebrh);
    ebCnt_ += ebrh->size();
  }
  // edm::Handle<HFDigiCollection> hfdg;
  // iEvent.getByLabel(hcalDigiLabel_,hfdg);

  // if (!hfdg.isValid()) {
  //   // can't find it!
  //   edm::LogWarning("HcalHOTBPlotAnal::analyze") << "One of HF Digis/RecHits "
  //     "not found";
  // } else {
  //   algo_->processDigi(*hfdg);
  //   hfDigiCnt_ += hfdg->size();
  // }

  // edm::Handle<HFRecHitCollection> hfrh;
  // iEvent.getByLabel(hfRHLabel_,hfrh);
  // if (!hfrh.isValid()) {
  //   // can't find it!
  //   edm::LogWarning("HcalHOTBPlotAnal::analyze") << "One of HF Digis/RecHits not found";
  // } else {
  //   algo_->processRH(*hfrh,*hfdg);
  // }

  if (doCalib_) {
    // No rechits as of yet...
    edm::Handle<HcalCalibDigiCollection> calibdg;
    iEvent.getByLabel(hcalDigiLabel_,calibdg);
    if (!calibdg.isValid()) {
      edm::LogWarning("HcalHOTBPlotAnal::analyze") << "Hcal Calib Digis not "
	"found";
    } else {
      algo_->processDigi(*calibdg,calibFC2GeV_);
    }
  }

  algo_->fillTree();

}


// ------ method called once each job just after ending the event loop  ------
void 
HcalHOTBPlotAnal::endJob()
{
  std::cout << "avg number of HB digis: " << hbDigiCnt_/double(evtCnt_) << '\n'
	    << "avg number of HO digis: " << hoDigiCnt_/double(evtCnt_) << '\n'
	    << "avg number of HF digis: " << hfDigiCnt_/double(evtCnt_) << '\n'
	    << "avg number of EB recHits: " << ebCnt_/double(evtCnt_) << '\n';

  algo_->end();
}

//define this as a plug-in
DEFINE_FWK_MODULE(HcalHOTBPlotAnal);
