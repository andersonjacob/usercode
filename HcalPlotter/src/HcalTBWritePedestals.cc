// -*- C++ -*-
//
// Package:    HcalTBWritePedestals
// Class:      HcalTBWritePedestals
// 
/**\class HcalTBWritePedestals HcalTBWritePedestals.cc RecoTBCalo/HcalTBWritePedestals/src/HcalTBWritePedestals.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Phillip R. Dudero
//         Created:  Tue Jan 16 21:11:37 CST 2007
// $Id: HcalTBWritePedestals.cc,v 1.1 2011/08/10 19:55:40 andersj Exp $
//
//


// system include files
#include <memory>
#include <string>
#include <vector>
#include <map>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"

#include "TBDataFormats/HcalTBObjects/interface/HcalTBTriggerData.h"
#include "usercode/HcalPlotter/src/HcalQLPlotHistoMgr.h"

#include "TF1.h"
//
// class declaration
//

class HcalTBWritePedestals : public edm::EDAnalyzer {
public:
  explicit HcalTBWritePedestals(const edm::ParameterSet&);
  ~HcalTBWritePedestals();
  
private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  void SetEventType(const HcalTBTriggerData& trigd);
  virtual void processDigi(const HBHEDigiCollection& hbhedigic,
			   const edm::EventSetup& eventSetup);
  virtual void processDigi(const HODigiCollection& hodigic,
   			   const edm::EventSetup& eventSetup);
  virtual void processDigi(const HFDigiCollection& hfdigic,
  			   const edm::EventSetup& eventSetup);

  // ----------member data ---------------------------
  edm::InputTag hcalDigiLabel_, hcalTrigLabel_;

  int hbDigiCnt_;
  int hoDigiCnt_;
  int hfDigiCnt_;
  int evtCnt_;

  edm::Service<TFileService> fs;
  HcalQLPlotHistoMgr::EventType triggerID_;
  HcalQLPlotHistoMgr *histos_;

  std::map< HcalDetId, std::vector<double>[4] > pedMap_;
};

//
// constants, enums and typedefs
//

namespace HcalTBWritePedestalImpl {
  template<class digic, class const_iter>
  inline void processDigi(const digic& digis, const edm::EventSetup& eventSetup,
			  std::map< HcalDetId, std::vector<double>[4] >& pedMap,
			  HcalQLPlotHistoMgr::EventType triggerID,
			  HcalQLPlotHistoMgr * histos) {
    edm::ESHandle<HcalDbService> conditions;
    eventSetup.get<HcalDbRecord>().get(conditions);

    const_iter it;
    char hname[100];
    TString detName;

    for (it  = digis.begin(); 
	 it != digis.end();
	 ++it) {
      if (triggerID == HcalQLPlotHistoMgr::PEDESTAL) {
	HcalDetId id (it->id());
	HcalElectronicsId eid (it->elecId());
	switch(id.subdet()){
	case HcalBarrel:
	  detName = "HB"; break;
	case HcalEndcap:
	  detName = "HE"; break;
	case HcalOuter:
	  detName = "HO"; break;
	case HcalForward:
	  detName = "HF"; break;
	default:
	  detName = "";
	}
	const HcalQIECoder* channelCoder = conditions->getHcalCoder(id);
	const HcalQIEShape* shape = conditions->getHcalShape (channelCoder);
	HcalCoderDb coder (*channelCoder, *shape);

	CaloSamples tool;
	coder.adc2fC(*it, tool);

	TH1* phist=histos->GetAHistogram(id,eid,HcalQLPlotHistoMgr::PULSE,
					 triggerID);
	if (phist) {
	  for (int bin=0; bin<tool.size(); ++bin) 
	    phist->Fill(bin*1.0,tool[bin]);	
	}

	std::map< HcalDetId, std::vector<double>[4] >::iterator peds = 
	  pedMap.find(id);
	if (peds == pedMap.end()) {
	  std::vector<double> holder[4];
	  std::map< HcalDetId, std::vector<double>[4] >::value_type p(id, holder);
	  peds = pedMap.insert(p).first;
	}
	for (int ts = 0; ts < tool.size(); ++ts) {
	  sprintf(hname, "ADC_%s_%d_%d_%d_cap_%d", detName.Data(), id.ieta(), 
		  id.iphi(), id.depth(), (*it)[ts].capid());
	  phist = histos->GetAHistogramImpl(hname, HcalQLPlotHistoMgr::ADC, 
					    triggerID);
	  if (phist)
	    phist->Fill(tool[ts]);
	  peds->second[(*it)[ts].capid()].push_back(tool[ts]);
       }
  
      }
    }
  }
}

//
// static data member definitions
//

//
// constructors and destructor
//
HcalTBWritePedestals::HcalTBWritePedestals(const edm::ParameterSet& iConfig) :
  hcalDigiLabel_(iConfig.getUntrackedParameter<edm::InputTag>("hcalDigiTag")),
  hcalTrigLabel_(iConfig.getUntrackedParameter<edm::InputTag>("hcalTrigTag"))
{
  hbDigiCnt_ = 0;
  hoDigiCnt_ = 0;
  hfDigiCnt_  = 0;
  evtCnt_ = 0;

  triggerID_=HcalQLPlotHistoMgr::UNKNOWN;
  //edm::Service<TFileService>  fs;
  histos_ = new HcalQLPlotHistoMgr(*fs, iConfig.getParameter<edm::ParameterSet>("HistoParameters"));
}


HcalTBWritePedestals::~HcalTBWritePedestals()
{
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
HcalTBWritePedestals::analyze(const edm::Event& iEvent, 
			  const edm::EventSetup& iSetup)
{
  evtCnt_++;
  // Step A/C: Get Inputs and process (repeatedly)
  edm::Handle<HcalTBTriggerData> trig;
  iEvent.getByLabel(hcalTrigLabel_,trig);
  if (!trig.isValid()) {
    edm::LogError("HcalTBWritePedestals::analyze") << "No Trigger Data found, "
      "skip event";
    return;
  } else {
    SetEventType(*trig);
  }
  edm::Handle<HBHEDigiCollection> hbhedg;
  iEvent.getByLabel(hcalDigiLabel_,hbhedg);
  if (!hbhedg.isValid()) {
    edm::LogWarning("HcalTBWritePedestals::analyze") << "One of HBHE "
      "Digis not found";
  } else {
    processDigi(*hbhedg, iSetup);
    hbDigiCnt_ += hbhedg->size();
  }

  edm::Handle<HODigiCollection> hodg;
  iEvent.getByLabel(hcalDigiLabel_,hodg);
  if (!hodg.isValid()) {
    // can't find it!
    edm::LogWarning("HcalTBWritePedestals::analyze") << "One of HO Digis "
      "not found";
  } else {
    processDigi(*hodg, iSetup);
    hoDigiCnt_ += hodg->size();
  }

  edm::Handle<HFDigiCollection> hfdg;
  iEvent.getByLabel(hcalDigiLabel_,hfdg);
  if (!hfdg.isValid()) {
    // can't find it!
    edm::LogWarning("HcalTBWritePedestals::analyze") << "One of HF Digis "
      "not found";
  } else {
    processDigi(*hfdg, iSetup);
    hfDigiCnt_ += hfdg->size();
  }

}

void HcalTBWritePedestals::SetEventType(const HcalTBTriggerData& trigd)
{
  if( trigd.wasInSpillPedestalTrigger()  ||
      trigd.wasOutSpillPedestalTrigger() ||
      trigd.wasSpillIgnorantPedestalTrigger() )
                                triggerID_=HcalQLPlotHistoMgr::PEDESTAL;
  if( trigd.wasLEDTrigger() )   triggerID_=HcalQLPlotHistoMgr::LED;
  if( trigd.wasLaserTrigger() ) triggerID_=HcalQLPlotHistoMgr::LASER;
  if( trigd.wasBeamTrigger() )  triggerID_=HcalQLPlotHistoMgr::BEAM;

  if( triggerID_ == HcalQLPlotHistoMgr::UNKNOWN) {
    edm::LogError("HcalQLPlotAnalAlgos::begin") <<
      "Trigger Type unrecognized, aborting";
    std::exception e;
    throw e;
  }
}

void HcalTBWritePedestals::processDigi(const HBHEDigiCollection& hbhedigic,
				       const edm::EventSetup& eventSetup) {
  HcalTBWritePedestalImpl::processDigi<HBHEDigiCollection, 
    HBHEDigiCollection::const_iterator>(hbhedigic,eventSetup,
					pedMap_, triggerID_,
					histos_);
}

void HcalTBWritePedestals::processDigi(const HODigiCollection& hodigic,
				       const edm::EventSetup& eventSetup) {
  HcalTBWritePedestalImpl::processDigi<HODigiCollection, 
    HODigiCollection::const_iterator>(hodigic,eventSetup,
				      pedMap_, triggerID_,
				      histos_);
}

void HcalTBWritePedestals::processDigi(const HFDigiCollection& hfdigic,
				       const edm::EventSetup& eventSetup) {
  HcalTBWritePedestalImpl::processDigi<HFDigiCollection, 
    HFDigiCollection::const_iterator>(hfdigic,eventSetup,
				      pedMap_, triggerID_,
				      histos_);
}

// ------ method called once each job just after ending the event loop  ------
void 
HcalTBWritePedestals::endJob()
{
  std::cout << "avg number of HB digis: " << hbDigiCnt_/double(evtCnt_) << '\n'
	    << "avg number of HO digis: " << hoDigiCnt_/double(evtCnt_) << '\n'
	    << "avg number of HF digis: " << hfDigiCnt_/double(evtCnt_) << '\n';

  std::ofstream pedFile("new_pedestals.txt", std::ios::out | std::ios::trunc);
  std::ofstream gainFile("new_gains.txt", std::ios::out | std::ios::trunc);

  pedFile << "#U ADC\n"
	  << "# eta   phi   dep   det   cap1   cap2   cap3   cap4   HcalDetId\n";
  gainFile << "#U ADC\n"
	   << "# eta   phi   dep   det   cap1   cap2   cap3   cap4   HcalDetId\n";

  std::map< HcalDetId, std::vector<double>[4] >::iterator ped;
  TString idStr = "  ";
  TH1 * phist = 0;
  char hname[100];
  double avgPed = 0.;

  for (ped = pedMap_.begin(); ped != pedMap_.end(); ++ped) {
    
    pedFile << ped->first.ieta() << "   "
	    << ped->first.iphi() << "   "
	    << ped->first.depth() << "   ";
    gainFile << ped->first.ieta() << "   "
	     << ped->first.iphi() << "   "
	     << ped->first.depth() << "   ";
    switch (ped->first.subdet()) {
    case HcalBarrel: 
      pedFile << "HB   "; 
      gainFile << "HB   "; 
      idStr = "HB";
      break;
    case HcalEndcap: 
      pedFile << "HE   ";
      gainFile << "HE   ";
      idStr = "HE";
      break;
    case HcalForward: 
      pedFile << "HF   "; 
      gainFile << "HF   "; 
      idStr = "HF";
      break;
    case HcalOuter: 
      pedFile << "HO   "; 
      gainFile << "HO   "; 
      idStr = "HO";
      break;
    default: 
      pedFile << "   ";
      gainFile << "   ";
      idStr = "";
    }
    pedFile << std::setprecision(4);
    gainFile << std::setprecision(4);
    for (int cap = 0; cap < 4; ++cap) {
      std::vector<double> & peds = ped->second[cap];
      double pedSum = 0;
      int pedCnt = 0;
      std::vector<double>::iterator p;

      sprintf(hname, "ADC_%s_%d_%d_%d_cap_%d", idStr.Data(), ped->first.ieta(), 
	      ped->first.iphi(), ped->first.depth(), cap);
      phist = histos_->GetAHistogramImpl(hname, HcalQLPlotHistoMgr::ADC, 
					 HcalQLPlotHistoMgr::PEDESTAL, true);

      double minimum(-0.5);
      double maximum(1e5);
      int maxBin(0);

      if (phist) {
	maxBin = phist->GetMaximumBin();
	minimum = phist->GetBinLowEdge(maxBin-1);
	maximum = phist->GetBinLowEdge(maxBin+2);
      }

      for (p = peds.begin(); p != peds.end(); ++p) {
	if ((*p > minimum) && (*p < maximum)) {
	  pedSum += *p;
	  ++pedCnt;
	}
      }

      avgPed = pedSum/pedCnt;
      if (phist) {
      	// int maxBin = phist->GetMaximumBin();
      	TF1 * ped = new TF1("ped", "gaus", phist->GetBinLowEdge(1),
      			    phist->GetBinLowEdge(maxBin+2));
      	phist->Fit(ped, "RQ");

      	// int minBin = (phist->GetMaximumBin()-2 > 0) 
      	//   ? phist->GetMaximumBin()-2 : 1;
      	// double sum = 0.;
      	// double sumw = 0.;
      	// for (int bin = minBin; bin < minBin+5; ++bin) {
      	//   sum += phist->GetBinCenter(bin)*phist->GetBinContent(bin);
      	//   sumw += phist->GetBinContent(bin);
      	// }
      	// if (sumw > 0.)
      	//   avgPed = sum/sumw;
      	avgPed = ped->GetParameter(1);
      }
      pedFile << avgPed << "   ";
      gainFile << 1.0 << "   ";
    }
    pedFile << std::hex << ped->first.rawId() << std::dec << '\n';
    gainFile << std::hex << ped->first.rawId() << std::dec << '\n';
  }

  pedFile.close();
  gainFile.close();

}

//define this as a plug-in
DEFINE_FWK_MODULE(HcalTBWritePedestals);
