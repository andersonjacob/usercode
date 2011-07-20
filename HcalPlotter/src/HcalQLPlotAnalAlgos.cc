// -*- C++ -*-
//
// Package:    HcalQLPlotAnal
// Class:      HcalQLPlotAnal
// 
/**\class HcalQLPlotAnal HcalQLPlotAnal.cc MyEDProducts/HcalPlotter/src/HcalQLPlotAnal.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Phillip R. Dudero
//         Created:  Tue Jan 16 21:11:37 CST 2007
// $Id: HcalQLPlotAnalAlgos.cc,v 1.5 2010/06/21 07:47:58 mansj Exp $
//
//


// system include files
#include <memory>
#include <math.h>
#include <fstream>

// user include files
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "andersj/HcalPlotter/src/HcalQLPlotAnalAlgos.h"
#include "DataFormats/HcalDetId/interface/HcalElectronicsId.h"
#include "DataFormats/HcalRecHit/interface/HcalCalibRecHit.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"

#include "TH1.h"
#include "TH2F.h"
#include "TH1F.h"
//
// constants, enums and typedefs
//

//
// static data member definitions
//
int const HcalQLPlotAnalAlgos::maxDim;
//
// constructors and destructor
//
HcalQLPlotAnalAlgos::HcalQLPlotAnalAlgos(const edm::ParameterSet& histoParams) :
  qadc(0), 
  dataTree(0), HBLego(0), HOLego(0), EBLegoFine(0), EBieta(0),
  EBiphi(0),
  VMBadc_(0.), VMFadc_(0), VMadc_(0.),
  S1adc_(0.), S2adc_(0.), S3adc_(0.), S4adc_(0.),
  BH1adc_(0.), BH2adc_(0.), BH3adc_(0.), BH4adc_(0.),
  maxEtaHO(0), maxPhiHO(0), maxEtaHB(0), maxPhiHB(0),
  HOE1(0.), HBE1(0.)
{
  char name[100];
  triggerID_=HcalQLPlotHistoMgr::UNKNOWN;
  //edm::Service<TFileService>  fs;
  histos_ = new HcalQLPlotHistoMgr(*fs,histoParams);
  dataTree = fs->make<TTree>("dataTree", "dataTree");
  dataTree->Branch("triggerID", &triggerID_, "triggerID/I");
  dataTree->Branch("VMBadc", &VMBadc_, "VMBadc/D");
  dataTree->Branch("VMFadc", &VMFadc_, "VMFadc/D");
  dataTree->Branch("VMadc", &VMadc_, "VMadc/D");
  dataTree->Branch("S1adc", &S1adc_, "S1adc/D");
  dataTree->Branch("S2adc", &S2adc_, "S2adc/D");
  dataTree->Branch("S3adc", &S3adc_, "S3adc/D");
  dataTree->Branch("S4adc", &S4adc_, "S4adc/D");
  dataTree->Branch("BH1adc", &BH1adc_, "BH1adc/D");
  dataTree->Branch("BH2adc", &BH2adc_, "BH2adc/D");
  dataTree->Branch("BH3adc", &BH3adc_, "BH3adc/D");
  dataTree->Branch("BH4adc", &BH4adc_, "BH4adc/D");
  dataTree->Branch("maxEtaHO", &maxEtaHO, "maxEtaHO/I");
  dataTree->Branch("maxPhiHO", &maxPhiHO, "maxPhiHO/I");
  dataTree->Branch("maxEtaHB", &maxEtaHB, "maxEtaHB/I");
  dataTree->Branch("maxPhiHB", &maxPhiHB, "maxPhiHB/I");
  dataTree->Branch("HOE1", &HOE1, "HOE1/D");
  dataTree->Branch("HBE1", &HBE1, "HBE1/D");
  sprintf(name, "HBE[%d][%d][%d]/D", maxDim, maxDim, maxDepth);
  dataTree->Branch("HBE", HBE, name);
  sprintf(name, "HOE[%d][%d]/D", maxDim, maxDim);
  dataTree->Branch("HOE", HOE, name);

  EBLego = fs->make<TH2F>("eblego", "EBLego", 33, -16.5, 16.5, 72, 0.5, 72.5);
  HBLego = fs->make<TH2F>("hblego", "HBLego", 33, -16.5, 16.5, 72, 0.5, 72.5);
  HOLego = fs->make<TH2F>("holego", "HOLego", 33, -16.5, 16.5, 72, 0.5, 72.5);
  EBLegoFine = fs->make<TH2F>("eblegofine", "EBLegoFine", 100, -100.5, -0.5,
			      20, 0.5, 20.5);

  EBieta = fs->make<TH1F>("ebieta", "ebieta", 100, -100.5, -0.5);
  EBiphi = fs->make<TH1F>("ebiphi", "ebiphi", 20, 0.5, 20.5);
}


//
// member functions
//

void HcalQLPlotAnalAlgos::end(void)
{
  std::ofstream pedFile("new_pedestals.txt", std::ios::out | std::ios::trunc);
  std::ofstream gainFile("new_gains.txt", std::ios::out | std::ios::trunc);

  pedFile << "#U ADC\n"
	  << "# eta   phi   dep   det   cap1   cap2   cap3   cap4   HcalDetId\n";
  gainFile << "#U ADC\n"
	   << "# eta   phi   dep   det   cap1   cap2   cap3   cap4   HcalDetId\n";

  std::map< HcalDetId, std::vector<int> >::iterator ped;

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
      break;
    case HcalEndcap: 
      pedFile << "HE   ";
      gainFile << "HE   ";
      break;
    case HcalForward: 
      pedFile << "HF   "; 
      gainFile << "HF   "; 
      break;
    case HcalOuter: 
      pedFile << "HO   "; 
      gainFile << "HO   "; 
      break;
    default: 
      pedFile << "   ";
      gainFile << "   ";
    }
    pedFile << std::setprecision(4);
    gainFile << std::setprecision(4);
    for (int cap = 0; cap < 4; ++cap) {
      pedFile << ped->second[cap]/double(ped->second[cap+4]) << "   ";
      gainFile << 1.0 << "   ";
    }
    pedFile << std::hex << ped->first.rawId() << std::dec << '\n';
    gainFile << std::hex << ped->first.rawId() << std::dec << '\n';
  }

  pedFile.close();
  gainFile.close();
}

void HcalQLPlotAnalAlgos::SetEventType(const HcalTBTriggerData& trigd)
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

void HcalQLPlotAnalAlgos::processBeamCounters() 
{
  VMBadc_ = 0;
  VMFadc_ = 0;
  VMadc_ = 0;
  S1adc_ = 0;
  S2adc_ = 0;
  S3adc_ = 0;
  S4adc_ = 0;
  BH1adc_ = 0;
  BH2adc_ = 0;
  BH3adc_ = 0;
  BH4adc_ = 0;
  if (qadc) {
    VMBadc_ = qadc->VMBadc();
    VMFadc_ = qadc->VMFadc(); 
    VMadc_ = qadc->VMadc();
    S1adc_ = qadc->S1adc();
    S2adc_ = qadc->S2adc();
    S3adc_ = qadc->S3adc();
    S4adc_ = qadc->S4adc();
    BH1adc_ = qadc->BH1adc();
    BH2adc_ = qadc->BH2adc();
    BH3adc_ = qadc->BH3adc();
    BH4adc_ = qadc->BH4adc();
  }
  //dataTree->Fill();
}

void HcalQLPlotAnalAlgos::processRH(const HBHERecHitCollection& hbherhc,
				    const HBHEDigiCollection& hbhedgc)
{
  HBHERecHitCollection::const_iterator it;

  std::map<int, double> towerE;

  for (int e = 0; e < maxDim; ++e)
    for (int p = 0; p < maxDim; ++p)
      for (int d = 0; d < maxDepth; ++d)
      HBE[e][p][d] = 0.;

  for (it  = hbherhc.begin(); 
       it != hbherhc.end();
       ++it) {
    HcalDetId id (it->id());
    HcalElectronicsId eid;
    HBHEDigiCollection::const_iterator dit = hbhedgc.find(id);
    if (dit != hbhedgc.end())
      eid = dit->elecId();
    else {
      edm::LogWarning("HcalQLPlotAnalAlgos::processRH") <<
	"No digi found for id" << id;
      continue;
    }
    int tmpPhi = id.iphi()%10;
    int tmpDepth = id.iphi()/10;
    if (triggerID_ == HcalQLPlotHistoMgr::BEAM) {
      HBLego->Fill(id.ieta(),tmpPhi,it->energy());
      int tmpIndex = id.ieta()*100+tmpPhi;
      std::map<int, double>::iterator tower = towerE.find(tmpIndex);
      if (tower == towerE.end()) {
	tower = towerE.insert(std::pair<int, double>(tmpIndex, 0.)).first;
      }
      tower->second += it->energy();
    }
    if ( (id.ieta() < maxDim) && (id.ieta() >= 0) &&
	 (tmpPhi < maxDim) && (tmpPhi >= 0) &&
	 (tmpDepth < maxDepth) && (tmpDepth >= 0) ) {
      HBE[id.ieta()][tmpPhi][tmpDepth] += it->energy();
    }
    // TH1* ehist=histos_->GetAHistogram(id,eid,HcalQLPlotHistoMgr::ENERGY,triggerID_);
    // if (ehist){
    //   ehist->Fill(it->energy());
    // }

    // TH1* thist=histos_->GetAHistogram(id,eid,HcalQLPlotHistoMgr::TIME,triggerID_);
    // if (thist){
    //   thist->Fill(it->time());
    // }
  }

  maxEtaHB = 0;
  maxPhiHB = 0;
  HBE1 = 0.;
  std::map<int, double>::iterator tower;
  for (tower = towerE.begin(); tower != towerE.end(); ++tower) {
    if (tower->second > HBE1) {
      maxEtaHB = tower->first/100;
      maxPhiHB = tower->first%100;
      HBE1 = tower->second;
    }
  }
}

void HcalQLPlotAnalAlgos::processRH(const HORecHitCollection& horhc,
				    const HODigiCollection& hodgc)
{
  HORecHitCollection::const_iterator it;

  for (int e = 0; e < maxDim; ++e)
    for (int p = 0; p < maxDim; ++p)
      HOE[e][p] = 0.;

  maxEtaHO = 0;
  maxPhiHO = 0;
  HOE1 = 0.;

  for (it  = horhc.begin(); 
       it != horhc.end();
       ++it) {
    HcalDetId id (it->id());
    HcalElectronicsId eid;
    HODigiCollection::const_iterator dit = hodgc.find(id);
    if (dit != hodgc.end())
      eid = dit->elecId();
    else {
      edm::LogWarning("HcalQLPlotAnalAlgos::processRH") <<
	"No digi found for id" << id;
      continue;
    }
    if (triggerID_ == HcalQLPlotHistoMgr::BEAM) {
      HOLego->Fill(id.ieta(), id.iphi(), it->energy());
      if (it->energy() > HOE1) {
	maxEtaHO = id.ieta();
	maxPhiHO = id.iphi();
	HOE1 = it->energy();
      }
    }
    if ( (id.ieta() < maxDim) && (id.ieta() >= 0) &&
	 (id.iphi() < maxDim) && (id.iphi() >= 0) ) {
      HOE[id.ieta()][id.iphi()] += it->energy();
    }

//     TH1* thist=histos_->GetAHistogram(id,eid,HcalQLPlotHistoMgr::TIME,triggerID_);
//     if (thist){
//       thist->Fill(it->time());
//     }
  }
}

void HcalQLPlotAnalAlgos::processRH(const HFRecHitCollection& hfrhc,
				    const HFDigiCollection& hfdgc)
{
  HFRecHitCollection::const_iterator it;

  for (it  = hfrhc.begin(); 
       it != hfrhc.end();
       ++it) {
    HcalDetId id (it->id());
    HcalElectronicsId eid;
    HFDigiCollection::const_iterator dit = hfdgc.find(id);
    if (dit != hfdgc.end())
      eid = dit->elecId();
    else {
      edm::LogWarning("HcalQLPlotAnalAlgos::processRH") <<
	"No digi found for id" << id;
      continue;
    }

//     TH1* ehist=histos_->GetAHistogram(id,eid,HcalQLPlotHistoMgr::ENERGY,triggerID_);
//     if (ehist){
//       ehist->Fill(it->energy());
//     }

//     TH1* thist=histos_->GetAHistogram(id,eid,HcalQLPlotHistoMgr::TIME,triggerID_);
//     if (thist){
//       thist->Fill(it->time());
//     }
  }
}

void HcalQLPlotAnalAlgos::processRH(const EcalRecHitCollection& ebrhc) 
{
  EcalRecHitCollection::const_iterator it;

  if (triggerID_ == HcalQLPlotHistoMgr::BEAM)
    for (it = ebrhc.begin(); it != ebrhc.end(); ++it) {
      EBDetId id(it->id());
      EBLego->Fill(id.tower_ieta(), id.tower_iphi(), it->energy());
      EBieta->Fill(id.ieta(), it->energy());
      EBiphi->Fill(id.iphi(), it->energy());
      EBLegoFine->Fill(id.ieta(), id.iphi(), it->energy());
    }
}

void HcalQLPlotAnalAlgos::processDigi(const HBHEDigiCollection& hbhedigic)
{
  HBHEDigiCollection::const_iterator it;

  for (it  = hbhedigic.begin(); 
       it != hbhedigic.end();
       ++it) {
    HcalDetId id (it->id());
    HcalElectronicsId eid (it->elecId());

    TH1* phist=histos_->GetAHistogram(id,eid,HcalQLPlotHistoMgr::PULSE,triggerID_);
    if (phist){
      for (int bin=0; bin<it->size(); bin++)
	phist->Fill(bin*1.0,(*it)[bin].nominal_fC());
    }

    if (triggerID_ == HcalQLPlotHistoMgr::PEDESTAL) {
      phist=histos_->GetAHistogram(id,eid,HcalQLPlotHistoMgr::ADC,triggerID_);
      if (phist){
	for (int bin=0; bin<it->size(); bin++)
	  phist->Fill((*it)[bin].adc());
      }
      std::map< HcalDetId, std::vector<int> >::iterator peds = pedMap_.find(id);
      if (peds == pedMap_.end()) {
	std::vector<int> holder;
	for (int i=0; i<8; ++i) holder.push_back(0);
	peds = pedMap_.insert(std::pair< HcalDetId, std::vector<int> >(id, holder)).first;
      }
      for (int ts = 0; ts < it->size(); ++ts) {
	peds->second[(*it)[ts].capid()] += (*it)[ts].adc();
	peds->second[(*it)[ts].capid()+4] += 1;
      }
    }
  }
}

void HcalQLPlotAnalAlgos::processDigi(const HODigiCollection& hodigic)
{
  HODigiCollection::const_iterator it;

//   for (int i = 0; i<10; ++i) 
//     HO_3_3_digi[i] = 0.;
  for (it  = hodigic.begin(); 
       it != hodigic.end();
       ++it) {
    HcalDetId id (it->id());
    HcalElectronicsId eid (it->elecId());

    TH1* phist=histos_->GetAHistogram(id,eid,HcalQLPlotHistoMgr::PULSE,triggerID_);
    if (phist){
      for (int bin=0; bin<it->size(); bin++) {
	phist->Fill(bin*1.0,(*it)[bin].nominal_fC());
// 	if ( (bin < 10) && (id.ieta() == 3) && (id.iphi() == 3) ) 
// 	  HO_3_3_digi[bin] = (*it)[bin].nominal_fC();
      }
    }

    if (triggerID_ == HcalQLPlotHistoMgr::PEDESTAL) {
      phist=histos_->GetAHistogram(id,eid,HcalQLPlotHistoMgr::ADC,triggerID_);
      if (phist){
	for (int bin=0; bin<it->size(); bin++)
	  phist->Fill((*it)[bin].adc());
      }
      std::map< HcalDetId, std::vector<int> >::iterator peds = pedMap_.find(id);
      if (peds == pedMap_.end()) {
	std::vector<int> holder;
	for (int i=0; i<8; ++i) holder.push_back(0);
	peds = pedMap_.insert(std::pair< HcalDetId, std::vector<int> >(id, holder)).first;
      }
      for (int ts = 0; ts < it->size(); ++ts) {
	peds->second[(*it)[ts].capid()] += (*it)[ts].adc();
	peds->second[(*it)[ts].capid()+4] += 1;
      }
    }
  }
}

void HcalQLPlotAnalAlgos::processDigi(const HFDigiCollection& hfdigic)
{
  HFDigiCollection::const_iterator it;

  for (it  = hfdigic.begin(); 
       it != hfdigic.end();
       ++it) {
    HcalDetId id (it->id());
    HcalElectronicsId eid (it->elecId());

    TH1* phist=histos_->GetAHistogram(id,eid,HcalQLPlotHistoMgr::PULSE,triggerID_);
    if (phist){
      for (int bin=0; bin<it->size(); bin++)
	phist->Fill(bin*1.0,(*it)[bin].nominal_fC());
    }

    if (triggerID_ == HcalQLPlotHistoMgr::PEDESTAL) {
      phist=histos_->GetAHistogram(id,eid,HcalQLPlotHistoMgr::ADC,triggerID_);
      if (phist){
	for (int bin=0; bin<it->size(); bin++)
	  phist->Fill((*it)[bin].adc());
      }
      std::map< HcalDetId, std::vector<int> >::iterator peds = pedMap_.find(id);
      if (peds == pedMap_.end()) {
	std::vector<int> holder;
	for (int i=0; i<8; ++i) holder.push_back(0);
	peds = pedMap_.insert(std::pair< HcalDetId, std::vector<int> >(id, holder)).first;
      }
      for (int ts = 0; ts < it->size(); ++ts) {
	peds->second[(*it)[ts].capid()] += (*it)[ts].adc();
	peds->second[(*it)[ts].capid()+4] += 1;
      }
    }
  }
}

HcalCalibRecHit HcalQLPlotAnalAlgos::recoCalib(const HcalCalibDataFrame& cdigi,
					       double calibFC2GeV)
{
  double nominal_ped = (cdigi[0].nominal_fC() + cdigi[1].nominal_fC())/2.0;

  double totamp = 0.0;
  double maxA = -1e99;
  int    maxI = -1;
  for (int i=0; i<cdigi.size(); i++) {
    double ampl = (cdigi[i].nominal_fC()-nominal_ped)*calibFC2GeV;
    totamp += ampl;

    if (ampl > maxA) {
      maxA = ampl;
      maxI = i;
    }
  }

  maxA=fabs(maxA);
  float t0 = (maxI > 0) ? (fabs((cdigi[maxI-1].nominal_fC()-nominal_ped))*calibFC2GeV):0.0;
  float t2 = fabs((cdigi[maxI+1].nominal_fC()-nominal_ped)*calibFC2GeV);    
  float wpksamp = (maxA + 2.0*t2) / (t0 + maxA + t2);
  float time = (maxI - cdigi.presamples() + wpksamp)*25.0;

  return HcalCalibRecHit(cdigi.id(),totamp,time);    
}

void HcalQLPlotAnalAlgos::processDigi(const HcalCalibDigiCollection& calibdigic,
				      double calibFC2GeV)
{
  HcalCalibDigiCollection::const_iterator it;

  for (it  = calibdigic.begin(); 
       it != calibdigic.end();
       ++it) {
    HcalCalibDetId     id (it->id());
    HcalElectronicsId eid (it->elecId());

    TH1* phist=histos_->GetAHistogram(id,eid,HcalQLPlotHistoMgr::PULSE,triggerID_);
    if (phist){
      for (int bin=0; bin<it->size(); bin++)
	phist->Fill(bin*1.0,(*it)[bin].nominal_fC());
    }

    // HACK-reco the calib digi into a rechit:
    //
    HcalCalibRecHit rh = recoCalib(*it, calibFC2GeV);

    TH1* ehist=histos_->GetAHistogram(id,eid,HcalQLPlotHistoMgr::ENERGY,triggerID_);
    if (ehist){
      ehist->Fill(rh.amplitude());
    }

    TH1* thist=histos_->GetAHistogram(id,eid,HcalQLPlotHistoMgr::TIME,triggerID_);
    if (thist){
      thist->Fill(rh.time());
    }
  }
}

