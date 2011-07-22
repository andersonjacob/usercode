// -*- C++ -*-

#ifndef PlotAllAnalAlgos_included
#define PlotAllAnalAlgos_included 1

// user include files
#include <map>
#include <vector>

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalCalibDataFrame.h"
#include "TBDataFormats/HcalTBObjects/interface/HcalTBTriggerData.h"
#include "TBDataFormats/HcalTBObjects/interface/HcalTBBeamCounters.h"
#include "andersj/HcalPlotter/src/HcalQLPlotHistoMgr.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "TTree.h"

//
// class declaration
//

class TH2F;
class TH1F;

class HcalQLPlotAnalAlgos {
public:
  HcalQLPlotAnalAlgos(
		   const edm::ParameterSet& histoParams);

  void end(void);
  void SetEventType(const HcalTBTriggerData& trigd) ;
  void setBeamCounters( HcalTBBeamCounters const& counters ) {
    qadc = &counters;
  }
  void setHBTableEtaPhi( double eta, double phi);
  void processBeamCounters();
  void processRH(const HBHERecHitCollection& hbherhc,
		 const HBHEDigiCollection& hbhedgc);
  void processRH(const HORecHitCollection& horhc,
		 const HODigiCollection& hodgc);
  void processRH(const HFRecHitCollection& hfrhc,
		 const HFDigiCollection& hfdgc);
  void processRH(const EcalRecHitCollection & ebrhc);
  void processDigi(const HBHEDigiCollection& hbhedigic);
  void processDigi(const HODigiCollection& hodigic);
  void processDigi(const HFDigiCollection& hfdigic);
  void processDigi(const HcalCalibDigiCollection& calibdigic,double calibFC2GeV);
  void fillTree() { dataTree->Fill(); }

private:
  HcalCalibRecHit recoCalib(const HcalCalibDataFrame& cdigi,
			    double calibFC2GeV);

  // ----------member data ---------------------------

  HcalQLPlotHistoMgr::EventType triggerID_;
  HcalQLPlotHistoMgr *histos_;

  edm::Service<TFileService> fs;

  HcalTBBeamCounters const * qadc;
  static int const maxDim = 11;
  static int const maxDepth = 5;
  static int const ebMaxPhi = 21;
  static int const ebMaxEta = 86;

  TTree * dataTree;
  TH2F * HBLego, * EBLego, * HOLego, * EBLegoFine;
  TH1F * EBieta, * EBiphi;
  TH1F * HBEnergy;
  double VMBadc_, VMFadc_, VMadc_;
  double S1adc_, S2adc_, S3adc_, S4adc_;
  double BH1adc_, BH2adc_, BH3adc_, BH4adc_;
  double HBTableEta, HBTablePhi;
  int ietaTable, iphiTable;
  int maxEtaHO, maxPhiHO, maxEtaHB, maxPhiHB, maxEtaEB, maxPhiEB;
  int NHBdigis, NHOdigis, NEBrecHits;
  double HOE1, HBE1, EBE1;
  double HOE[maxDim][maxDim];
  double HBE[maxDim][maxDim][maxDepth];
  double EBE[ebMaxEta][ebMaxPhi];
  double HBE9, HOE9, EBE25, EBE81, EBE225;
  std::map< HcalDetId, std::vector< int > > pedMap_;
};

#endif // HcalQLPlotAnalAlgos_included
