// -*- C++ -*-
//
// Package:    PedTextWriter
// Class:      PedTextWriter
// 
/**\class PedTextWriter PedTextWriter.cc usercode/PedTextWriter/plugins/PedTextWriter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  "Jacob Anderson"
//         Created:  Mon, 22 Jul 2013 18:43:43 GMT
// $Id$
//
//


// system include files
#include <memory>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"

//
// class declaration
//

class PedTextWriter : public edm::EDAnalyzer {
public:
  typedef std::map< HcalDetId, std::vector<float> > mapPerCap;
  explicit PedTextWriter(const edm::ParameterSet&);
  ~PedTextWriter();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  void writeFile(char const * fname, mapPerCap& theMap);
  // ----------member data ---------------------------

  edm::InputTag hcalDigiLabel_;

  mapPerCap pedMap_, pedWidthMap_, gainMap_;

};

//
// constants, enums and typedefs
//

namespace PedTextWriterImpl {
  template<class digic, class const_iter>
  void processDigi(digic const& digis, HcalDbService const* conditions,
		   PedTextWriter::mapPerCap& pedMap,
		   PedTextWriter::mapPerCap& pedWidthMap,
		   PedTextWriter::mapPerCap& gainMap) {

    const_iter it;
    for (it = digis.begin(); it != digis.end(); ++it) {
      HcalDetId id(it->id());

      HcalPedestal const * pedestal = conditions->getPedestal(id);
      HcalPedestalWidth const * pedestalWidth = conditions->getPedestalWidth(id);
      HcalGain const * gain = conditions->getGain(id);

      for (int cap = 0; cap < 4; ++cap) {
	pedMap[id].push_back(pedestal->getValue(cap));
	pedWidthMap[id].push_back(pedestalWidth->getWidth(cap));
	gainMap[id].push_back(gain->getValue(cap));
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
PedTextWriter::PedTextWriter(const edm::ParameterSet& iConfig) :
  hcalDigiLabel_(iConfig.getUntrackedParameter<edm::InputTag>("hcalDigiLabel"))
{
   //now do what ever initialization is needed

}


PedTextWriter::~PedTextWriter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
PedTextWriter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   ESHandle<HcalDbService> conditions;
   iSetup.get<HcalDbRecord>().get(conditions);

   Handle<HBHEDigiCollection> hbhedigis;
   iEvent.getByLabel(hcalDigiLabel_, hbhedigis);
   if (!hbhedigis.isValid()) {
   } else {
     PedTextWriterImpl::processDigi<HBHEDigiCollection, 
				    HBHEDigiCollection::const_iterator>
       (*hbhedigis, conditions.product(), pedMap_, pedWidthMap_, gainMap_);
   }

   Handle<HODigiCollection> hodigis;
   iEvent.getByLabel(hcalDigiLabel_, hodigis);
   if (!hodigis.isValid()) {
   } else {
     PedTextWriterImpl::processDigi<HODigiCollection, 
				    HODigiCollection::const_iterator>
       (*hodigis, conditions.product(), pedMap_, pedWidthMap_, gainMap_);
   }

   Handle<HFDigiCollection> hfdigis;
   iEvent.getByLabel(hcalDigiLabel_, hfdigis);
   if (!hfdigis.isValid()) {
   } else {
     PedTextWriterImpl::processDigi<HFDigiCollection, 
				    HFDigiCollection::const_iterator>
       (*hfdigis, conditions.product(), pedMap_, pedWidthMap_, gainMap_);
   }

   Handle<ZDCDigiCollection> zdcdigis;
   iEvent.getByLabel(hcalDigiLabel_, zdcdigis);
   if (!zdcdigis.isValid()) {
   } else {
     PedTextWriterImpl::processDigi<ZDCDigiCollection, 
				    ZDCDigiCollection::const_iterator>
       (*zdcdigis, conditions.product(), pedMap_, pedWidthMap_, gainMap_);
   }

}


// ------------ method called once each job just before starting event loop  ------------
void 
PedTextWriter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PedTextWriter::endJob() 
{
  std::cout << "digi size: " << pedMap_.size() << '\n';

  writeFile("pedestals.txt", pedMap_);
  writeFile("pedestalWidths.txt", pedWidthMap_);
  writeFile("gains.txt", gainMap_);
}

void PedTextWriter::writeFile(char const * fname, mapPerCap& theMap) {

  std::ofstream outf(fname, std::ios::out | std::ios::trunc);

  outf << "#U ADC\n"
       << "# eta   phi   dep   det   cap1   cap2   cap3   cap4   HcalDetId\n";

  mapPerCap::const_iterator it;

  for (it = theMap.begin(); it != theMap.end(); ++it) {
    outf << it->first.ieta() << "   "
	 << it->first.iphi() << "   "
	 << it->first.depth() << "   ";
    switch (it->first.subdet()) {
    case HcalBarrel:
      outf << "HB   ";
      break;
    case HcalEndcap:
      outf << "HE   ";
      break;
    case HcalForward:
      outf << "HF   ";
      break;
    case HcalOuter:
      outf << "HO   ";
      break;
    default :
      std::cout << "unknown sub detector " << it->first << "\n";
    }
    outf << std::setprecision(4);
    for (int cap = 0; cap<4; ++cap) {
      outf << it->second[cap] << "   ";
    }
    outf << std::hex << it->first.rawId() << std::dec << "\n";
  }

  outf.close();
}

// ------------ method called when starting to processes a run  ------------
/*
void 
PedTextWriter::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
PedTextWriter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
PedTextWriter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
PedTextWriter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PedTextWriter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PedTextWriter);
