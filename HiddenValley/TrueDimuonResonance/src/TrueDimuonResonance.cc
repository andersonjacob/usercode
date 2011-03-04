// -*- C++ -*-
//
// Package:    TrueDimuonResonance
// Class:      TrueDimuonResonance
// 
/**\class TrueDimuonResonance TrueDimuonResonance.cc HiddenValley/TrueDimuonResonance/src/TrueDimuonResonance.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  "Jacob Anderson"
//         Created:  Thu Mar  3 13:49:33 CST 2011
// $Id: TrueDimuonResonance.cc,v 1.2 2011/03/03 21:43:55 andersj Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"


//
// class declaration
//

class TrueDimuonResonance : public edm::EDFilter {
   public:
      explicit TrueDimuonResonance(const edm::ParameterSet&);
      ~TrueDimuonResonance();

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void printParticleDecay(reco::Candidate const * part, 
				      std::string prefix) ;
      
      // ----------member data ---------------------------
  edm::InputTag mGenName;
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
TrueDimuonResonance::TrueDimuonResonance(const edm::ParameterSet& iConfig) :
  mGenName(iConfig.getParameter<edm::InputTag>("genParticlesName"))
{
   //now do what ever initialization is needed

}


TrueDimuonResonance::~TrueDimuonResonance()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
TrueDimuonResonance::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<reco::GenParticleCollection> genParticles;
   iEvent.getByLabel(mGenName, genParticles);
   if (!genParticles.isValid()) return false;

   reco::GenParticleCollection::const_iterator particle;
   reco::Candidate const * daughter;
   //reco::Candidate::const_iterator dauItr;

   int mucnt;
   for (particle = genParticles->begin(); particle != genParticles->end();
	++particle) {
     mucnt = 0;
     // for (dauItr = particle->begin(); dauItr != particle->end();
     // 	  ++dauItr) {
     for(unsigned int i=0; i < particle->numberOfDaughters(); ++i) {
       daughter = particle->daughter(i);
       if ( (abs(daughter->pdgId()) == 13) && 
       	    ((daughter->status() == 1) || (daughter->status() == 3)) )
       	 ++mucnt;
     }
     if (mucnt == 2) {
       printParticleDecay(&(*particle), "");
       return true;
     }
   }

   return false;
}

// ------------ method called once each job just before starting event loop  ------------
void 
TrueDimuonResonance::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TrueDimuonResonance::endJob() {
}


void TrueDimuonResonance::printParticleDecay(reco::Candidate const * part, 
						 std::string prefix) {
  std::cout << prefix << " id: " << part->pdgId()
	    << " status: " << part->status()
	    << " nDau: " << part->numberOfDaughters() << '\n';
  reco::Candidate const * daughter;
  for (unsigned int dau = 0; dau < part->numberOfDaughters(); ++dau) {
    daughter = part->daughter(dau);
    printParticleDecay(daughter,  "  " + prefix);
  }
}
//define this as a plug-in
DEFINE_FWK_MODULE(TrueDimuonResonance);
