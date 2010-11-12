// -*- C++ -*-
#include <vector>
#include <list>
#include <stdint.h>
#include "FWCore/Framework/interface/EventSetup.h"

class TMultiGraph;

class HOAcceptance {
 public:
  HOAcceptance();
  virtual ~HOAcceptance() { };
  bool isChannelDead(uint32_t id);
  bool isChannelSiPM(uint32_t id);
  bool inGeomAccept(double eta, double phi, double delta_eta = 0.,
			   double delta_phi = 0.);
  bool inNotDeadGeom(double eta, double phi, double delta_eta = 0.,
			    double delta_phi = 0.);
  bool inSiPMGeom(double eta, double phi, double delta_eta = 0.,
			 double delta_phi = 0.);
  void initIds(edm::EventSetup const& eSetup);
  bool Inited() { return inited; }
  TMultiGraph * graphDeadRegions() { return graphRegions(deadRegions); }
  TMultiGraph * graphSiPMRegions() { return graphRegions(SiPMRegions); }

 private:

  struct deadRegion {
    deadRegion( double eMin = 0., double eMax = 0., 
		double pMin = 0., double pMax = 0. ) :
      etaMin(eMin), etaMax(eMax), phiMin(pMin), phiMax(pMax) { }
    deadRegion( deadRegion const& other ) :
      etaMin(other.etaMin), etaMax(other.etaMax), 
      phiMin(other.phiMin), phiMax(other.phiMax) { }
    double etaMin;
    double etaMax;
    double phiMin;
    double phiMax;
    bool operator== ( deadRegion const& other) {
      return ((other.etaMin==etaMin) && (other.etaMax==etaMax) &&
	      (other.phiMin==phiMin) && (other.phiMax==phiMax));
    }
  };

  struct deadIdRegion {
    deadIdRegion( int eMin = 0, int eMax = 0, int pMin = 0, int pMax = 0 ) :
      etaMin(eMin), etaMax(eMax), phiMin(pMin), phiMax(pMax) { }
    deadIdRegion( deadIdRegion const& other ) :
      etaMin(other.etaMin), etaMax(other.etaMax), 
      phiMin(other.phiMin), phiMax(other.phiMax) { }
    int etaMin;
    int etaMax;
    int phiMin;
    int phiMax;
    bool operator== ( deadIdRegion const& other ) { 
      return ((other.etaMin==etaMin) && (other.etaMax==etaMax) &&
	      (other.phiMin==phiMin) && (other.phiMax==phiMax));
    }
    bool sameEta (deadIdRegion const& other) {
      return ((other.etaMin==etaMin) && (other.etaMax==etaMax));
    }
    bool samePhi (deadIdRegion const& other) {
      return ((other.phiMax==phiMax) && (other.phiMin==phiMin));
    }
    bool adjacentEta (deadIdRegion const& other) {
      return ( (other.etaMin-1 == etaMax) || 
	       (etaMin-1 == other.etaMax ) );
    }
    bool adjacentPhi (deadIdRegion const& other) {
      return ( (other.phiMin-1 == phiMax) ||
	       (phiMin-1 == other.phiMax) );
    }
    void merge (deadIdRegion const& other);
  };

  void buildDeadAreas();
  void buildSiPMAreas();
  void mergeRegionLists(std::list<deadIdRegion>& didregions);
  void convertRegions(std::list<deadIdRegion> const& idregions,
			     std::vector<deadRegion>& regions);
  TMultiGraph * graphRegions(std::vector<deadRegion> const& regions);

  std::vector<uint32_t> deadIds;
  std::vector<deadRegion> deadRegions;
  std::vector<uint32_t> SiPMIds;
  std::vector<deadRegion> SiPMRegions;
  bool inited;
  double const twopi;
  // int const etaBounds;
  // int const phiSectors;
  // std::vector<double> etaMin;
  // std::vector<double> etaMax;
  // std::vector<double> phiMinR0;
  // std::vector<double> phiMaxR0;
  // std::vector<double> phiMinR12;
  // std::vector<double> phiMaxR12;
};
