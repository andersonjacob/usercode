// -*- C++ -*-
#include <vector>
#include <stdint.h>

class HcalAcceptanceId {
 public:
  static bool isChannelDead(uint32_t id);
  static bool inGeomAccept(double eta, double phi, double delta_eta = 0.,
			   double delta_phi = 0.);
  static bool isNotDeadGeom(double eta, double phi, double delta_eta = 0.,
			    double delta_phi = 0.);
  static void initDeadIds();

 private:

  static void buildDeadAreas();

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

  static std::vector<uint32_t> deadIds;
  static std::vector<deadRegion> deadRegions;
  static bool inited;
  static int const etaBounds;
  static double const etaMin[];
  static double const etaMax[];
  static double const twopi;
  static int const phiSectors;
  static double const phiMinR0[];
  static double const phiMaxR0[];
  static double const phiMinR12[];
  static double const phiMaxR12[];
};
