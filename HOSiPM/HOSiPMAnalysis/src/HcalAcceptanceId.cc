#include "HOSiPM/HOSiPMAnalysis/interface/HcalAcceptanceId.h"

#include <iostream>
#include <algorithm>
#include <list>
#include <cmath>

#include "TTree.h"

#include "DataFormats/HcalDetId/interface/HcalDetId.h"

std::vector<uint32_t> HcalAcceptanceId::deadIds(86);
std::vector<HcalAcceptanceId::deadRegion> HcalAcceptanceId::deadRegions(3);
bool HcalAcceptanceId::inited = false;
int const HcalAcceptanceId::etaBounds = 5;
double const HcalAcceptanceId::etaMin[etaBounds] = 
  {-1.262, -0.861, -0.307, 0.341, 0.885};
double const HcalAcceptanceId::etaMax[etaBounds] = 
  {-0.885, -0.341,  0.307, 0.861, 1.262};
double const HcalAcceptanceId::twopi = 2.*3.14159265358979323846;
int const HcalAcceptanceId::phiSectors = 12;
double const HcalAcceptanceId::phiMinR0[phiSectors] = {-0.16172,
						        0.3618786,
						        0.8854773,
						        1.409076116,
						        1.932674892,
						        2.456273667,
						        2.979872443,
						        3.503471219,
						        4.027069994,
						        4.55066877,
						        5.074267545,
						        5.597866321 };
double const HcalAcceptanceId::phiMaxR0[phiSectors] = { 0.317395374,
							0.84099415,
							1.364592925,
							1.888191701,
							2.411790477,
							2.935389252,
							3.458988028,
							3.982586803,
							4.506185579,
							5.029784355,
							5.55338313,
							6.076981906 };
double const HcalAcceptanceId::phiMinR12[phiSectors] = {-0.166264081,
							 0.357334694,
							 0.88093347,
							 1.404532245,
							 1.928131021,
							 2.451729797,
							 2.975328572,
							 3.498927348,
							 4.022526123,
							 4.546124899,
							 5.069723674,
							 5.59332245 };
double const HcalAcceptanceId::phiMaxR12[phiSectors] = { 0.34398862,
							 0.867587396,
							 1.391186172,
							 1.914784947,
							 2.438383723,
							 2.961982498,
							 3.485581274,
							 4.00918005,
							 4.532778825,
							 5.056377601,
							 5.579976376,
							 6.103575152 };

bool HcalAcceptanceId::isChannelDead(uint32_t id) {
  if (!inited) initDeadIds();
  std::vector<uint32_t>::const_iterator found = 
    std::find(deadIds.begin(), deadIds.end(), id);
  if ((found != deadIds.end()) && (*found == id)) return true;
  else return false;
}

bool HcalAcceptanceId::inGeomAccept(double eta, double phi, 
				    double delta_eta, double delta_phi)
{
  for (int ieta = 0; ieta<etaBounds; ++ieta) {
    if ( (eta > etaMin[ieta]+delta_eta) &&
	 (eta < etaMax[ieta]-delta_eta) ) {
      for (int iphi = 0; iphi<phiSectors; ++iphi) {
	double const * mins =  ((ieta == 2) ? phiMinR0 : phiMinR12);
	double const * maxes = ((ieta == 2) ? phiMaxR0 : phiMaxR12);
	while (phi < mins[0]) 
	  phi += twopi;
	while (phi > mins[0]+twopi)
	  phi -= twopi;
	if ( ( phi > mins[iphi] + delta_phi ) &&
	     ( phi < maxes[iphi] - delta_phi ) ) {
	  return true;
	}
      }
      return false;
    }
  }
  return false;
}

bool HcalAcceptanceId::isNotDeadGeom(double eta, double phi, 
				     double delta_eta, double delta_phi) {
  if (!inited) initDeadIds();
  int ieta = int(eta/0.087) + ((eta>0) ? 1 : -1);
  double const * mins = ((std::abs(ieta) > 4) ? phiMinR12 : phiMinR0);
  while (phi < mins[0])
    phi += twopi;
  while (phi > mins[0]+twopi)
    phi -= twopi;
  std::vector<deadRegion>::const_iterator region;
  for (region = deadRegions.begin(); region != deadRegions.end(); ++region) {
    if ( (phi < region->phiMax + delta_phi) && 
	 (phi > region->phiMin - delta_phi) &&
	 (eta < region->etaMax + delta_eta) &&
	 (eta > region->etaMin - delta_eta) )
      return false;
  }
  return true;
}

void HcalAcceptanceId::initDeadIds() {
  deadIds.clear();
  TTree * deads = new TTree("deads", "deads");
  deads->ReadFile("HOdeadnessChannels.txt", "ieta/I:iphi/I:deadness/D");
  int ieta, iphi;
  double deadness;
  deads->SetBranchAddress("ieta", &ieta);
  deads->SetBranchAddress("iphi", &iphi);
  deads->SetBranchAddress("deadness", &deadness);
  deads->Print();
  for (int i=0; i<deads->GetEntries(); ++i) {
    deads->GetEntry(i);
    if (deadness > 0.4) {
      HcalDetId did(HcalOuter,ieta,iphi,4);
      deadIds.push_back(did.rawId());
    }
  }
  std::sort(deadIds.begin(), deadIds.end());
  delete deads;
  buildDeadAreas();
  inited = true;
}

void HcalAcceptanceId::buildDeadAreas() {
  std::vector<uint32_t>::iterator did;
  std::list<deadIdRegion> didregions;
  for (did = deadIds.begin(); did != deadIds.end(); ++did) {
    HcalDetId tmpId(*did);
    didregions.push_back( deadIdRegion( tmpId.ieta(), tmpId.ieta(), 
					tmpId.iphi(), tmpId.iphi() ) );
  }
  // std::cout << "dead regions: " << didregions.size() << '\n';
  std::list<deadIdRegion>::iterator curr;
  std::list<deadIdRegion> list2;
  unsigned int startSize;
  do {
    startSize = didregions.size();
    curr = didregions.begin();
    while (curr != didregions.end()) {
      deadIdRegion merger(*curr);
      curr = didregions.erase(curr);
      // std::cout << "merger: (" << merger.etaMin << ',' << merger.phiMin << ") ("
      // 		<< merger.etaMax << ',' << merger.phiMax << ")\n";
      // std::cout << "dead regions: " << didregions.size() << '\n';
      while (curr != didregions.end()) {
	if ( (merger.sameEta(*curr)) && (merger.adjacentPhi(*curr)) ) {
	  merger.merge(*curr);
	  curr = didregions.erase(curr);
	} else ++curr;
      }
      // std::cout << "merger: (" << merger.etaMin << ',' << merger.phiMin << ") ("
      // 		<< merger.etaMax << ',' << merger.phiMax << ")\n";
      list2.push_back(merger);
      curr = didregions.begin();
    }
    // std::cout << "dead regions: " << didregions.size() << '\n'
    // 	      << "list2 regions: " << list2.size() << '\n';
    curr = list2.begin();
    while (curr != list2.end()) {
      deadIdRegion merger(*curr);
      curr = list2.erase(curr);
      while (curr != list2.end()) {
	if ( (merger.samePhi(*curr)) && (merger.adjacentEta(*curr)) ) {
	  merger.merge(*curr);
	  curr = list2.erase(curr);
	} else ++curr;
      }
      didregions.push_back(merger);
      curr = list2.begin();
    }
    // std::cout << "dead regions: " << didregions.size() << '\n'
    // 	      << "list2 regions: " << list2.size() << '\n'
    // 	      << "start size: " << startSize << '\n';
  } while (startSize > didregions.size());

  double e1, e2;
  double eMin,eMax,pMin,pMax;
  double const offset = 2.;
  std::cout << "dead regions: " << didregions.size() << '\n';
  for (curr = didregions.begin(); curr != didregions.end(); ++curr) {
    std::cout << "region boundaries: ieta,iphi\n"
	      << "              min: " << curr->etaMin << ',' << curr->phiMin << '\n'
	      << "              max: " << curr->etaMax << ',' << curr->phiMax << '\n';
    e1 = (std::abs(curr->etaMin)-1)*0.087*
      (-(curr->etaMin<0) + (curr->etaMin>0));
    e2 = std::abs(curr->etaMin)*0.087*
      (-(curr->etaMin<0) + (curr->etaMin>0));
    eMin = std::min(e1,e2);
    e1 = (std::abs(curr->etaMax)-1)*0.087*
      (-(curr->etaMax<0) + (curr->etaMax>0));
    e2 = std::abs(curr->etaMax)*0.087*
      (-(curr->etaMax<0) + (curr->etaMax>0));
    eMax = std::max(e1,e2);
    double const * mins = ((std::abs(curr->etaMin)>4) ? phiMinR12 : phiMinR0);
    pMin = (curr->phiMin-1)*0.087 + mins[0] + 0.087*offset;
    pMax = curr->phiMax*0.087 + mins[0] + 0.087*offset;
    while (pMax < mins[0])
      pMax += twopi;
    while (pMax > mins[0]+twopi) 
      pMax -= twopi;
    while (pMin < mins[0])
      pMin += twopi;
    while (pMin > mins[0]+twopi) 
      pMin -= twopi;
    deadRegion tmp(eMin, eMax, pMin, pMax);
    deadRegions.push_back(tmp);
    std::cout << "                 : eta,phi\n"
	      << "              min: " << tmp.etaMin << ',' << tmp.phiMin << '\n'
	      << "              max: " << tmp.etaMax << ',' << tmp.phiMax << '\n';
  }
}

void HcalAcceptanceId::deadIdRegion::merge (deadIdRegion const& other) {
  etaMin = std::min(etaMin, other.etaMin);
  etaMax = std::max(etaMax, other.etaMax);
  phiMin = std::min(phiMin, other.phiMin);
  phiMax = std::max(phiMax, other.phiMax);
}
