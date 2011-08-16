// -*- C++ -*-

#include <vector>
#include <algorithm>
#include <limits>

#ifndef QIEBINEDGES_H
#define QIEBINEDGES_H

class QIEBinEdges {
 public:
  QIEBinEdges(int n, int const * binEdges) {
    for(int i=0; i<n; ++i) {
      edges.push_back(binEdges[i]);
    }
    sort(edges.begin(),edges.end());
  }

  QIEBinEdges() {
  }

  int findEdge(double value) {
    std::vector<int>::iterator ub = 
      upper_bound(edges.begin(), edges.end(), int(value));
    return (*(--ub));
  }

  int binWidth(double value) {
    std::vector<int>::iterator ub =
      upper_bound(edges.begin(), edges.end(), int(value));
    if (ub == edges.end()) return numeric_limits<int>::max();
    int upper = *ub;
    int lower = *(--ub);
    return upper-lower;
  }

 private:
  std::vector<int> edges;
  
};
#endif
