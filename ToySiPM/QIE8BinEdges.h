// -*- c++ -*-

#include <vector>
#include <algorithm>

#ifndef QIEBINEDGES_H
#define QIEBINEDGES_H

struct QIEBinEdges {
  QIEBinEdges(int n, int * binEdges) {
    for(int i=0; i<n; ++i)
      edges.push_back(binEdges[i]);
    sort(edges.begin(),edges.end());
  }
  QIEBinEdges() {
  }

  int findEdge(double value) {
    std::vector<int>::iterator ub = upper_bound(edges.begin(), edges.end(), int(value));
    return *(--iterator);
  }

private:
  std::vector<int> edges;
  
};
#endif
