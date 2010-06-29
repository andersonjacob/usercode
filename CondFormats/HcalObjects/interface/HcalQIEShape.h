#ifndef HcalQIEShape_h
#define HcalQIEShape_h

/** 
\class HcalQIEData
\author Fedor Ratnikov (UMd)
POOL object to store QIE basic shape
$Author: ratnikov
$Date: 2010/05/24 15:33:34 $
$Revision: 1.1 $
*/

#include <vector>
#include <algorithm>

// 128 QIE channels
class HcalQIEShape {
 public:
  HcalQIEShape();
  ~HcalQIEShape();
  float lowEdge (unsigned fAdc) const;
  float highEdge (unsigned fAdc) const;
  float center (unsigned fAdc) const;
  bool setLowEdges (const float fValue [64]);
  unsigned range (unsigned fAdc) const {return (fAdc >> 5) & 0x3;}
  unsigned local (unsigned fAdc) const {return fAdc & 0x3f;}
 protected:
 private:
  void expand ();
  bool setLowEdge (float fValue, unsigned fAdc);
  float mValues [257];
};

#endif
