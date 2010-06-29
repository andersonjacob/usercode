#include "CalibFormats/HcalObjects/interface/HcalTPGCoder.h"

std::vector<unsigned short> HcalTPGCoder::getLinearizationLUT(HcalDetId id) const {
  std::vector<unsigned short> lut(256);
  for (unsigned char i=0; i<256; ++i) lut[i]=adc2Linear(i,id);
  return lut;
}
