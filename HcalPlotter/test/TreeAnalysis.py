import sys, os
sys.path.append(os.environ['HOME']+'/pyroot')
del sys
del os

from optparse import OptionParser

parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False,
                  help='no X11 windows')
parser.add_option('-o', dest='outputFile', default='analysis.root',
                  help='output filename')
(opts, args) = parser.parse_args()

import root_logon
#import pyroot_fwlite.py

from ROOT import TFile, TH1F, TCanvas, TLorentzVector, kRed, kBlue, TMath, \
     TBrowser, TH2F, TTree
from array import array

maxDim = 11
maxDepth = 5
ebMaxPhi = 21
ebMaxEta = 86

def eta2ieta(eta):
    ieta = int(eta/0.087)
    if (ieta > 0):
        ieta += 1
    else:
        ieta -= 1
    return ieta

def phi2iphi(phi):
    iphi = int(phi/0.087) + 1
    return iphi

def HcalIndex(ieta, iphi, depth = 0):
    if (ieta < maxDim) and (iphi < maxDim) and (depth < maxDepth) and \
       (ieta >= 0) and (iphi >= 0) and (depth >= 0):
        if depth < 1:
            index = ieta*maxDim + iphi
        else:
            index = ieta*maxDim*maxDepth + iphi*maxDepth + depth
        ## print '({0},{1}) depth {2} => index: {3}'.format(ieta,iphi,
        ##                                                  depth,index)
        return index
    else:
        return -1

def EcalIndex(ieta, iphi):
    if (ieta < ebMaxEta) and (iphi < ebMaxPhi) and \
       (ieta >= 0) and (iphi >= 0):
        return ieta*ebMaxPhi + iphi
    else:
        return -1
    
def HcalEnergyAround(hits, ieta, iphi, depth = 0, radius = 1):
    energy = 0.
    ## print "len hits:",len(hits)
    for e in range(ieta-radius,ieta+radius+1):
        for p in range(iphi-radius,iphi+radius+1):
            index = HcalIndex(e,p,depth)
            if (index >= 0) and (index < len(hits)):
                energy += hits[index]
    return energy

def EcalEnergyAround(hits, ieta, iphi, radius = 2):
    energy = 0.
    for e in range(ieta-radius,ieta+radius+1):
        for p in range(iphi-radius,iphi+radius+1):
            index = EcalIndex(e,p)
            ## print '({0},{1}) => index: {2}'.format(e,p,index)
            if (index >= 0) and (index < len(hits)):
                energy += hits[index]
    return energy

def qualityCuts(event):
    if (event.triggerID != 4):
        return False
    if (event.VMBadc > 50):
        return False
    return True

myBlue = kBlue + 2

EvtN = 0

inFile = TFile(args[0])

outFile = TFile(opts.outputFile, 'recreate')

HBHist = TH1F("HBHist", "HB Energy", 100, 0., 1000.)
HOHist = TH1F("HOHist", "HO Energy", 100, 0., 500.)
EBHist = TH1F("EBHist", "EB Energy", 100, 0., 120.)
VMBHist = TH1F("VMBHist", "Back Muon Veto", 100, 0., 500.)

dataTree = inFile.Get("plotanal/dataTree");

for event in dataTree:
    EvtN += 1
    ## if EvtN > 1:
    ##     break
    if not qualityCuts(event):
        continue
    ieta = eta2ieta(event.HBTableEta)
    iphi = phi2iphi(event.HBTablePhi)

    ecalXtalieta = ieta*5-2
    ecalXtaliphi = ebMaxPhi - (iphi*5-2)

    if EvtN%500 == 1:
        print 'record:',EvtN,
        print '(ieta,iphi):', '({0},{1})'.format(ieta,iphi),
        print 'Xtal (ieta,iphi): ({0},{1})'.format(ecalXtalieta,ecalXtaliphi)
    
    VMBHist.Fill(event.VMBadc)

    EB25 = EcalEnergyAround(event.EBE, ecalXtalieta, ecalXtaliphi)
    EBHist.Fill(EB25)
    ## rowCnt = 0
    ## for hit in event.EBE:
    ##     print '{0:0.2f} '.format(hit),
    ##     rowCnt += 1
    ##     if (rowCnt >= ebMaxPhi):
    ##         print
    ##         rowCnt = 0
    ## print "EB25: {0:0.4f}".format(EB25)

    HB9 = HcalEnergyAround(event.HBE, ieta, iphi, depth=1)
    HBHist.Fill(HB9)
    ## rowCnt = 0
    ## dCnt = 0
    ## for hit in event.HBE:
    ##     print '{0:0.3f}\t'.format(hit),
    ##     rowCnt += 1
    ##     if (rowCnt >= maxDepth):
    ##         print
    ##         rowCnt = 0
    ##         dCnt += 1
    ##     if (dCnt >= maxDim):
    ##         print '\n'
    ##         dCnt = 0
    ## print "HB9: {0:0.4f}".format(HB9)
            
    HO9 = HcalEnergyAround(event.HOE, ieta, iphi)
    HOHist.Fill(HO9)
    ## rowCnt = 0
    ## for hit in event.HOE:
    ##     print '{0:0.3f}\t'.format(hit),
    ##     rowCnt += 1
    ##     if (rowCnt >= maxDim):
    ##         print
    ##         rowCnt = 0
    ## print 'HO9: {0:0.4f}'.format(HO9)    

outFile.Write()
outFile.Close()
inFile.Close()

f = TFile(outFile.GetName())
HBHist = f.Get('HBHist')
HOHist = f.Get('HOHist')
EBHist = f.Get('EBHist')
VMBHist = f.Get('VMBHist')

HBHist.Draw()
c2 = TCanvas("c2", "c2")
VMBHist.Draw()
