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

from tbRoutines import *

def passesCuts(event):
    #is beam
    if (event.triggerID != 4):
        return False
    # event is complete
    if (event.NHBdigis != 72) or (event.NHOdigis != 34):
        return False
    # is a pion
    if (event.VMBadc > 50):
        return False
    return True

myBlue = kBlue + 2

EvtN = 0

inFile = TFile(args[0])

outFile = TFile(opts.outputFile, 'recreate')

HBHist = TH1F("HBHist", "HB Energy", 100, 0., 1000.)
HOHist = TH1F("HOHist", "HO Energy", 100, 0., 500.)
EBHist = TH1F("EBHist", "EB Energy", 60, 0., 120.)
BarrelHist = TH1F("BarrelHist", "HB+EB Energy", 100, 0., 200.)
VMBHist = TH1F("VMBHist", "Back Muon Veto", 100, 0., 500.)

dataTree = inFile.Get("plotanal/dataTree");

for event in dataTree:
    EvtN += 1
    ## if EvtN > 1:
    ##     break
    if not passesCuts(event):
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

    EB81 = EcalEnergyAround(event.EBE, ecalXtalieta, ecalXtaliphi, radius=4)
    EBHist.Fill(EB81)
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

    BarrelHist.Fill(HB9/5.183+EB81)
    
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
BarrelHist = f.Get("BarrelHist")
VMBHist = f.Get('VMBHist')

c1 = TCanvas("c1", "HB Hist")
HBHist.Draw()
c2 = TCanvas("c2", "Barrel Hist")
BarrelHist.Draw()
