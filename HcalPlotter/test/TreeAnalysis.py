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

def energyAround(hits, ieta, iphi, depth = 0, radius = 1):
    energy = 0.
    for e in range(ieta-radius,ieta+radius+1):
        if (e>0) and (e<len(hits)):
            for p in range(iphi-radius,iphi+radius+1):
                if (p > 0) and (p < len(hits[e])):
                    if depth > 0:
                        if depth < len(hits[e][p]):
                            energy += hits[ieta][iphi][depth]
                    else:
                        energy += hits[ieta][iphi]
    return energy

myBlue = kBlue + 2

EvtN = 0

inFile = TFile(args[0])

outFile = TFile(opts.outputFile, 'recreate')

HBHist = TH1F("HBHist", "HB Energy", 100, 0., 500.)
HOHist = TH1F("HOHist", "HO Energy", 100, 0., 500.)
EBHist = TH1F("EBHist", "EB Energy", 100, 0., 150.)

dataTree = inFile.Get("plotanal/dataTree");

for event in dataTree:
    EvtN += 1
    if EvtN > 10:
        break
    ieta = eta2ieta(event.HBTableEta)
    iphi = phi2iphi(event.HBTablePhi)
    if EvtN%500 == 1:
        print 'record:',EvtN,
        print '(ieta,iphi):', '({0},{1})'.format(ieta,iphi)

    ecalXtalieta = ieta*5-2
    ecalXtaliphi = iphi*5-2

    print event.HOE, len(event.HOE), event.HOE[2][4]
    
    EB25 = energyAround(event.EBE, ecalXtalieta, ecalXtaliphi, radius=2)
    EBHist.Fill(EB25)
    HB9 = energyAround(event.HBE, ieta, iphi, depth=1)
    HBHist.Fill(HB9)
    HO9 = energyAround(event.HOE, ieta, iphi)
    HOHist.Fill(HO9)
    

outFile.Write()
outFile.Close()
inFile.Close()

f = TFile(outFile.GetName())
HBHist = f.Get('HBHist')
HOHist = f.Get('HOHist')
EBHist = f.Get('EBHist')

HBHist.Draw()
