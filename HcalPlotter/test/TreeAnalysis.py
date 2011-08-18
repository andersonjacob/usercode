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
parser.add_option('-n', action='store_true', dest='noStartup', default=False,
                  help='do not run statup scripts')
parser.add_option('-c', '--calib', dest='initCalib', default='',
                  help='starting calibration')
parser.add_option('-d', '--depths', action='store_true', dest='depths',
                  default=False, help='use 4 HB depths')
parser.add_option('--beam', dest='beamE', type='float', default=100.,
                  help='beam energy')
parser.add_option('--mipHO', dest='mipHO', type='float', default=1.0,
                  help='fC/mip in HO')
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
    if (event.NHBdigis != 68) or (event.NHOdigis != 33) or \
           (event.NEBrecHits != 1700):
        return False
    ## # is a pion
    ## if (event.VMBadc > 50):
    ##     return False
    return True

myBlue = kBlue + 2

EvtN = 0

inFile = TFile(args[0])

outFile = TFile(opts.outputFile, 'recreate')

HBHist = TH1F("HBHist", "HB Energy", 100, 5., opts.beamE*1.6)
HBHist.GetXaxis().SetTitle("HB energy (GeV)")
HOHist = TH1F("HOHist", "HO Energy", 200, 0., opts.beamE*1.5)
HOHist.GetXaxis().SetTitle("HO energy (mips)")
EBHist = TH1F("EBHist", "EB Energy", 100, 5., opts.beamE*1.3)
EBHist.GetXaxis().SetTitle("EB energy (GeV)")
BarrelHist = TH1F("BarrelHist", "HB+EB Energy", 100, 5., opts.beamE*2)
BarrelHist.GetXaxis().SetTitle("HB+EB Energy (GeV)")
BarrelvHO = TH2F("BarrelvHO", "HO v. EB+HB", 100, 5., opts.beamE*2,
                 200, 2., opts.beamE*1.5);
BarrelvHO.GetXaxis().SetTitle("EB+HB Energy (GeV)");
BarrelvHO.GetYaxis().SetTitle("HO Energy (mips)");
VMBHist = TH1F("VMBHist", "Back Muon Veto", 100, 0., 500.)

dataTree = inFile.Get("plotanal/dataTree");

calibConst = [1.0]*maxDim*maxDim*maxDepth
calibConstHO = [1.]*maxDim*maxDim

if len(opts.initCalib) > 0:
    calibConst = loadCalibration(opts.initCalib)

HBdepths = 2
if (opts.depths):
    HBdepths = maxDepth

for event in dataTree:
    EvtN += 1
    ## if EvtN > 2:
    ##     break
    ieta = eta2ieta(event.HBTableEta)
    iphi = phi2iphi(event.HBTablePhi)

    ecalXtalieta = ieta*5-2
    ecalXtaliphi = ebMaxPhi - (iphi*5-2)

    if EvtN%500 == 1:
        print 'record:',EvtN,
        print '(ieta,iphi):', '({0},{1})'.format(ieta,iphi),
        print 'Xtal (ieta,iphi): ({0},{1})'.format(ecalXtalieta,ecalXtaliphi)
    
    if not passesCuts(event):
        continue

    VMBHist.Fill(event.VMBadc)

    EB81 = EcalEnergyAround(event.EBE, ecalXtalieta, ecalXtaliphi,
                            radius=4)*calibConst[0]

    ## if (EB81 > 1.0):
    ##     continue
    
    EBHist.Fill(EB81)

    ## rowCnt = 0
    ## for hit in event.EBE:
    ##     print '{0:0.2f} '.format(hit),
    ##     rowCnt += 1
    ##     if (rowCnt >= ebMaxPhi):
    ##         print
    ##         rowCnt = 0
    ## print "EB25: {0:0.4f}".format(EB25)

    HB9 = 0.
    for d in range(1, HBdepths):
        HB9 += HcalEnergyAround(event.HBE, ieta, iphi, depth=d,
                                calib=calibConst)
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

    BarrelHist.Fill(HB9 + EB81)
    
    HOmip = HcalEnergyAround(event.HOE, ieta, iphi, calib=calibConstHO,
                             radius = 0)/opts.mipHO
    HOHist.Fill(HOmip)

    BarrelvHO.Fill(HB9+EB81, HOmip)
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
BarrelvHO = f.Get('BarrelvHO')
VMBHist = f.Get('VMBHist')

c1 = TCanvas("c1", "HB Energy")
HBHist.Draw()
c2 = TCanvas("c2", "Barrel Energy")
BarrelHist.Draw()
c3 = TCanvas('c3', 'HO MIPs')
c3.SetLogy()
HOHist.Draw()
c4 = TCanvas('c4', 'EB Energy')
EBHist.Draw()
c5 = TCanvas('c5', 'HO v. Barrel')
c5.SetLogz()
BarrelvHO.Draw('boxcolz')

