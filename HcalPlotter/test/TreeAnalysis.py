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
parser.add_option('--ho', dest='initCalibHO', default='',
                  help='starting HO calibration')
parser.add_option('-d', '--depths', action='store_true', dest='depths',
                  default=False, help='use 4 HB depths')
parser.add_option('--beam', dest='beamE', type='float', default=100.,
                  help='beam energy')
parser.add_option('--mipHO', dest='mipHO', type='float', default=1.0,
                  help='fC/mip in HO')
parser.add_option('--fCpe', dest='fCpe', type='float', default=0.,
                  help='fC/pe in HO')
parser.add_option('--cells', dest='NP', type='int', default=0,
                  help='number of HO SiPM pixels illuminated')
parser.add_option('--phi', type='int', default=0, dest='phi',
                  help='override table iphi position')
parser.add_option('--eta', type='int', default=0, dest='eta',
                  help='override table ieta position')
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
    if (event.NHBdigis != NHB) or (event.NHOdigis != NHO) \
           or (event.NEBrecHits != NEB):
    # if (event.NHBdigis != 72) or (event.NHOdigis != 34):
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
HOHistPE = TH1F("HOHistPE", "HO Energy", 270, 100, 2800.)
HOHistPE.GetXaxis().SetTitle("HO energy (pe)")
EBHist = TH1F("EBHist", "EB Energy", 100, 5., opts.beamE*1.3)
EBHist.GetXaxis().SetTitle("EB energy (GeV)")
BarrelHist = TH1F("BarrelHist", "EB+HB Energy", 100, opts.beamE/3., opts.beamE*2)
BarrelHist.GetXaxis().SetTitle("EB+HB Energy (GeV)")
AllCaloHist = TH1F("AllCaloHist", "EB+HB+HO Energy", 100, opts.beamE/3., opts.beamE*2)
AllCaloHist.GetXaxis().SetTitle("EB+HB+HO Energy (GeV)")
BarrelvHO = TH2F("BarrelvHO", "HO v. EB+HB", 100, 5., opts.beamE*2,
                 200, 2., opts.beamE*1.5);
BarrelvHO.GetXaxis().SetTitle("EB+HB Energy (GeV)");
BarrelvHO.GetYaxis().SetTitle("HO Energy (mips)");
VMBHist = TH1F("VMBHist", "Back Muon Veto", 100, 0., 500.)
HOCorr = TH1F("HOCorr", "HO Energy (corrected)", 200, 0., opts.beamE*1.5)
HOCorr.GetXaxis().SetTitle("corrected HO energy (mips)")
HOCorrPE = TH1F("HOCorrPE", "HO Energy (corrected)", 270, 100., 2800.)
HOCorrPE.GetXaxis().SetTitle("corrected HO energy (pe)")
BarrelvHOCorr = TH2F("BarrelvHOCorr", "HO (corrected) v. EB+HB", 100, 5.,
                     opts.beamE*2, 200, 2., opts.beamE*1.5);
BarrelvHOCorr.GetXaxis().SetTitle("EB+HB Energy (GeV)")
BarrelvHOCorr.GetYaxis().SetTitle("HO Energy (mips)")

dataTree = inFile.Get("plotanal/dataTree")
NHB = int(dataTree.GetMaximum('NHBdigis'))
NHO = int(dataTree.GetMaximum('NHOdigis'))
NEB = int(dataTree.GetMaximum('NEBrecHits'))

calibConst = [1.0]*maxDim*maxDim*maxDepth
calibConstHO = [1.]*maxDim*maxDim

if len(opts.initCalib) > 0:
    calibConst = loadCalibration(opts.initCalib)
if len(opts.initCalibHO) > 0:
    calibConstHO = loadCalibration(opts.initCalibHO)

HBdepths = 2
if (opts.depths):
    HBdepths = maxDepth

ieta = opts.eta
iphi = opts.phi
ecalXtalieta = ieta*5-2
ecalXtaliphi = ebMaxPhi - (iphi*5-2)

(ie,ip,xie,xip) = findBeamCaloCoords(dataTree)

if (ieta == 0):
    ieta = ie
    ecalXtalieta = xie
    print 'hb ieta:',ieta,'eb ieta:',ecalXtalieta

if (iphi == 0):
    iphi = ip
    ecalXtaliphi = xip
    print 'hb iphi:',iphi,'eb iphi:',ecalXtaliphi

for event in dataTree:
    EvtN += 1
    ## if EvtN > 2:
    ##     break

    if (ieta == 0) or (iphi==0):
        ieta = eta2ieta(event.HBTableEta)
        iphi = phi2iphi(event.HBTablePhi)
        ecalXtalieta = ieta*5-2
        ecalXtaliphi = ebMaxPhi - (iphi*5-2)

    if EvtN%5000 == 1:
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
    
    HOmip = HcalEnergyAround(event.HOE, ieta, iphi, radius = 0)
    HO9 = HcalEnergyAround(event.HOE, ieta, iphi, calib=calibConstHO,
                           radius = 1)
    HOHist.Fill(HOmip/opts.mipHO)
    HOHistPE.Fill(HOmip/opts.fCpe)

    AllCaloHist.Fill(EB81+HB9+HO9)

    if (opts.NP > 0) and (opts.fCpe != 0):
        pes = HOmip/opts.fCpe
        newHO = correctSaturation(pes, opts.NP)*opts.fCpe
        HOCorr.Fill(newHO/opts.mipHO)
        HOCorrPE.Fill(newHO/opts.fCpe)
        BarrelvHOCorr.Fill(HB9+EB81, newHO/opts.mipHO)

    BarrelvHO.Fill(HB9+EB81, HOmip/opts.mipHO)
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
AllCaloHist = f.Get("AllCaloHist")
BarrelvHO = f.Get('BarrelvHO')
VMBHist = f.Get('VMBHist')
HOCorr = f.Get('HOCorr')
BarrelvHOCorr = f.Get('BarrelvHOCorr')

c4 = TCanvas('c4', 'EB Energy')
EBHist.Draw()
c1 = TCanvas("c1", "HB Energy")
HBHist.Draw()
c2 = TCanvas("c2", "Barrel Energy")
BarrelHist.Draw()
c3 = TCanvas('c3', 'HO MIPs')
c3.SetLogy()
HOHist.Draw()
c8 = TCanvas('c8', 'All Calo Energy')
AllCaloHist.Draw()
c5 = TCanvas('c5', 'HO v. Barrel')
c5.SetRightMargin(0.15)
c5.SetLogz()
BarrelvHO.Draw('boxcolz')

if (opts.NP > 0) and (opts.fCpe != 0):
    c6 = TCanvas('c6', 'corrected HO MIPs')
    c6.SetLogy()
    HOCorr.Draw()
    c7 = TCanvas('c7', 'corrected HO v. Barrel')
    c7.SetRightMargin(0.15)
    c7.SetLogz()
    BarrelvHOCorr.Draw('boxcolz')
