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
parser.add_option('--doi', dest='doi', default=0, type='int',
                  help='HB depth of interest')
parser.add_option('--doiMIP', dest='doiMIP', type='float', default=1.0,
                  help='mip in fC for doi')
parser.add_option('--doiScale', dest='doiScale', type='float', default=1.0,
                  help='scale factor for the doi')
parser.add_option('--beam', dest='beamE', type='float', default=100.,
                  help='beam energy')
parser.add_option('--mipHO', dest='mipHO', type='float', default=1.0,
                  help='fC/mip in HO')
parser.add_option('--fCpe', dest='fCpe', type='float', default=3.,
                  help='fC/pe in HB')
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
EBHist = TH1F("EBHist", "EB Energy", 100, 5., opts.beamE*1.3)
EBHist.GetXaxis().SetTitle("EB energy (GeV)")
BarrelHist = TH1F("BarrelHist", "EB+HB Energy", 100, opts.beamE/3., opts.beamE*2)
BarrelHist.GetXaxis().SetTitle("EB+HB Energy (GeV)")
AllCaloHist = TH1F("AllCaloHist", "EB+HB+HO Energy", 100, opts.beamE/3., opts.beamE*2.)
AllCaloHist.GetXaxis().SetTitle("EB+HB+HO Energy (GeV)")
HBDoiHist = TH1F('HBDoiHist', 'HB depth {0} energy'.format(opts.doi), 100,
                 0., opts.beamE*30.*1.25*3./opts.fCpe)
HBDoiHist.GetXaxis().SetTitle("HB depth {0} energy (pe)".format(opts.doi))
HBDoiHist.GetXaxis().SetNdivisions(508)
HBDoiMipHist = TH1F('HBDoiMipHist', 'HB depth {0} energy'.format(opts.doi), 100,
                 0., opts.beamE*2.75)
HBDoiMipHist.GetXaxis().SetTitle("HB depth {0} energy (mip)".format(opts.doi))
HBDoiMipHist.GetXaxis().SetNdivisions(508)
HBDoiAdcHist = TH2F('HBDoiAdcHist', 'HB ADC entries',128, -0.5, 127.5, 10, -0.5, 9.5)
HBDoiAdcHist.GetXaxis().SetTitle('HB ADC codes')
HBDoiAdcHist.GetYaxis().SetTitle('time slice')
## BarrelvHO = TH2F("BarrelvHO", "HO v. EB+HB", 100, 5., opts.beamE*2,
##                  200, 2., opts.beamE*1.5);
## BarrelvHO.GetXaxis().SetTitle("EB+HB Energy (GeV)");
## BarrelvHO.GetYaxis().SetTitle("HO Energy (mips)");
## VMBHist = TH1F("VMBHist", "Back Muon Veto", 100, 0., 500.)
## HOCorr = TH1F("HOCorr", "HO Energy (corrected)", 200, 0., opts.beamE*1.5)
## HOCorr.GetXaxis().SetTitle("corrected HO energy (mips)")
## BarrelvHOCorr = TH2F("BarrelvHOCorr", "HO (corrected) v. EB+HB", 100, 5.,
##                      opts.beamE*2, 200, 2., opts.beamE*1.5);
## BarrelvHOCorr.GetXaxis().SetTitle("EB+HB Energy (GeV)");
## BarrelvHOCorr.GetYaxis().SetTitle("HO Energy (mips)");

dataTree = inFile.Get("plotanal/dataTree");

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

    if EvtN%500 == 1:
        print 'record:',EvtN,
        print '(ieta,iphi):', '({0},{1})'.format(ieta,iphi),
        print 'Xtal (ieta,iphi): ({0},{1})'.format(ecalXtalieta,ecalXtaliphi)
    
    if not passesCuts(event):
        continue

    EB81 = EcalEnergyAround(event.EBE, ecalXtalieta, ecalXtaliphi,
                            radius=4)*calibConst[0]

    ## if (EB81 > 1.0):
    ##     continue
    
    EBHist.Fill(EB81)

    HB9 = 0.
    for d in range(1, HBdepths):
        HB9 += HcalEnergyAround(event.HBE, ieta, iphi, depth=d,
                                calib=calibConst)
    if opts.doi > 0:
        HBdoi = HcalEnergyAround(event.HBE, ieta, iphi, depth=opts.doi)
        HBDoiHist.Fill(HBdoi/opts.fCpe*opts.doiScale)
        HBDoiMipHist.Fill(HBdoi/opts.doiMIP)
        for ts in range(0, 10):
            di = digiIndex(ts, opts.doi)
            if (di >= 0):
                HBDoiAdcHist.Fill(event.HBadc[di],ts)
        
    HBHist.Fill(HB9)

    BarrelHist.Fill(HB9 + EB81)
    
    HOmip = HcalEnergyAround(event.HOE, ieta, iphi, radius = 0)
    HO9 = HcalEnergyAround(event.HOE, ieta, iphi, calib=calibConstHO,
                           radius = 1)
    HOHist.Fill(HOmip/opts.mipHO)

    AllCaloHist.Fill(EB81+HB9+HO9)

    ## if (opts.NP > 0) and (opts.fCpe != 0):
    ##     pes = HOmip/opts.fCpe
    ##     newHO = correctSaturation(pes, opts.NP)*opts.fCpe
    ##     HOCorr.Fill(newHO/opts.mipHO)
    ##     BarrelvHOCorr.Fill(HB9+EB81, newHO/opts.mipHO)

    ## BarrelvHO.Fill(HB9+EB81, HOmip/opts.mipHO)

outFile.Write()
outFile.Close()
inFile.Close()

f = TFile(outFile.GetName())
HBHist = f.Get('HBHist')
HOHist = f.Get('HOHist')
EBHist = f.Get('EBHist')
BarrelHist = f.Get("BarrelHist")
AllCaloHist = f.Get("AllCaloHist")
HBDoiHist = f.Get('HBDoiHist')
HBDoiMipHist = f.Get('HBDoiMipHist')
HBDoiAdcHist = f.Get('HBDoiAdcHist')

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
c5 = TCanvas('c5', 'HB Depth {0}'.format(opts.doi))
c5.SetLogy()
HBDoiHist.Draw()
c6 = TCanvas('c6', 'HB Depth {0} mips'.format(opts.doi))
c6.SetLogy()
HBDoiMipHist.Draw()
c7 = TCanvas('c7', 'HB Depth {0} adc'.format(opts.doi))
c7.SetRightMargin(0.15)
c7.SetLogz()
HBDoiAdcHist.Draw('colz')
proj = HBDoiAdcHist.ProjectionX('_projX', 5, 6)
c9 = TCanvas('c9', 'projection')
proj.Draw()

print 'depth\t126\t127'
print '{0}\t{1:0.0f}\t{2:0.0f}'.format(opts.doi,
                                       proj.GetBinContent(proj.GetNbinsX()-1),
                                       proj.GetBinContent(proj.GetNbinsX()))
