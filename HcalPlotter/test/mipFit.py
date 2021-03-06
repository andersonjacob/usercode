import sys, os
sys.path.append(os.environ['HOME']+'/pyroot')
del sys
del os

from optparse import OptionParser

parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False,
                  help='no X11 windows')
parser.add_option('-n', action='store_true', dest='noStartup', default=False,
                  help='do not run statup scripts')
parser.add_option('-c', '--calib', dest='initCalib', default='',
                  help='starting calibration')
parser.add_option('-o', dest='outputFile', default='analysis.root',
                  help='output filename')
parser.add_option('-d', '--depth', type='int', dest='depth',
                  default=1, help='use which HB depth')
parser.add_option('-s', '--sipm', action='store_true', dest='sipm',
                  default=False, help='sipm S/N assumption')
parser.add_option('--phi', type='int', default=0, dest='phi',
                  help='override table iphi position')
parser.add_option('--eta', type='int', default=0, dest='eta',
                  help='override table ieta position')
parser.add_option('--hb', action='store_true', dest='hb', default=False,
                  help='do HB instead of HO')
(opts, args) = parser.parse_args()

import root_logon
#import pyroot_fwlite.py

from ROOT import TFile, TH1F, TCanvas, TLorentzVector, kRed, kBlue, TMath, \
     TTree, gDirectory, gROOT, TMath, RooWorkspace, RooRealVar
## from array import array

from tbRoutines import *
import re

myBlue = kBlue + 2

EvtN = 0

inFile = TFile(args[0])

outFile = TFile(opts.outputFile, 'recreate')

dataTree = inFile.Get("plotanal/dataTree");

HBDigis = int(dataTree.GetMaximum('NHBdigis'))
HODigis = int(dataTree.GetMaximum('NHOdigis'))

dataTree.SetEstimate(dataTree.GetEntries())

for event in dataTree:
    break

ieta = eta2ieta(event.HBTableEta)
if (opts.eta > 0):
    ieta = opts.eta
iphi = phi2iphi(event.HBTablePhi)
if (opts.phi > 0):
    iphi = opts.phi

qualCut = '(NHOdigis=={0})&&(VMBadc>50.)'.format(HODigis)
#qualCut = '(NHOdigis=={0})&&(VMBadc>50.)&&(VMFadc>100.)&&(VMFadc<1500)&&(S3adc>0.)'.format(HODigis)
#qualCut = '(NHOdigis=={0})'.format(HODigis)
pedCut = '(triggerID==1)&&(NHOdigis=={0})'.format(HODigis)
sigCut = '(triggerID==4)&&{0}'.format(qualCut)

HOTower = 'HOE[{0}][{1}]'.format(ieta,iphi)
if opts.hb:
    HOTower = 'HBE[{0}][{1}][{2}]'.format(ieta,iphi,opts.depth)

print '{0}'.format(HOTower),\
      '(eta,phi): ({0:0.2f},{1:0.2f})'.format(event.HBTableEta,
                                              event.HBTablePhi)

dataTree.Draw('{0}'.format(HOTower), pedCut, 'goff')

havePeds = False
if (dataTree.GetSelectedRows() > 1):
    minPed = TMath.MinElement(dataTree.GetSelectedRows(), dataTree.GetV1())
    maxPed = TMath.MaxElement(dataTree.GetSelectedRows(), dataTree.GetV1())
    pedRms = TMath.RMS(dataTree.GetSelectedRows(), dataTree.GetV1())
    havePeds = True
else:
    minPed = -6.
    maxPed = 10.
    pedRms = 1.2

minPed = int(minPed) - 1.5
maxPed = int(maxPed) + 1.5

minMip = int(minPed) - 0.5
binWidth = 1.
if opts.sipm:
    maxMip = int(maxPed + pedRms*80) + 0.5
    binWidth = 5.
else :
    maxMip = int(maxPed + pedRms*20) + 0.5
while int((maxMip-minMip)/binWidth)*binWidth < (maxMip-minMip):
    maxMip += 1
Nbins = int((maxMip - minMip)/binWidth + 0.5)


print 'ped min: {0:0.2f} max: {1:0.2f}'.format(minPed,maxPed), \
      'sig min: {0:0.1f} max: {1:0.1f}'.format(minMip, maxMip)

dataTree.Draw('{0}>>ped_hist({1},{2:0.1f},{3:0.1f})'.format(HOTower,
                                                            int(maxPed-minPed),
                                                            minPed,maxPed),
              pedCut, 'goff')
dataTree.Draw('{0}>>sig_hist({1},{2:0.1f},{3:0.1f})'.format(HOTower,Nbins,
                                                            minMip,maxMip),
              sigCut, 'goff')
ped_hist = gDirectory.Get('ped_hist')
ped_hist.SetLineColor(myBlue)

sig_hist = gDirectory.Get('sig_hist')

gROOT.ProcessLine('.L langaus.C+')
from ROOT import langaupedfit, langaupro, Double, Long, preFitHisto
from array import array

fpar = array('d', [0.]*7)
fparerr = array('d', [0.]*7)
fpar[6] = pedRms

chisqr = Double(0.)
ndf = Long(0)

if havePeds:
    fit = langaupedfit(sig_hist, ped_hist, fpar, fparerr, chisqr, ndf)
    #fit = preFitHisto(sig_hist, fpar, fparerr, chisqr, ndf)
else:
    fit = preFitHisto(sig_hist, fpar, fparerr, chisqr, ndf)

#maxx = Double(0.)
#fwhm = Double(0.)

#langaupro(fpar,maxx, fwhm)

maxx = fit.GetMaximumX(fpar[1]-2.*fpar[0], fpar[1]+2.*fpar[0])

c1 = TCanvas('c1', 'pedestal')
ped_hist.Draw()
ped = ped_hist.GetFunction('ped')
if ped:
    ped.Draw('same')

c2 = TCanvas('c2', 'signal')
sig_hist.Draw()
#c2.Update()
fit.Draw('same')

print 'for {0} MIP MPV: {1:0.2f}'.format(HOTower,maxx-fpar[5]),\
      'S/N: {0:0.2f}'.format((maxx-fpar[5])/fpar[6])
