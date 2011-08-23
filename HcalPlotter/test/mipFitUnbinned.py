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

from ROOT import TFile, TH1F, TCanvas, kRed, kBlue, TMath, \
     TTree, gDirectory, gROOT, RooWorkspace, RooRealVar, RooLandau, \
     RooFit, RooFFTConvPdf, RooAddPdf, RooConstVar, RooGaussian, kDashed, \
     RooArgSet, RooAbsReal
## from array import array

from tbRoutines import *
from pedRoutines import fillDataSet, fitPed
import re


def findMax(func, x, mpv, width):
    p = mpv - 0.1*width
    precision = abs(0.000001*mpv)
    step = 0.05*width
    lold = -9999.
    l = -9999.
    maxCalls = 10000
    i = 0

    norm = RooArgSet(x)

    while (abs(step) > precision) and (i < maxCalls):
        i += 1
        lold = l

        x.setVal(p + step)
        l = func.getVal(norm)
        if (l < lold):
            step = -step/10
        p += step

    if (i < maxCalls):
        maxx =  x.getVal()
    else:
        return -1

    return maxx
        

## gROOT.ProcessLine('.L langaus.C+')
## from ROOT import langaupro, Double, Long, langaupedfit, preFitHisto
## from array import array

myBlue = kBlue + 2

EvtN = 0

inFile = TFile(args[0])

outFile = TFile(opts.outputFile, 'recreate')

dataTree = inFile.Get("plotanal/dataTree");

HBDigis = 68
HODigis = 33

dataTree.SetEstimate(dataTree.GetEntries())

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


qualCut = '(NHOdigis=={0})&&(VMBadc>50.)'.format(HODigis)
#qualCut = '(NHOdigis=={0})&&(VMBadc>50.)&&(VMFadc>100.)&&(VMFadc<1500)&&(S3adc>0.)'.format(HODigis)
#qualCut = '(NHOdigis=={0})'.format(HODigis)
pedCut = '(triggerID==1)&&(NHOdigis=={0})'.format(HODigis)
sigCut = '(triggerID==4)&&{0}'.format(qualCut)

HOTower = 'HOE[{0}][{1}]'.format(ieta,iphi)
if opts.hb:
    HOTower = 'HBE[{0}][{1}][{2}]'.format(ieta,iphi,opts.depth)

## print '{0}'.format(HOTower),\
##       '(eta,phi): ({0:0.2f},{1:0.2f})'.format(event.HBTableEta,
##                                               event.HBTablePhi)

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
    maxMip = int(maxPed + pedRms*100) + 0.5
    binWidth = 5.
else :
    maxMip = int(maxPed + pedRms*20) + 0.5
while int((maxMip-minMip)/binWidth)*binWidth < (maxMip-minMip):
    maxMip += 1
Nbins = int((maxMip - minMip)/binWidth + 0.5)


print 'ped min: {0:0.2f} max: {1:0.2f}'.format(minPed,maxPed), \
      'sig min: {0:0.1f} max: {1:0.1f}'.format(minMip, maxMip)

ws = RooWorkspace('ws')
x = RooRealVar('x', 'energy', minMip, maxMip, 'fC')

dataTree.Draw('{0}>>ped_hist({1},{2:0.1f},{3:0.1f})'.format(HOTower,
                                                            int(maxPed-minPed),
                                                            minPed,maxPed),
              pedCut, 'goff')
ped_hist = gDirectory.Get('ped_hist')
ped_hist.SetLineColor(myBlue)

xfped = x.frame(minPed, maxPed, int(maxPed-minPed))

if havePeds:
    pedDS = fillDataSet(dataTree.GetV1(), x, dataTree.GetSelectedRows())
    getattr(ws, 'import')(pedDS)

    fitPed(ped_hist, ws, x.GetName())
    ws.var('pedMean').setConstant(True)
    ws.var('pedWidth').setConstant(True)

    pedDS.plotOn(xfped)
    ws.pdf('ped').plotOn(xfped)

    c1 = TCanvas('c1', 'pedestal')
    xfped.Draw()

dataTree.Draw('{0}>>sig_hist({1},{2:0.1f},{3:0.1f})'.format(HOTower,Nbins,
                                                            minMip,maxMip),
              sigCut, 'goff')
sig_hist = gDirectory.Get('sig_hist')

if (sig_hist.GetEntries > 0):
    ## fpar = array('d', [0.]*7)
    ## fparerr = array('d', [0.]*7)
    ## fpar[6] = pedRms

    ## chisqr = Double(0.)
    ## ndf = Long(0)

    ## if havePeds:
    ##     fit = langaupedfit(sig_hist, ped_hist, fpar, fparerr, chisqr, ndf)
    ##     #fit = preFitHisto(sig_hist, fpar, fparerr, chisqr, ndf)
    ## else:
    ##     fit = preFitHisto(sig_hist, fpar, fparerr, chisqr, ndf)

    sigDS = fillDataSet(dataTree.GetV1(), x, dataTree.GetSelectedRows(),
                        'sigDS')
    mpv = RooRealVar('mpv', 'mpv',
                     sig_hist.GetBinCenter(sig_hist.GetMaximumBin()),
                     ## fpar[1],
                     x.getMin(), x.getMax(), 'fC')
    width = RooRealVar('width', 'width', sig_hist.GetRMS()/5., #fpar[0],
                       0., 50., 'fC')
    lnd = RooLandau('lnd', 'lnd', x, mpv, width)
    sigma = RooRealVar('sigma', '#sigma', #fpar[3],
                       sig_hist.GetRMS()/5.,
                       0., 50., 'fC')
    mean = RooConstVar('mean', 'mean', 0.)
    res = RooGaussian('res', 'res', x, mean, sigma)

    lxg = RooFFTConvPdf('lxg', 'lxg', x, lnd, res)

    xf = x.frame(RooFit.Bins(Nbins))
    sigDS.plotOn(xf)

    fped = RooRealVar('fped', 'f_{ped}', 0., 0., 1.)
    lxgplus = RooAddPdf('lxgplus', 'lxgplus', ws.pdf('ped'), lxg, fped)
    if havePeds:
        ws.pdf('ped').plotOn(xf, RooFit.LineColor(kRed),
                             RooFit.LineStyle(kDashed),
                             RooFit.Normalization(ped_hist.GetEntries(),
                                                  RooAbsReal.Raw))
    else:
        fped.setVal(0.)
        fped.setConstant(True)

    lxgplus.fitTo(sigDS, RooFit.Minos(False))
    lxgplus.plotOn(xf)
    if (fped.getVal() > 0.1):
        lxgplus.plotOn(xf, RooFit.LineColor(kRed+2),
                       RooFit.LineStyle(kDashed),
                       RooFit.Components('ped'))
        lxgplus.plotOn(xf, RooFit.LineColor(myBlue),
                       RooFit.LineStyle(kDashed),
                       RooFit.Components('lxg'))
        
    #lxgplus.paramOn(xf)
    
    c2 = TCanvas('c2', 'signal')
    xf.Draw()
    ## fpar[0] = width.getVal()
    ## fpar[1] = mpv.getVal()
    ## fpar[2] = sig_hist.GetEntries()*(1-fped.getVal())
    ## fpar[3] = sigma.getVal()
    ## fpar[4] = sig_hist.GetEntries()*fped.getVal()
    ## fpar[5] = ws.var('pedMean').getVal()
    ## fpar[6] = ws.var('pedWidth').getVal()
    
    ## maxx = Double(0.)
    ## fwhm = Double(0.)

    maxx = findMax(lxg, x, mpv.getVal(), width.getVal())

    if havePeds:
        maxx -= ws.var('pedMean').getVal()

    ## langaupro(fpar,maxx, fwhm)

    print 'for {0} MIP MPV: {1:0.2f}'.format(HOTower,maxx), \
          'S/N: {0:0.2f}'.format(maxx/ws.var('pedWidth').getVal())
          #'FWHM: {0:0.2f}'.format(fwhm),\




## c1 = TCanvas('c1', 'pedestal')
## ped_hist.Draw()
## ped = ped_hist.GetFunction('ped')
## if ped:
##     ped.Draw('same')

## c2 = TCanvas('c2', 'signal')
## sig_hist.Draw()
## #c2.Update()
## fit.Draw('same')

