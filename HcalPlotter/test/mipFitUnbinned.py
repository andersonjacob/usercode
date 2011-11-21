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
parser.add_option('--digi', action='store_true', dest='digi', default=False,
                  help='do HO digi instead of rechit')
(opts, args) = parser.parse_args()

import root_logon
#import pyroot_fwlite.py

from ROOT import TFile, TH1F, TCanvas, kRed, kBlue, TMath, \
     TTree, gDirectory, gROOT, RooWorkspace, RooRealVar, RooLandau, \
     RooFit, RooFFTConvPdf, RooAddPdf, RooConstVar, RooGaussian, kDashed, \
     RooArgSet, RooAbsReal, RooArgList, RooDataHist, RooHistPdf
## from array import array

from tbRoutines import *
from pedRoutines import fillDataSet, fitPed, findOnePe
import re
from math import sqrt


def findMax(func, x):

    norm = RooArgList(x)
    pdfF1 = func.asTF(norm)

    maxx1 = pdfF1.GetMaximumX()

    print 'maxx: {0:0.4f}'.format(maxx1)
                
    return maxx1

def makeHistPdf(hist, ws, x):
    theVars = RooArgList(x)
    v = RooArgSet(x)
    dataHist = RooDataHist(hist.GetName()+'_dh', 'dataHist', theVars, hist)
    hpdf = RooHistPdf(hist.GetName()+'_pdf', 'hist pdf', v, dataHist)
    getattr(ws, 'import')(hpdf)
    return hpdf

## gROOT.ProcessLine('.L langaus.C+')
## from ROOT import langaupro, Double, Long, langaupedfit, preFitHisto
## from array import array

myBlue = kBlue + 2

EvtN = 0

inFile = TFile(args[0])

outFile = TFile(opts.outputFile, 'recreate')

dataTree = inFile.Get("plotanal/dataTree");

HBDigis = int(dataTree.GetMaximum('NHBdigis'))
HODigis = int(dataTree.GetMaximum('NHOdigis'))

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


qualCut = '(NHOdigis=={0})&&(VMBadc>360)'.format(HODigis)
#qualCut = '(NHOdigis=={0})&&(VMBadc>50.)&&(VMFadc>100.)&&(VMFadc<1500)&&(S3adc>0.)'.format(HODigis)
#qualCut = '(NHOdigis=={0})'.format(HODigis)
pedCut = '(triggerID==1)&&(NHOdigis=={0})'.format(HODigis)
sigCut = '(triggerID==4)&&{0}'.format(qualCut)

HOTower = 'HOE[{0}][{1}]'.format(ieta,iphi)
if opts.hb:
    HOTower = 'HBE[{0}][{1}][{2}]'.format(ieta,iphi,opts.depth)
if opts.digi:
    #HOTower = '(HODigi[2]+HODigi[3]+HODigi[4]+HODigi[5]+HODigi[6]+HODigi[7])'
    HOTower = '(HODigi[3]+HODigi[4]+HODigi[5]+HODigi[6])'
    if opts.hb:
        HOTower = '(HBDigi[3][{0}]+HBDigi[4][{0}]+HBDigi[5][{0}]+HBDigi[6][{0}])'.format(opts.depth)

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
while int((maxPed-minPed)/2.)*2. < (maxPed-minPed):
    maxPed += 1.

minMip = int(minPed) - 0.5
binWidth = 1.
if opts.sipm:
    # maxMip = int(maxPed + pedRms*80) + 0.5
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
    if pedDS.numEntries() < 5:
        havePeds = False
    else:

        if opts.sipm:
            findOnePe(ped_hist, ws, x.GetName())
            ws.var('fped').setConstant(True)
            ws.var('peMean').setConstant(True)
        else:
            fitPed(ped_hist, ws, x.GetName())
        ws.var('pedMean').setConstant(True)
        ws.var('pedWidth').setConstant(True)

        pedDS.plotOn(xfped)
        if opts.sipm:
            ws.pdf('pedPlusOne').plotOn(xfped)
            ws.pdf('pedPlusOne').paramOn(xfped)
        else:
            ws.pdf('ped').plotOn(xfped)
            ws.pdf('ped').paramOn(xfped)

        c1 = TCanvas('c1', 'pedestal')
        xfped.Draw()

        makeHistPdf(ped_hist, ws, x)
        savePedWidth = ws.var('pedWidth').getVal()

    ws.Print()


dataTree.Draw('{0}>>sig_hist({1},{2:0.1f},{3:0.1f})'.format(HOTower,Nbins,
                                                            minMip,maxMip),
              sigCut, 'goff')
sig_hist = gDirectory.Get('sig_hist')
x.setRange(minMip, maxMip)

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

    c2 = TCanvas('c2', 'signal')

    sigDS = fillDataSet(dataTree.GetV1(), x, dataTree.GetSelectedRows(),
                        'sigDS')
    mpv = RooRealVar('mpv', 'mpv',
                     sig_hist.GetBinCenter(sig_hist.GetMaximumBin()),
                     #50.,
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

    fmip = RooRealVar('fmip', 'f_{mip}', 0.95, 0., 1.)
    if havePeds:
        if opts.sipm:
            lxgplus = RooAddPdf('lxgplus', 'lxgplus', lxg,
                                ws.pdf('pedPlusOne'), fmip)
        else:
            lxgplus = RooAddPdf('lxgplus', 'lxgplus', lxg, ws.pdf('ped'), fmip)
        ws.pdf('ped').plotOn(xf, RooFit.LineColor(kRed),
                             RooFit.LineStyle(kDashed),
                             RooFit.Normalization(sig_hist.GetEntries()/3,
                                                  RooAbsReal.Raw))
        fitter = lxgplus
    else:
        fitter = lxg
        fmip.setVal(1.0)
        fmip.setConstant(True)

    xf.Draw()
    # fitter.Print('v')
    # raise 'Stop here'
    fr = fitter.fitTo(sigDS, RooFit.Minos(False), RooFit.Save(True))
    fitter.plotOn(xf)
    if (fmip.getVal() < 0.9):
        lxgplus.plotOn(xf, RooFit.LineColor(kRed+2),
                       RooFit.LineStyle(kDashed),
                       RooFit.Components('ped*'))
        lxgplus.plotOn(xf, RooFit.LineColor(myBlue),
                       RooFit.LineStyle(kDashed),
                       RooFit.Components('lxg'))
        
    #lxgplus.paramOn(xf)
    
    ## c2 = TCanvas('c2', 'signal')
    xf.Draw()
    c2.Modified()
    c2.Update()
    ## fpar[0] = width.getVal()
    ## fpar[1] = mpv.getVal()
    ## fpar[2] = sig_hist.GetEntries()*(1-fped.getVal())
    ## fpar[3] = sigma.getVal()
    ## fpar[4] = sig_hist.GetEntries()*fped.getVal()
    ## fpar[5] = ws.var('pedMean').getVal()
    ## fpar[6] = ws.var('pedWidth').getVal()
    
    ## maxx = Double(0.)
    ## fwhm = Double(0.)

    maxx = findMax(lxg, x) #, mpv.getVal(), width.getVal())
    maxxErr = maxx*mpv.getError()/mpv.getVal()

    if havePeds:
        maxx -= ws.var('pedMean').getVal()
        maxxErr = sqrt(ws.var('pedMean').getError()**2 + maxxErr**2)

    ## langaupro(fpar,maxx, fwhm)

    fr.Print('v')
    print 'for {0} MIP MPV: {1:0.2f} +/- {2:0.2f}'.format(HOTower,maxx,
                                                          maxxErr), \
          'S/N: {0:0.2f}'.format(maxx/savePedWidth),
          #'FWHM: {0:0.2f}'.format(fwhm),\
    if opts.sipm and havePeds:
        onePE = ws.var('peMean').getVal() # - ws.var('pedMean').getVal()
        peErr = ws.var('peMean').getError()
                #sqrt(ws.var('peMean').getError()**2 +
                #     ws.var('pedMean').getError()**2)
        print 'fC/pe: {0:0.3f} +/- {1:0.3f}'.format(onePE,peErr)
    else:
        print




## c1 = TCanvas('c1', 'pedestal')
## ped_hist.Draw()
## ped = ped_hist.GetFunction('ped')
## if ped:
##     ped.Draw('same')

## c2 = TCanvas('c2', 'signal')
## sig_hist.Draw()
## #c2.Update()
## fit.Draw('same')

