import sys, os
sys.path.append(os.environ['HOME']+'/pyroot')
del sys
del os

#gStyle.SetOptTitle(0)
from array import array
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-e", "--energy", metavar="TOWERCNT", dest='towers',
                  type='int', help="number of towers 1,9,25,49",
                  default=9)
parser.add_option("-m", "--mip", metavar="MIPE", dest='mip',
                  type='float', help='RecHit GeV/MIP', default=0.81)
(options, args) = parser.parse_args()

import pyroot_logon
import root_logon

from ROOT import gROOT

gROOT.ProcessLine('.L truncRMS.cc+')

from ROOT import TFile, TTree, gPad, TFitter, kBlue, kGreen, \
     Double, Long, \
     setData, hookupMinuit, makeHOFunc

dataf = TFile(args[0])
hcalhits = TTree()
dataf.GetObject("hcalhits", hcalhits)

mipE = options.mip
lblHB = 'HB_E'
lblEB = 'EB_E'
lblHO = '(HO_E'
cuts = '(HB_E9 + EB_E9 > 5)'
if (options.towers > 1):
    lblHB += str(options.towers)
    lblEB += str(options.towers)
    lblHO += str(options.towers)
lblHO += '/' + str(mipE) + ')'

barrelH = TH1D('barrelH', 'barrelH', 100, 0., hcalhits.GetMaximum(lblHB))
hcalhits.Draw(lblEB + ' + ' + lblHB + '>>barrelH', cuts)
gPad.Update()

sume = 0.
sumw = 0.

for i in range(barrelH.GetMaximumBin()-2, barrelH.GetMaximumBin()+3):
    sume += barrelH.GetBinCenter(i)*barrelH.GetBinContent(i)
    sumw += barrelH.GetBinContent(i)

mpv = sume/sumw
truncRes = sqrt(1.2*1.2/mpv + 0.069*0.069) 
## setCutoff(truncRes)
## setRange(mpv-3*barrelH.GetRMS(),
##          mpv+4*barrelH.GetRMS())
#setMPV(mpv)


#'((' + lblHB + '+' + lblEB + ') < ' + str(mpv*(1.+5.*truncRes)) + ') &&' +
dataCut = '((' + lblHB + '+' + lblEB + ') > ' + str(mpv*(1.-2.*truncRes)) + ')'
#dataCut = '((' + lblHB + '+' + lblEB + ') > ' + str(mpv*(0.8)) + ')'

#dataCut = '((' + lblEB + '+' + lblHB + ') > 10)'
print "cuts:",cuts,"dataCut:",dataCut,"lblHO:",lblHO
hcalhits.Draw(lblHO + ':' + lblHB + ':' + lblEB, cuts + ' && ' + dataCut)
gPad.Update()
setData(hcalhits.GetSelectedRows(), hcalhits.GetV2(), hcalhits.GetV1(),
        hcalhits.GetV3())

lastRms = truncRes*mpv

print mpv, truncRes, Ndata, lastRms

minner = TFitter(7)
pl = array('d', [1.])
errdef = array('d', [0.001])
minner.ExecuteCommand('SET PRINT', pl, 1)
minner.ExecuteCommand('SET ERR', errdef, 1)

#setPrintLevel(minner, 0)

strat = array('d', [1.])
#minner.ExecuteCommand('SET STR', strat, 1)

hookupMinuit(minner)
#minner.SetFCN(truncRMS)

minner.SetParameter(0, 'HO_wgt', 0., 0.1, -0.5, 2.0)
minner.FixParameter(0)
minner.SetParameter(1, 'EB_wgt', 1.2, 0.1, 0.75, 2.0)
minner.SetParameter(2, 'MPV', mpv, mpv*truncRes,
                    mpv*(1-3*truncRes), mpv*(1+3*truncRes))

arglist = array('d', [1000., 0.1])
#minuit = minner.GetMinuit()
#minuit.Migrad()
#minuit.mnhess()
minner.ExecuteCommand('SIMPLEX', arglist, 2)
arglist[1] = 1.
minner.ExecuteCommand('MINIMIZE', arglist, 2)
minner.ExecuteCommand('HESSE', arglist, 1)

## ho_wgt = minner.GetParameter(0)
## eb_wgt = minner.GetParameter(1)
## hcalhits.Draw(str(eb_wgt) + '*' + lblEB + '+' + lblHB + '>>htemp', cuts)
## sume = 0.;
## sumw = 0.;
## for i in range(htemp.GetMaximumBin()-2, htemp.GetMaximumBin()+3):
##     sume += htemp.GetBinCenter(i)*htemp.GetBinContent(i)
##     sumw += htemp.GetBinContent(i)
## mpv = sume/sumw
## print mpv

#minner.FixParameter(2)
minner.ReleaseParameter(0)
minner.SetParameter(0, 'HO_wgt', 1.0, 0.05, 0., 2.5)

#minner.FixParameter(1)
#minner.ExecuteCommand('SIMPLEX', arglist, 2)
#minner.ExecuteCommand('HESSE', arglist, 1)
#minner.ExecuteCommand('MINIMIZE', arglist, 2)
#minner.ExecuteCommand('HESSE', arglist, 1)
#minner.ReleaseParameter(1)

arglist[0] = 1500.
arglist[1] = 0.01
minner.ExecuteCommand('SIMPLEX', arglist, 2)
minner.ExecuteCommand('HESSE', arglist, 1)

## ho_wgt = minner.GetParameter(0)
## eb_wgt = minner.GetParameter(1)
## hcalhits.Draw(str(eb_wgt) + '*' + lblEB + '+' + lblHB + '+' + \
##               str(ho_wgt) + '*' + lblHO + '>>htemp2',
##               cuts)
## sume = 0.;
## sumw = 0.;
## for i in range(htemp2.GetMaximumBin()-2, htemp2.GetMaximumBin()+3):
##     sume += htemp2.GetBinCenter(i)*htemp2.GetBinContent(i)
##     sumw += htemp2.GetBinContent(i)
## mpv = sume/sumw
#setMPV(mpv)

## print mpv

arglist[1] = 1
minOut =  minner.ExecuteCommand('MINIMIZE', arglist, 2)
if (minOut == 0):
    minner.ExecuteCommand('HESSE', arglist, 1)
amin = Double(1.0)
edm = Double(0.1)
errdef = Double(0.01)
nvpar = Long(2)
nparx = Long(2)
minner.GetStats(amin, edm, errdef, nvpar, nparx)
print amin,edm,errdef,nvpar,nparx
## if (edm > 0.001):
##     minner.ExecuteCommand('MIGRAD', arglist, 2)
##     minner.ExecuteCommand('HESSE', arglist, 1)

ho_wgt = minner.GetParameter(0)
eb_wgt = minner.GetParameter(1)
mpv = minner.GetParameter(2)
finalRes = sqrt(1.2*1.2/mpv + 0.069*0.069) 

print "eb_wgt:",eb_wgt,"ho_wgt:",ho_wgt,"mpv:",mpv,"finalResTB06:",finalRes


hcalhits.Draw(lblEB + '+' + lblHB + '>>simpleBE', cuts)
maxHist = simpleBE.GetMaximum()
gPad.Update()

hcalhits.Draw(str(eb_wgt) + '*' + lblEB + '+' + lblHB + '>>optBE',
              cuts)
optBE.SetLineColor(kBlue)
if (optBE.GetMaximum() > maxHist):
    maxHist = optBE.GetMaximum()
gPad.Update()
## hcalhits.Draw(lblEB + '+' + lblHB + '+' + lblHO + '>>simpleSum')
## simpleSum.SetLineColor(kRed)
## if (simpleSum.GetMaximum() > maxHist):
##     maxHist = simpleSum.GetMaximum()
## gPad.Update()
hcalhits.Draw(str(eb_wgt) + '*' + lblEB + '+' + lblHB + '+' + \
              str(ho_wgt) + '*' + lblHO + '>>optSum',
              cuts)
optSum.SetLineColor(kGreen)
gPad.Update()
if (optSum.GetMaximum() < maxHist):
    optSum.SetMaximum(maxHist*1.05)

print "maximum:",maxHist
optSum.Draw('')
#simpleSum.Draw('same')
optBE.Draw('same')
#simpleBE.Draw('same')
gPad.SetLogy()
gPad.Update()

hofunc = makeFuncHO()
hofunc.SetParameter(0, eb_wgt)
hofunc.SetParameter(1, mpv)
c2 = TCanvas("c2","c2")
hofunc.Draw()
