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
parser.add_option('-R', '--refine', dest='refine', action='store_true',
                  default=False, help='refine by cutting off low tail and '+ \
                  'lowering "significance" threshold')
parser.add_option('-o', dest='outputFile', default='analysis.root',
                  help='output filename')
parser.add_option('-d', '--depths', action='store_true', dest='depths',
                  default=False, help='use 4 HB depths')
parser.add_option('--beam', dest='beamE', type='float', default=100.,
                  help='beam energy')
parser.add_option('--phi', type='int', default=0, dest='phi',
                  help='override table iphi position')
parser.add_option('--eta', type='int', default=0, dest='eta',
                  help='override table ieta position')
(opts, args) = parser.parse_args()

import root_logon
#import pyroot_fwlite.py

from ROOT import TFile, TH1F, TCanvas, TLorentzVector, kRed, kBlue, TMath, \
     TBrowser, TH2F, TTree, TH1, gDirectory
## from array import array

from tbRoutines import *
from minRoutine import runMinimization
import re
from math import sqrt

def passesCuts(event):
    #is beam
    if (event.triggerID != 4):
        return False
    # event is complete
    ## if (event.NHBdigis != 68) or (event.NHOdigis != 33) or \
    ##        (event.NEBrecHits != 1700):
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

dataTree = inFile.Get("plotanal/dataTree");

calibData = {}
calibDataHO = {}

trivialConst = [1.0]*maxDim*maxDim*maxDepth
trivialConst[0] = 1.0
if len(opts.initCalib) > 0:
    trivialConst = loadCalibration(opts.initCalib)

trivialConstHO = [1.]*maxDim*maxDim

HBdepths = 2
if (opts.depths):
    HBdepths = maxDepth

beamE = opts.beamE

TH1.SetDefaultBufferSize(10000)
prelim = TH1F('prelimE', 'prelim energy', 100, 0., -1.)
passed = 0
radius = 1

resolution = beamE*sqrt(1/beamE + 0.01)

ieta = opts.eta
iphi = opts.phi
ecalXtalieta = ieta*5-2
ecalXtaliphi = ebMaxPhi - (iphi*5-2)

(ie,ip,xie,xip) = findBeamCaloCoords(dataTree)

print 'assumed resolution: {0:0.2f}'.format(resolution)

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
    ## if EvtN > 100:
    ##     break
    if (ieta==0) or (iphi==0):
        ieta = eta2ieta(event.HBTableEta)
        iphi = phi2iphi(event.HBTablePhi)


    if EvtN%500 == 1:
        print 'record:',EvtN,
        print '(ieta,iphi):', '({0},{1})'.format(ieta,iphi),
        print 'Xtal (ieta,iphi): ({0},{1})'.format(ecalXtalieta,ecalXtaliphi)

    if not passesCuts(event):
        continue
    # don't optimize the long low tails
    HB9 = 0.
    for d in range(1, HBdepths):
        HB9 += HcalEnergyAround(event.HBE, ieta, iphi, depth=d,
                                calib=trivialConst)

    EB81 = EcalEnergyAround(event.EBE, ecalXtalieta, ecalXtaliphi, radius=4)\
           *trivialConst[0]

    if opts.refine and ((HB9+EB81) < beamE - 2.*resolution):
        # print 'HB9:',HB9
        continue

    prelim.Fill(HB9+EB81)

    passed += 1
    if (event.NEBrecHits > 1):
        if not ('EB81' in calibData.keys()):
            calibData['EB81'] = []
        calibData['EB81'].append(EB81)

    for e in range(ieta-radius, ieta+radius+1):
        for p in range(iphi-radius, iphi+radius+1):
    ## for e in range(6, 10):
    ##     for p in range(2, 6):
            for d in range(1, HBdepths):
                tmpIndex = HcalIndex(e,p,d)
                if not isInstrumented('HB', e, p, d, isHPD=(not opts.depths)):
                    # print 'hb_{0}_{1}_{2}'.format(e,p,d),opts
                    tmpIndex = -1
                if (tmpIndex > 0):
                    Ename = 'hb_{0}_{1}_{2}'.format(e,p,d)
                    if not (Ename in calibData.keys()):
                        calibData[Ename] = []
                    calibData[Ename].append(event.HBE[tmpIndex]*\
                                            trivialConst[tmpIndex])
            ## tmpIndex = HcalIndex(e,p)
            ## if (tmpIndex > 0):
            ##     Ename = 'ho_{0}_{1}'.format(e,p)
            ## if not (Ename in calibDataHO.keys()):
            ##     calibDataHO[Ename] = []
            ## calibDataHO[Ename].append(event.HOE[tmpIndex]*\
            ##                         trivialConstHO[tmpIndex])

## print 'events passed:', passed,\
##       'events kept:', len(calibData[calibData.keys()[0]])

prelim.BufferEmpty()
#prelim.Draw()

minner = runMinimization(calibData, beamE)

calibConst = list(trivialConst)
calibConstHO = list(trivialConstHO)

for par in range(1, minner.GetNumberTotalParameters()):
    parName = minner.GetParName(par)
    parVal = minner.GetParameter(par)
    parErr = minner.GetParError(par)
    parSignif = parVal/parErr
    digits = re.findall('\d+', parName)
    calIndex = -1
    initConst = 1.0
    if (parName[:2] == 'EB'):
        calIndex = 0
        print 'Ecal calib:',
    elif (parName[:2] == 'hb'):
        calieta = int(digits[0])
        caliphi = int(digits[1])
        caldepth = int(digits[2])
        calIndex = HcalIndex(calieta,caliphi,caldepth)
        print 'HB ({0},{1},{2}):'.format(calieta,caliphi,caldepth),
    elif (parName[:2] == 'ho'):
        calieta = int(digits[0])
        caliphi = int(digits[1])
        calIndex = HcalIndex(calieta,caliphi)
        print 'HO ({0},{1},4):'.format(calieta,caliphi),
    if calIndex >= 0:
        parSignif = abs(1.0-parVal)/parErr
        parVal *= trivialConst[calIndex]
        #parErr *= trivialConst[calIndex]
        print 'initial value: {0:0.3f}'.format(trivialConst[calIndex]),
        if (parSignif > 0.5):
            calibConst[calIndex] = parVal
        elif opts.refine and (parSignif > 0.25):
            calibConst[calIndex] = parVal
        print 'new value: {0:0.3f} "significance": {1:0.4f}'.format(parVal, parSignif)

BarrelBefore = TH1F("BarrelBefore", "Barrel Before", 100, 0., -1.0)
BarrelAfter = TH1F("BarrelAfter", "Barrel After", 100, 0., -1.0)
BarrelBefore.SetLineColor(myBlue)

EvtN = 0
for event in dataTree:
    EvtN += 1
    ## if EvtN > 1:
    ##     break
    if (ieta == 0) or (iphi==0):
        ieta = eta2ieta(event.HBTableEta)
        iphi = phi2iphi(event.HBTablePhi)

    if EvtN%500 == 1:
        print 'record:',EvtN,
        print '(ieta,iphi):', '({0},{1})'.format(ieta,iphi),
        print 'Xtal (ieta,iphi): ({0},{1})'.format(ecalXtalieta,ecalXtaliphi)
    
    if not passesCuts(event):
        continue

    EB81 = EcalEnergyAround(event.EBE, ecalXtalieta, ecalXtaliphi, radius=4)
    sumBefore = EB81*trivialConst[0]
    sumAfter = EB81*calibConst[0]

    for d in range(1, HBdepths):
        sumAfter += HcalEnergyAround(event.HBE, ieta, iphi, depth=d,
                                     radius=radius, calib=calibConst)
        sumBefore += HcalEnergyAround(event.HBE, ieta, iphi, depth=d,
                                      radius=radius, calib=trivialConst)
    
    BarrelBefore.Fill(sumBefore)
    BarrelAfter.Fill(sumAfter)

BarrelBefore.BufferEmpty()
c1 = TCanvas('c1', 'prelim')
prelim.Draw()
BarrelAfter.BufferEmpty()
c2 = TCanvas('c2', 'before')
BarrelBefore.Draw()
c3 = TCanvas('c3', 'after')
BarrelAfter.Draw()

storageFile = 'hb_calib_{0}_{1}.pkl'.format(ieta, iphi)
storeCalibration(calibConst, storageFile)
print 'file {0} written'.format(storageFile)
