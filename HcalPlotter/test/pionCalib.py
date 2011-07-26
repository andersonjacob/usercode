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
(opts, args) = parser.parse_args()

import root_logon
#import pyroot_fwlite.py

from ROOT import TFile, TH1F, TCanvas, TLorentzVector, kRed, kBlue, TMath, \
     TBrowser, TH2F, TTree
## from array import array

from tbRoutines import *
from minRoutine import runMinimization
import re

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

dataTree = inFile.Get("plotanal/dataTree");

calibData = {}
calibDataHO = {}

trivialConst = [100./518.3]*maxDim*maxDim*maxDepth
trivialConst[0] = 1.4
if len(opts.initCalib) > 0:
    trivialConst = loadCalibration(opts.initCalib)

trivialConstHO = [1./2.]*maxDim*maxDim

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
    
    EB81 = EcalEnergyAround(event.EBE, ecalXtalieta, ecalXtaliphi, radius=4)\
           *trivialConst[0]

    if not ('EB81' in calibData.keys()):
        calibData['EB81'] = []
    calibData['EB81'].append(EB81)

    for e in range(ieta-1, ieta+2):
        for p in range(iphi-1, iphi+2):
            tmpIndex = HcalIndex(e,p,1)
            if (tmpIndex > 0):
                Ename = 'hb_{0}_{1}_{2}'.format(e,p,1)
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

minner = runMinimization(calibData, 100.)

calibConst = list(trivialConst)
calibConstHO = list(trivialConstHO)

for par in range(1, minner.GetNumberTotalParameters()):
    parName = minner.GetParName(par)
    parVal = minner.GetParameter(par)
    digits = re.findall('\d+', parName)
    calIndex = -1
    initConst = 1.0
    if (parName[:2] == 'EB'):
        calIndex = 0
        print 'Ecal calib:',
        if calIndex >= 0:
            parVal *= trivialConst[calIndex]
            calibConst[calIndex] = parVal
            print '{0:0.3f} index: {1}'.format(parVal, calIndex)
    elif (parName[:2] == 'hb'):
        calieta = int(digits[0])
        caliphi = int(digits[1])
        caldepth = int(digits[2])
        calIndex = HcalIndex(calieta,caliphi,caldepth)
        print 'HB ({0},{1},{2}):'.format(calieta,caliphi,caldepth),
        if calIndex > 0:
            parVal *= trivialConst[calIndex]
            calibConst[calIndex] = parVal
            print '{0:0.3f} index: {1}'.format(parVal, calIndex)
    elif (parName[:2] == 'ho'):
        calieta = int(digits[0])
        caliphi = int(digits[1])
        calIndex = HcalIndex(calieta,caliphi)
        print 'HO ({0},{1},4):'.format(calieta,caliphi),
        if calIndex > 0:
            parVal *= trivialConstHO[calIndex]
            calibConstHO[calIndex] = parVal
            print '{0:0.3f} index: {1}'.format(parVal, calIndex)

BarrelBefore = TH1F("BarrelBefore", "Barrel Before", 100, 0., 200.)
BarrelAfter = TH1F("BarrelAfter", "Barrel After", 100, 0., 200.)
BarrelBefore.SetLineColor(myBlue)

EvtN = 0
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
    
    EB81 = EcalEnergyAround(event.EBE, ecalXtalieta, ecalXtaliphi, radius=4)
    sumBefore = EB81*trivialConst[0]
    sumAfter = EB81*calibConst[0]

    sumAfter += HcalEnergyAround(event.HBE, ieta, iphi, depth=1,
                                 calib=calibConst)
    sumBefore += HcalEnergyAround(event.HBE, ieta, iphi, depth=1,
                                  calib=trivialConst)
    
    BarrelBefore.Fill(sumBefore)
    BarrelAfter.Fill(sumAfter)

c1 = TCanvas('c1', 'after')
BarrelAfter.Draw()
c2 = TCanvas('c2', 'before')
BarrelBefore.Draw()

storageFile = 'hb_calib_{0}_{1}.pkl'.format(ieta, iphi)
storeCalibration(calibConst, storageFile)
print 'file {0} written'.format(storageFile)
