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

from tbRoutines import *
from minRoutine import runMinimization
import re

centralXtalEta = 32
centralXtalPhi = 7

def passesCuts(event):
    #is beam
    if (event.triggerID != 4):
        return False
    # event is complete
    if (event.NHBdigis != 68):
        return False
    # is an electron
    if (event.VMBadc > 50):
        return False
    if (event.HBE9 > 40):
        return False
    ## if (event.maxEtaEB != centralXtalEta):
    ##     return False
    ## if (event.maxPhiEB != centralXtalPhi):
    ##     return False
    return True

myBlue = kBlue + 2

EvtN = 0

inFile = TFile(args[0])

outFile = TFile(opts.outputFile, 'recreate')

dataTree = inFile.Get("plotanal/dataTree")

calibData = {}

trivialConst = [1.18]*ebMaxPhi*ebMaxEta
if len(opts.initCalib) > 0:
    trivialConst = loadCalibration(opts.initCalib)

calibRadius = 3

for event in dataTree:
    EvtN += 1
    ## if EvtN > 1:
    ##     break
    ieta = eta2ieta(event.HBTableEta)
    iphi = phi2iphi(event.HBTablePhi)

    ecalXtalieta = centralXtalEta
    ecalXtaliphi = centralXtalPhi
    ## ecalXtalieta = ieta*5-2
    ## ecalXtaliphi = ebMaxPhi - (iphi*5-2)

    if EvtN%500 == 1:
        print 'record:',EvtN,
        print '(ieta,iphi):', '({0},{1})'.format(ieta,iphi),
        print 'Xtal (ieta,iphi): ({0},{1})'.format(ecalXtalieta,ecalXtaliphi)

    if not passesCuts(event):
        continue

    EB81 = EcalEnergyAround(event.EBE, ecalXtalieta, ecalXtaliphi, radius=4,
                            calib = trivialConst)
    if EB81 < 75:
        continue

    for e in range(ecalXtalieta-calibRadius, ecalXtalieta+calibRadius+1):
        for p in range(ecalXtaliphi-calibRadius, ecalXtaliphi+calibRadius+1):
            tmpIndex = EcalIndex(e,p)
            if (tmpIndex > 0):
                Ename = 'eb_{0}_{1}'.format(e,p)
                if not (Ename in calibData.keys()):
                    calibData[Ename] = []
                calibData[Ename].append(event.EBE[tmpIndex]*\
                                        trivialConst[tmpIndex])
                

minner = runMinimization(calibData, 80., 0.0001)

calibConst = list(trivialConst)

for par in range(1,minner.GetNumberTotalParameters()):
    parName = minner.GetParName(par)
    parVal = minner.GetParameter(par)
    parErr = minner.GetParError(par)
    calIndex = -1
    if (parName[:2] == 'eb'):
        digits = re.findall('\d+', parName)
        ## print parName,'digits:',digits
        calieta = int(digits[0])
        caliphi = int(digits[1])
        calIndex = EcalIndex(calieta,caliphi)
        print 'EB ({0},{1}):'.format(calieta,caliphi),
        if calIndex > 0:
            parVal *= trivialConst[calIndex]
            parErr *= trivialConst[calIndex]
            if (parVal/parErr > 2.0):
                calibConst[calIndex] = parVal
            print '{0:0.3f} "significance": {1:0.4f}'.format(parVal,
                                                             parVal/parErr)

BarrelBefore = TH1F("BarrelBefore", "Barrel Before", 200, 0., 100.)
BarrelAfter = TH1F("BarrelAfter", "Barrel After", 200, 0., 100.)
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

    ecalXtalieta = event.maxEtaEB
    ecalXtaliphi = event.maxPhiEB

    if EvtN%500 == 1:
        print 'record:',EvtN,
        print '(ieta,iphi):', '({0},{1})'.format(ieta,iphi),
        print 'Xtal (ieta,iphi): ({0},{1})'.format(ecalXtalieta,ecalXtaliphi)
    
    sumBefore = EcalEnergyAround(event.EBE, ecalXtalieta, ecalXtaliphi,
                                 radius=2, calib=trivialConst)
    
    sumAfter = EcalEnergyAround(event.EBE, ecalXtalieta, ecalXtaliphi,
                                radius=2, calib=calibConst)

    BarrelBefore.Fill(sumBefore)
    BarrelAfter.Fill(sumAfter)

c1 = TCanvas('c1', 'after')
BarrelAfter.Draw()
c2 = TCanvas('c2', 'before')
BarrelBefore.Draw()

storageFile = 'eb_calib_{0}_{1}.pkl'.format(centralXtalEta,centralXtalPhi)
storeCalibration(calibConst, storageFile)
print 'file {0} written'.format(storageFile)
