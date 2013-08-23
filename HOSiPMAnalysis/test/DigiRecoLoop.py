def unpackHcalIndex(idx):
    det  = (idx>>28)&15
    depth= (idx>>26)&3
    depth+=1
    lay  = (idx>>21)&31
    #lay+=1
    z    = (idx>>20)&1
    eta  = (idx>>10)&1023
    phi  = (idx&1023)

    z = -1 if z==0 else 1

    # if (det==2) and (eta==16):
    #     lay += 8

    return (det,z,eta,phi,depth,lay)

from math import exp, atan, sin, cos
def eta2theta(eta):
    return 2*atan(exp(-eta))

def findMax(func, x):

    norm = RooArgList(x)
    pdfF1 = func.asTF(norm)

    maxx1 = pdfF1.GetMaximumX()

    print 'max {1}: {0:0.4f}'.format(maxx1, x.GetName())
                
    return maxx1

def fillDataSet(data, x, N, dsName = 'ds'):
    cols = RooArgSet(x)
    ds = RooDataSet(dsName, dsName, cols)
    #ds.Print()
    print 'length data:', N
    for datum in range(0,N):
        if (data[datum] < x.getMax()) and (data[datum] > x.getMin()):
            x.setVal(data[datum])
            ds.add(cols)
    ds.Print()
    return ds

#2345678911234567892123456789312345678941234567895123456789612345678971234567898
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False,
                  help='no X11 windows')
parser.add_option('-o', dest='outputFile', default='DigiRecoHistograms.root',
                  help='output filename')
parser.add_option('--ieta', dest='ieta', default=1, type='int',
                  help='target ieta')
parser.add_option('--iphi', dest='iphi', default=1, type='int',
                  help='target iphi')
parser.add_option('--testNumbering', dest='testNumbers', default=False,
                  action='store_true', help='use SimHit test number unpacker')
(opts, args) = parser.parse_args()

print 'opts:',opts,'\nargs:',args

try:
    import pyroot_logon
except ImportError as e:
    print "no pyroot_logon.py file to customize ROOT defaults"

from DataFormats.FWLite import Events, Handle
from ROOT import *
from array import array
from math import sqrt
import os



files = []
if len(args) > 0:
    files = args
print 'len(files)',len(files),'first:',files[0]
events = Events (files)

simHitHandle = Handle('vector<PCaloHit>')
horecoHandle = Handle('edm::SortedCollection<HORecHit,edm::StrictWeakOrdering<HORecHit> >')
simDigiHandle = Handle('edm::SortedCollection<HODataFrame,edm::StrictWeakOrdering<HODataFrame> >')

simHitLabel = ('g4SimHits','HcalHits')
horecoLabel = ('horeco')
simDigiLabel = ('simHcalUnsuppressedDigis')

EvtN = 0

outFile = TFile(opts.outputFile, 'recreate')
HOSimHitSum = TH1F("HOSimHitSum", 'Sum HO Sim Hits', 40, 0., 0.02)
EnergyVsEta = TH2F("EnergyVsEta", "Energy vs Eta", 31, -15.5, 15.5, 50, 0., 0.05)
SimTargetEt = TH1F("SimTargetEt", "SimTargetEt", 40, 0., 0.02)

HOTargetRecHit = TH1F("HOTargetRecHit", "HO (target)", 40, -1., 5.)
HOPedestalRecHit = TH1F("HOPedestalRecHit", "HO (target)", 400, -1., 5.)
HOPedestalRecHit.SetLineColor(kRed)

HOTargetDigi = TH1F("HOTargetDigi", "HO (target)", 30, 0., 300.)
HOPedestalDigi = TH1F("HOPedestalDigi", "HO (target)", 300, 0., 300.)
HOPedestalDigi.SetLineColor(kRed)

ietaPed = -1*opts.ieta
iphiPed = opts.iphi + 36
if iphiPed > 72:
    iphiPed -= 72

data = []

linEta = opts.ieta if (opts.ieta > 0) else opts.ieta+1
targetEta = (linEta-1)*0.087 + 0.0435

print 'linEta', linEta, 'eta:', targetEta, 'theta:', eta2theta(targetEta)

for event in events:
    ## if EvtN > 9:
    ##     break
    EvtN += 1
    if EvtN%500 == 1:
        print 'record:',EvtN,'Run:',event.object().id().run(),\
              'event:',event.object().id().event()

    event.getByLabel(simHitLabel, simHitHandle)
    simHits = simHitHandle.product()

    sumSim = 0.
    sum11 = 0.
    sumHB = 0.
    sumHBTarget = 0.
    sumHE = 0.
    for hit in simHits:
        if opts.testNumbers:
            (det, z, ieta, iphi, depth, layer) = unpackHcalIndex(hit.id())
        else:
            hid = HcalDetId(hit.id())
            det = hid.subdet()
            ieta = hid.ieta()
            iphi = hid.iphi()
        if (det == 3):
            sumSim += hit.energy()
            if (ieta == opts.ieta) and (iphi == opts.iphi):
                sum11 += hit.energy()
        elif (det == 1):
            sumHB += hit.energy()
            if (ieta == opts.ieta) and (iphi == opts.iphi):
                sumHBTarget += hit.energy()
        elif (det == 2):
            sumHE += hit.energy()
            if (ieta == opts.ieta) and (iphi == opts.iphi):
                sumHBTarget += hit.energy()
    HOSimHitSum.Fill(sumSim)
    SimTargetEt.Fill(sum11*sin(eta2theta(targetEta)))
    data.append(sum11)
    
    event.getByLabel(horecoLabel, horecoHandle)
    horeco = horecoHandle.product()

    for hohit in horeco:
        if (hohit.id().ieta() == opts.ieta) and (hohit.id().iphi() == opts.iphi):
            HOTargetRecHit.Fill(hohit.energy()*sin(eta2theta(targetEta)))
        if (hohit.id().ieta() == ietaPed) and (hohit.id().iphi() == iphiPed):
            HOPedestalRecHit.Fill(hohit.energy()*sin(eta2theta(-targetEta)))

    event.getByLabel(simDigiLabel, simDigiHandle)
    simDigis = simDigiHandle.product()

    for simDigi in simDigis:
        if (simDigi.id().ieta() == opts.ieta) and (simDigi.id().iphi()==opts.iphi):
            HOTargetDigi.Fill(simDigi[5].nominal_fC()+simDigi[6].nominal_fC())
        if (simDigi.id().ieta() == ietaPed) and (simDigi.id().iphi()==iphiPed):
            HOPedestalDigi.Fill(simDigi[5].nominal_fC()+simDigi[6].nominal_fC())
        

print 'total records processed:',EvtN

HOSimHitSum.Print()
# HOSimHitSum.Draw()
# gPad.Update()
# gPad.WaitPrimitive()

c1 = TCanvas('c1', 'Sim hits')
SimTargetEt.Draw()
SimTargetEt.Draw()

c2 = TCanvas('c2', 'Rec hits')
HOPedestalRecHit.Draw()
HOTargetRecHit.Draw('same')
HOTargetRecHit.Print()

c3 = TCanvas('c3', 'Digis')
HOPedestalDigi.Draw()
HOTargetDigi.Draw('same')
HOTargetDigi.Print()

outFile.cd()
outFile.Write()

thews = RooWorkspace('theWS', 'theWS')
simE = thews.factory('simE[0., 0.05]')

ds = fillDataSet(data, simE, len(data))
getattr(thews, 'import')(ds)

mpv = thews.factory('mpv[0.005, 0., 5.]')
width = thews.factory('width[0.0003, 0., 5.]')
lnd = thews.factory('RooLandau::mip(simE, mpv, width)')

sigma = thews.factory('sigma[0.0001, 0., 5.]')
res = thews.factory('RooGaussian::res(simE, 0., sigma)')

if opts.ieta > 4:
    mpv.setVal(0.003)
    sigma.setVal(0.00001)
    width.setVal(0.0001)

lndgaus = thews.factory('RooFFTConvPdf::lndgaus(simE, mip, res)')

# fitter = lnd
fitter = lndgaus
fitter.Print()

fitter.fitTo(ds)

c4 = TCanvas('c4', 'Fit')
xf = simE.frame(RooFit.Bins(50), RooFit.Range(0., 0.02))
xf.SetName('SimHitEFit')
ds.plotOn(xf)
fitter.plotOn(xf)
xf.Draw()

Emax = findMax(fitter, simE)
print 'max ET:', Emax*sin(eta2theta(targetEta))

thews.Write()
xf.Write()
# outFile.Close()
