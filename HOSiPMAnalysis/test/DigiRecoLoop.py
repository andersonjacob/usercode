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

from DataFormats.FWLite import Events, Handle
from ROOT import *
from array import array
from math import sqrt


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
HOPedestalRecHit = TH1F("HOPedestalRecHit", "HO (target)", 40, -1., 5.)
HOPedestalRecHit.SetLineColor(kRed)

HOTargetDigi = TH1F("HOTargetDigi", "HO (target)", 50, 0., 300.)
HOPedestalDigi = TH1F("HOPedestalDigi", "HO (target)", 50, 0., 300.)
HOPedestalDigi.SetLineColor(kRed)

ietaPed = -1*opts.ieta
iphiPed = opts.iphi + 36
if iphiPed > 72:
    iphiPed -= 72

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
    linEta = (opts.ieta > 0) if opts.ieta else opts.ieta+1
    targetEta = (opts.ieta-1)*0.087 + 0.0435
    HOSimHitSum.Fill(sumSim)
    SimTargetEt.Fill(sum11*sin(eta2theta(targetEta)))
    
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
HOSimHitSum.Draw()
# gPad.Update()
# gPad.WaitPrimitive()

c2 = TCanvas('c2', 'Rec hits')
HOPedestalRecHit.Draw()
HOTargetRecHit.Draw('same')

c3 = TCanvas('c3', 'Digis')
HOPedestalDigi.Draw()
HOTargetDigi.Draw('same')

outFile.Write()
# outFile.Close()
