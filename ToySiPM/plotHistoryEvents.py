
from optparse import OptionParser

parser = OptionParser(usage='usage: %prog [options] INPUTFILENAME')
parser.add_option('--plotdir', dest='plotdir', metavar='PLOTDIR',
                  help='directory for plots', default='.')
parser.add_option('-b', action='store_true', dest='batch',
                  help='run in batch mode without graphics', default=False)

(opts,args) = parser.parse_args()

from ROOT import TTree,TFile,TH1F,TH1,gPad,gROOT,TCanvas
import sys
sys.path.append('/uscms/home/andersj/pyroot')
import pyroot_logon
import re

fname = args[0]
fin = open(fname, 'r')
firstLine = fin.readline()
fin.close()
floatGroup = r'([0-9\.]*)'
Test1 = float(re.search('Test1 ' + floatGroup, firstLine).group(1))
Test2 = float(re.search('Test2 ' + floatGroup, firstLine).group(1))
pescale = float(re.search('pescale ' + floatGroup, firstLine).group(1))
print floatGroup, Test1, Test2, pescale

gev2pe = 30.*pescale

Resps = TTree()
Resps.ReadFile(fname, 'TotPe/I:LivePxIntegral/F:NP/I:ped[3]/F:Layer/I:Response/F:GenPE/I:Response2/F:GenPE2/I:respForm[50]/F:genForm[50]/F')

depths = int(Resps.GetMaximum('Layer'))
BarrelPE = TH1F("BarrelPE", "Barrel PE", 100, 0., Test1*30.*pescale)
BarrelPE2 = TH1F("BarrelPE2", "Barrel PE", 100, 0., Test2*30.*pescale)
BarrelRes = TH1F("BarrelRes", "Barrel Resolution", 100, -0.3, 0.05)
BarrelRes2 = TH1F("BarrelRes2", "Barrel Resolution", 100, -0.3, 0.05)
BarrelSum = 0.;
BarrelSum2 = 0.;
BarrelDiff = 0.;
BarrelDiff2 = 0.;
for entry in Resps:
    if (entry.Layer > depths-1):
        BarrelPE.Fill(BarrelSum)
        BarrelSum = 0.
        BarrelPE2.Fill(BarrelSum2)
        BarrelSum2 = 0.
        BarrelRes.Fill(BarrelDiff/gev2pe/Test1)
        BarrelDiff = 0.
        if (Test2 > 0):
            BarrelRes2.Fill(BarrelDiff2/gev2pe/Test2)
        else:
            BarrelRes2.Fill(0.)
        BarrelDiff2 = 0.

    if (entry.Layer < depths):
        BarrelSum += entry.Response
        BarrelSum2 += entry.Response2
        BarrelDiff += entry.Response - entry.GenPE
        BarrelDiff2 += entry.Response2 - entry.GenPE2

c1 = TCanvas()
BarrelPE.Draw()
c2 = TCanvas()
BarrelRes.Draw()
