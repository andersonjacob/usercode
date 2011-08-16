
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

fname = args[0]

pes = TTree()
pes.ReadFile(fname, 'LayerPES[19]/F')

BarrelPE = TH1F("BarrelPE", "Barrel PE", 100, 0., 1500.*30.)
BarrelSum = 0.;
for entry in pes:
    BarrelSum = 0.
    for i in range(0,17):
        BarrelSum += entry.LayerPES[i]
    if (BarrelSum > 300):
        BarrelPE.Fill(BarrelSum)

c1 = TCanvas()
BarrelPE.Draw()
