import sys
import ROOT
import os
from ROOT import RooFit,RooStats
ROOT.gSystem.SetIncludePath('-I$ROOFITSYS/include')
if os.access('RooPowerFunction.cxx', os.R_OK):
    ROOT.gROOT.ProcessLine('.L RooPowerFunction.cxx+')

ROOT.gStyle.SetOptStat(111111)

################################################################################
def GetBayesianInterval(filename = "workspace.root",
                        wsname = 'myWS', interactive = False):
    # this function loads a workspace and computes
    # a Bayesian upper limit

    pInFile = ROOT.TFile(filename, "read")

    # load workspace
    pWs = pInFile.Get(wsname)
    if not pWs:
        print "workspace ", wsname, " not found" 
        return -1

    # printout workspace content
    pWs.Print()

    # load and print data from workspace
    data = pWs.data("data")
    data.Print()

    # load and print S+B Model Config
    pSbHypo = pWs.obj("SbHypo")
    pSbHypo.Print()

    # create RooStats Bayesian calculator and set parameters
    bCalc = RooStats.BayesianCalculator(data, pSbHypo)
    bCalc.SetName("myBC")
    bCalc.SetConfidenceLevel(0.95)
    bCalc.SetLeftSideTailFraction(0.0)
    #bcalc->SetIntegrationType("ROOFIT")

    # estimate credible interval
    # NOTE: unfortunate notation: the UpperLimit() name refers
    #       to the upper boundary of an interval,
    #       NOT to the upper limit on the parameter of interest
    #       (it just happens to be the same for the one-sided
    #       interval starting at 0)
    pSInt = bCalc.GetInterval()
    upper_bound = pSInt.UpperLimit()
    lower_bound = pSInt.LowerLimit()

    print "one-sided 95%.C.L. bayesian credible interval for xsec: ", "[" , lower_bound , ", " , upper_bound , "]"
    

    # make posterior PDF plot for POI
    c1 = ROOT.TCanvas("posterior")
    bCalc.SetScanOfPosterior(100)
    pPlot = bCalc.GetPosteriorPlot()
    pPlot.Draw()
    ROOT.gPad.Update()
    c1.SaveAs("bayesian_num_posterior.pdf")
    if interactive:
        raw_input("press <enter> to continue")



################################################################################
if __name__=="__main__":
    GetBayesianInterval(interactive = True)
