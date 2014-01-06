import sys
import ROOT
import os
from ROOT import RooFit,RooStats
ROOT.gSystem.SetIncludePath('-I$ROOFITSYS/include')
if os.access('RooPowerFunction.cxx', os.R_OK):
    ROOT.gROOT.ProcessLine('.L RooPowerFunction.cxx+')

ROOT.gStyle.SetOptStat(111111)

################################################################################
# this function loads a workspace and computes
# a Bayesian upper limit

def GetBayesianInterval(filename = "workspace.root",
                        wsname = 'myWS', interactive = False):

    pInFile = ROOT.TFile(filename, "read");

    # load workspace
    pWs = pInFile.Get("myWS");
    if not pWs:
        print "workspace ", wsname, " not found" 
        return -1

    # printout workspace content
    pWs.Print();

    # load and print data from workspace
    data = pWs.data("data");
    data.Print();
    
    # load and print S+B Model Config
    pSbHypo = pWs.obj("SbHypo")
    pSbHypo.Print();
    
    # create RooStats Bayesian calculator and set parameters
    # Metropolis-Hastings algorithm needs a proposal function
    sp = RooStats.SequentialProposal(10.0)
    
    mcmc = RooStats.MCMCCalculator( data, pSbHypo )
    mcmc.SetConfidenceLevel(0.95)
    mcmc.SetNumIters(100000)          # Metropolis-Hastings algorithm iterations
    mcmc.SetProposalFunction(sp)
    mcmc.SetNumBurnInSteps(500) # first N steps to be ignored as burn-in
    mcmc.SetLeftSideTailFraction(0.0)
    mcmc.SetNumBins(40) # for plotting only - does not affect limit calculation
    
    
    # estimate credible interval
    # NOTE: unfortunate notation: the UpperLimit() name refers
    #       to the upper boundary of an interval,
    #       NOT to the upper limit on the parameter of interest
    #       (it just happens to be the same for the one-sided
    #       interval starting at 0)
    pMcmcInt = mcmc.GetInterval()
    upper_bound = pMcmcInt.UpperLimit( pWs.var("xsec") )
    lower_bound = pMcmcInt.LowerLimit( pWs.var("xsec") )
    
    print "one-sided 95%.C.L. bayesian credible interval for xsec: " , "[" , lower_bound , ", " , upper_bound , "]" 
    
    # make posterior PDF plot for POI
    c1 = ROOT.TCanvas ("posterior")
    plot = RooStats.MCMCIntervalPlot(pMcmcInt)
    plot.Draw()
    c1.SaveAs("bayesian_mcmc_posterior.pdf")
    
    # make scatter plots to visualise the Markov chain
    c2 = ROOT.TCanvas("xsec_vs_beta_lumi")
    plot.DrawChainScatter( pWs.var("xsec"), pWs.var("lumi_beta"))
    c2.SaveAs("scatter_mcmc_xsec_vs_beta_lumi.pdf")
    
    c3 = ROOT.TCanvas ("xsec_vs_beta_efficiency")
    plot.DrawChainScatter( pWs.var("xsec"), pWs.var("eff_beta"))
    c3.SaveAs("scatter_mcmc_xsec_vs_beta_efficiency.pdf")
    
    c4 = ROOT.TCanvas ("xsec_vs_beta_nbkg")
    plot.DrawChainScatter( pWs.var("xsec"), pWs.var("nbkg_beta"))
    c4.SaveAs("scatter_mcmc_xsec_vs_beta_nbkg.pdf")
    
    ROOT.gPad.Update()
    if interactive:
        raw_input("press <enter> to continue")


################################################################################
if __name__=="__main__":
    GetBayesianInterval(interactive = True)


