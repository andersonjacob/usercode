import sys
import os
import ROOT
from ROOT import RooFit,RooStats
ROOT.gSystem.SetIncludePath('-I$ROOFITSYS/include')
if os.access('RooPowerFunction.cxx', os.R_OK):
    ROOT.gROOT.ProcessLine('.L RooPowerFunction.cxx+')

def MakeWorkspace(outfilename):
    # use this area to implement your model for a counting experiment
    print "building counting model..."

    #create workspace
    pWs = ROOT.RooWorkspace("myWS")

    #observable
    pWs.factory("n[0]") #this create a literal, real-valued container called n

    #signal yield
    # pWs.factory("nsig[0,0,100]") # name[value,min bound,max bound]
    # pWs.factory("lumi[0]") #luminosity
    pWs.factory("lumi_nom[20000, 10000, 30000]") #nominal lumi of 20 fb-1
    pWs.factory("lumi_kappa[1.026]")
    pWs.factory("lumi_beta[0,-5,5]")
    pWs.factory("RooPowerFunction::lumi_alpha(lumi_kappa, lumi_beta)")
    pWs.factory("prod::lumi(lumi_nom,lumi_alpha)") #the new lumi

    pWs.factory("Gaussian::lumi_constr(lumi_beta,lumi_glob[0,-5,5],1)")

    pWs.factory("xsec[0,0,0.1]") #cross-section (our new POI)

    # pWs.factory("eff[0]") #efficiency
    pWs.factory("eff_nom[0.1, 0.01, 0.20]") #nominal eff
    pWs.factory("eff_kappa[1.1]")
    pWs.factory("eff_beta[0,-5,5]")
    pWs.factory("RooPowerFunction::eff_alpha(eff_kappa, eff_beta)")
    pWs.factory("prod::eff(eff_nom,eff_alpha)") #the new eff

    pWs.factory("Gaussian::eff_constr(eff_beta,eff_glob[0,-5,5],1)")

    pWs.factory("prod::nsig(lumi,xsec,eff)") #new definition of nsig

    #background yield
    # pWs.factory("nbkg[10,0,100]")
    pWs.factory("nbkg_nom[10.,5.,15.]")
    pWs.factory("nbkg_kappa[1.1]")
    pWs.factory("nbkg_beta[0,-5,5]")
    pWs.factory("RooPowerFunction::nbkg_alpha(nbkg_kappa, nbkg_beta)")
    pWs.factory("prod::nbkg(nbkg_nom,lumi_alpha,nbkg_alpha)")

    pWs.factory("Gaussian::nbkg_constr(nbkg_beta, nbkg_glob[0,-5,5],1)")

    #full yield
    pWs.factory("sum::yield(nsig,nbkg)") # creates a function nsig+nbkg


    # Bayesian prior
    pWs.factory("Uniform::prior(xsec)")

    # define the physics model a Poisson distribution
    pWs.factory("Poisson::model_core(n,yield)")

    # model and systematics
    # pWs.factory("PROD::model(model_core,lumi_constr)") # all-cap PROD means PDF instead of only function
    pWs.factory("PROD::model(model_core,lumi_constr,eff_constr,nbkg_constr)")

    # create set of observables (needed for datasets and ModelConfig later)
    pObs = pWs.var("n") # get the pointer to the observable
    obs = ROOT.RooArgSet("observables")
    obs.add(pObs)

    # create the dataset
    pObs.setVal(11) # this is your observed data: you counted eleven events
    data = ROOT.RooDataSet("data", "data", obs)
    data.add(obs)

    # import dataset into workspace
    getattr(pWs,'import')(data) # we call it this way because "import" is a reserved word in python

    # create set of global observables (need to be defined as constants!)
    pWs.var("lumi_glob").setConstant(True)
    pWs.var("eff_glob").setConstant(True)
    pWs.var("nbkg_glob").setConstant(True)
    globalObs = ROOT.RooArgSet ("global_obs")
    globalObs.add( pWs.var("lumi_glob") )
    globalObs.add( pWs.var("eff_glob") )
    globalObs.add( pWs.var("nbkg_glob") )

    # create set of parameters of interest (POI)
    poi = ROOT.RooArgSet ("poi")
    poi.add( pWs.var("xsec") )

    # create set of nuisance parameters
    nuis = ROOT.RooArgSet("nuis")
    nuis.add( pWs.var("lumi_beta") )
    nuis.add( pWs.var("eff_beta") )
    nuis.add( pWs.var("nbkg_beta") )

    # fix all other variables in model:
    # everything except observables, POI, and nuisance parameters
    # must be constant
    pWs.var("lumi_nom").setConstant(True)
    pWs.var("eff_nom").setConstant(True)
    pWs.var("nbkg_nom").setConstant(True)
    pWs.var("lumi_kappa").setConstant(True)
    pWs.var("eff_kappa").setConstant(True)
    pWs.var("nbkg_kappa").setConstant(True)
    fixed = ROOT.RooArgSet ("fixed")
    fixed.add( pWs.var("lumi_nom") )
    fixed.add( pWs.var("eff_nom") )
    fixed.add( pWs.var("nbkg_nom") )
    fixed.add( pWs.var("lumi_kappa") )
    fixed.add( pWs.var("eff_kappa") )
    fixed.add( pWs.var("nbkg_kappa") )

    # create signal+background Model Config
    sbHypo = RooStats.ModelConfig("SbHypo")
    sbHypo.SetWorkspace( pWs )
    sbHypo.SetPdf( pWs.pdf("model") )
    sbHypo.SetObservables( obs )
    sbHypo.SetGlobalObservables( globalObs )
    sbHypo.SetParametersOfInterest( poi )
    sbHypo.SetNuisanceParameters( nuis )
    sbHypo.SetPriorPdf( pWs.pdf("prior") ) # this is optional

    # set parameter snapshot that corresponds to the best fit to data
    pNll = sbHypo.GetPdf().createNLL(data)
    pProfile = pNll.createProfile( globalObs ) # do not profile global observables
    pProfile.getVal() # this will do fit and set POI and nuisance parameters to fitted values
    pPoiAndNuisance = ROOT.RooArgSet("poiAndNuisance")
    pPoiAndNuisance.add(sbHypo.GetNuisanceParameters())
    pPoiAndNuisance.add(sbHypo.GetParametersOfInterest())
    sbHypo.SetSnapshot(pPoiAndNuisance)

    # import S+B ModelConfig into workspace
    getattr(pWs,'import')(sbHypo)

    # create background-only Model Config from the S+B one
    bHypo = RooStats.ModelConfig(sbHypo)
    bHypo.SetName("BHypo")
    bHypo.SetWorkspace(pWs)

    # set parameter snapshot for bHypo, setting xsec=0
    # it is useful to understand how this block of code works
    # but you can also use it as a recipe to make a parameter snapshot
    pNll = bHypo.GetPdf().createNLL( data )
    poiAndGlobalObs = ROOT.RooArgSet ("poiAndGlobalObs")
    poiAndGlobalObs.add( poi )
    poiAndGlobalObs.add( globalObs )
    pProfile = pNll.createProfile( poiAndGlobalObs ) # do not profile POI and global observables
    poi.first().setVal( 0 )  # set xsec=0 here
    pProfile.getVal() # this will do fit and set nuisance parameters to profiled values
    pPoiAndNuisance = ROOT.RooArgSet( "poiAndNuisance" )
    pPoiAndNuisance.add( nuis )
    pPoiAndNuisance.add( poi )
    bHypo.SetSnapshot(pPoiAndNuisance)

    # import model config into the workspace
    getattr(pWs,'import')(bHypo)

    #print the contents of the workspace
    pWs.Print()

    #save the workspace to a file
    pWs.SaveAs(outfilename)

    return 0

################################################################################
if __name__ == "__main__":
    MakeWorkspace("workspace.root")
