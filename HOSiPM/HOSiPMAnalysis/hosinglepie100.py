import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
process.source = cms.Source ("PoolSource",fileNames = readFiles,
                             secondaryFileNames = secFiles)

readFiles.extend( [
       '/store/user/andersj/SinglePiE100_GEN_SIM_RECO_20fCperMIP_HO0/SinglePiE100_GEN_SIM_RECO_20fCperMIP_HO0/4f297b0d30539a2aaa62ecc89bed624f/SinglePiE100_HO0_GEN_SIM_RECO_5.root',
       '/store/user/andersj/SinglePiE100_GEN_SIM_RECO_20fCperMIP_HO0/SinglePiE100_GEN_SIM_RECO_20fCperMIP_HO0/4f297b0d30539a2aaa62ecc89bed624f/SinglePiE100_HO0_GEN_SIM_RECO_4.root',
       '/store/user/andersj/SinglePiE100_GEN_SIM_RECO_20fCperMIP_HO0/SinglePiE100_GEN_SIM_RECO_20fCperMIP_HO0/4f297b0d30539a2aaa62ecc89bed624f/SinglePiE100_HO0_GEN_SIM_RECO_3.root',
       '/store/user/andersj/SinglePiE100_GEN_SIM_RECO_20fCperMIP_HO0/SinglePiE100_GEN_SIM_RECO_20fCperMIP_HO0/4f297b0d30539a2aaa62ecc89bed624f/SinglePiE100_HO0_GEN_SIM_RECO_2.root',
       '/store/user/andersj/SinglePiE100_GEN_SIM_RECO_20fCperMIP_HO0/SinglePiE100_GEN_SIM_RECO_20fCperMIP_HO0/4f297b0d30539a2aaa62ecc89bed624f/SinglePiE100_HO0_GEN_SIM_RECO_1.root' ] )

secFiles.extend( [
    ] )

process.demo = cms.EDAnalyzer('HOWithReco',
                              outfname = \
                              cms.untracked.string("SinglePiE100Reco.root"),
                              mipE = cms.untracked.double(1.0),
                              doFit = cms.untracked.bool(False)
                              )


process.p = cms.Path(process.demo)
