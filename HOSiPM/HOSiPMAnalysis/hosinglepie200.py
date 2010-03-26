import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
process.source = cms.Source ("PoolSource",fileNames = readFiles,
                             secondaryFileNames = secFiles)

readFiles.extend( [
       '/store/user/andersj/SinglePiE200_GEN_SIM_RECO_20fCperMIP_HO0/SinglePiE200_GEN_SIM_RECO_20fCperMIP_HO0/ef46547bea9a0af61bcb19c36d8df40b/SinglePiE200_HO0_GEN_SIM_RECO_5.root',
       '/store/user/andersj/SinglePiE200_GEN_SIM_RECO_20fCperMIP_HO0/SinglePiE200_GEN_SIM_RECO_20fCperMIP_HO0/ef46547bea9a0af61bcb19c36d8df40b/SinglePiE200_HO0_GEN_SIM_RECO_4.root',
       '/store/user/andersj/SinglePiE200_GEN_SIM_RECO_20fCperMIP_HO0/SinglePiE200_GEN_SIM_RECO_20fCperMIP_HO0/ef46547bea9a0af61bcb19c36d8df40b/SinglePiE200_HO0_GEN_SIM_RECO_3.root',
       '/store/user/andersj/SinglePiE200_GEN_SIM_RECO_20fCperMIP_HO0/SinglePiE200_GEN_SIM_RECO_20fCperMIP_HO0/ef46547bea9a0af61bcb19c36d8df40b/SinglePiE200_HO0_GEN_SIM_RECO_2.root',
       '/store/user/andersj/SinglePiE200_GEN_SIM_RECO_20fCperMIP_HO0/SinglePiE200_GEN_SIM_RECO_20fCperMIP_HO0/ef46547bea9a0af61bcb19c36d8df40b/SinglePiE200_HO0_GEN_SIM_RECO_1.root' ] )

secFiles.extend( [
    ] )

process.demo = cms.EDAnalyzer('HOWithReco',
                              outfname = \
                              cms.untracked.string("SinglePiE200Reco.root"),
                              mipE = cms.untracked.double(1.0),
                              doFit = cms.untracked.bool(False)
                              )


process.p = cms.Path(process.demo)
