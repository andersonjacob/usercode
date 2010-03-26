import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
process.source = cms.Source ("PoolSource",fileNames = readFiles,
                             secondaryFileNames = secFiles)

readFiles.extend( [
       '/store/user/andersj/SingleMuE50_AllHcalLayers_GEN_SIM/SingleMuE50_AllHcalLayers_GEN_SIM/9098ce6b71f4a04080d7fe88db3bfe5d/SingleMuE50_HO0_GEN_SIM_RECO_2.root',
       '/store/user/andersj/SingleMuE50_AllHcalLayers_GEN_SIM/SingleMuE50_AllHcalLayers_GEN_SIM/9098ce6b71f4a04080d7fe88db3bfe5d/SingleMuE50_HO0_GEN_SIM_RECO_1.root' ] );

secFiles.extend( [
    ] )

process.demo = cms.EDAnalyzer('HOWithReco',
                              outfname = \
                              cms.untracked.string("SingleMuE50Reco.root"),
                              mipE = cms.untracked.double(1.0),
                              doFit = cms.untracked.bool(False)
                              )


process.p = cms.Path(process.demo)
