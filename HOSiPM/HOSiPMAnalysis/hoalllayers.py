import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
process.source = cms.Source ("PoolSource",fileNames = readFiles,
                             secondaryFileNames = secFiles)

readFiles.extend( [
    '/store/user/andersj/SingleMuE50_AllHcalLayers_GEN_SIM/SingleMuE50_AllHcalLayers_GEN_SIM/171798086f1c004633947652b9c359dd/SingleMuE50_HO0_GEN_SIM_5.root',
    '/store/user/andersj/SingleMuE50_AllHcalLayers_GEN_SIM/SingleMuE50_AllHcalLayers_GEN_SIM/171798086f1c004633947652b9c359dd/SingleMuE50_HO0_GEN_SIM_4.root',
    '/store/user/andersj/SingleMuE50_AllHcalLayers_GEN_SIM/SingleMuE50_AllHcalLayers_GEN_SIM/171798086f1c004633947652b9c359dd/SingleMuE50_HO0_GEN_SIM_3.root',
    '/store/user/andersj/SingleMuE50_AllHcalLayers_GEN_SIM/SingleMuE50_AllHcalLayers_GEN_SIM/171798086f1c004633947652b9c359dd/SingleMuE50_HO0_GEN_SIM_2.root' ] )


secFiles.extend( [
    ] )

process.demo = cms.EDAnalyzer('HOAllLayers',
                              outfname = \
                              cms.untracked.string("SingleMuE50.root"),
                              #mipE = cms.untracked.double(4.45)
                              doFit = cms.untracked.bool(True)
                              )


process.p = cms.Path(process.demo)
