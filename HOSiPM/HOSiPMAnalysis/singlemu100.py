import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
#    fileNames = cms.untracked.vstring(
#        'file:../../HO/SinglePionE300_reco_1.root'
#    )
    fileNames = cms.untracked.vstring(
    #'file:/uscms_data/d2/andersj/HO/SingleMuonE150.root'
#    'file:../../muminE100Hot.root'
    'file:/uscms_data/d2/andersj/mu/muminE100Hot.root'
    )
)

process.demo = cms.EDAnalyzer('HOSiPMAnalysis',
                 outfname = cms.untracked.string("SingleMu100.root"),
                 #mipE = cms.untracked.double(4.45)
                 mipE = cms.untracked.double(1.0),
                 centralEta = cms.untracked.int32(8),
                 centralPhi = cms.untracked.int32(1)
                              
)


process.p = cms.Path(process.demo)
