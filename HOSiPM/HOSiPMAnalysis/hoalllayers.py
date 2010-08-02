import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
process.source = cms.Source ("PoolSource",fileNames = readFiles,
                             secondaryFileNames = secFiles)

readFiles.extend( [
    "file:/uscms_data/d1/andersj/upgradeMuons/digi_1_5.root",
    #"file:/uscms_data/d1/andersj/upgradeMuons/digi_2_5.root"
 ] )


secFiles.extend( [
    ] )

process.demo = cms.EDAnalyzer(
    'HOAllLayers',
    outfname = cms.untracked.string("mu100.root"),
    centralEta = cms.untracked.int32(3),
    centralPhi = cms.untracked.int32(1),
    findCenter = cms.untracked.bool(False),
    #mipE = cms.untracked.double(4.45)
    doFit = cms.untracked.bool(False)
    )

process.p = cms.Path(process.demo)
