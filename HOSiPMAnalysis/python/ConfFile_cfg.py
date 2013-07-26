import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.Geometry.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("EmptySource",
                            firstRun = cms.untracked.uint32(1),
                            numberEventsInRun = cms.untracked.uint32(1)
                            )

process.demo = cms.EDAnalyzer(
    'PedTextWriter',
    hcalDigiLabel = cms.untracked.InputTag('simHcalUnsuppressedDigis')
)


process.p = cms.Path(process.demo)
