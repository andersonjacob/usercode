import FWCore.ParameterSet.Config as cms

from SimCalorimetry.HcalSimProducers.hcalSimParameters_cfi import *
simHcalUnsuppressedDigis = cms.EDProducer("HcalDigiProducer",
    hcalSimParameters,
    doNoise = cms.bool(True),
    doHPDNoise = cms.bool(False),
    doTimeSlew = cms.bool(True),
    hitsProducer = cms.string('g4SimHits'),
    RelabelHits = cms.untracked.bool(True),
    RelabelRules = cms.untracked.PSet(
        # Eta1 = cms.untracked.vint32(1,1,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4),
        Eta1 = cms.untracked.vint32(1,2,2,2,2,3,3,3,3,4,4,4,4,4,4,4,4,5,5),
        Eta17 = cms.untracked.vint32(1,1,1,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4)
        )
)
