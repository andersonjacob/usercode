import FWCore.ParameterSet.Config as cms

HcalTopologyIdealEP = cms.ESProducer(
    "HcalTopologyIdealEP",
    SLHCMode=cms.untracked.bool(True)
    )
