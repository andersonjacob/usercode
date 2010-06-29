import FWCore.ParameterSet.Config as cms

# HCAL setup suitable for MC simulation and production
hcal_db_producer = cms.ESProducer("HcalDbProducer",
    dump = cms.untracked.vstring(''),
    file = cms.untracked.string('')
)

es_hardcode = cms.ESSource("HcalHardcodeCalibrations",
    toGet=cms.untracked.vstring("Pedestals","PedestalWidths","Gains","QIEData",
                                "GainWidths","ElectronicsMap",
                                "ChannelQuality","RespCorrs","ZSThresholds"),
                           H2Mode=cms.untracked.bool(False),
                           SLHCMode=cms.untracked.bool(True)
)


