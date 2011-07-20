import FWCore.ParameterSet.Config as cms

hcal_db_producer = cms.ESProducer("HcalDbProducer",
     dump = cms.untracked.vstring(''),
     file = cms.untracked.string('')
)
 
hcales_ascii = cms.ESSource("HcalTextCalibrations",
                        input = cms.VPSet(
    cms.PSet(
    object = cms.string('ElectronicsMap'),
    #file = cms.FileInPath('andersj/HcalPlotter/data/tb2010_map_HO_HPD.txt')
    file = cms.FileInPath('andersj/HcalPlotter/data/tb2011_map.txt')
    ),
    cms.PSet(
    object = cms.string('Pedestals'),
    file = cms.FileInPath('andersj/HcalPlotter/data/ped_tb2010_000495.txt')
    #file = cms.FileInPath('andersj/HcalPlotter/data/ped_tb2010_000495_HO_HPD.txt')
    ),
    cms.PSet(
    object = cms.string('Gains'),
    file = cms.FileInPath('andersj/HcalPlotter/data/gain_tb2010_muon1.txt')
    )
    )
                                )

hcales_hardcode = cms.ESSource("HcalHardcodeCalibrations",
                               toGet = cms.untracked.vstring('PedestalWidths', 'LutMetadata',
                                                             'GainWidths', 'LUTCorrs', 
                                                             'PFCorrs', 'QIEData',
                                                             'L1TriggerObjects','ZSThresholds','DcsValues',
                                                             'ChannelQuality','RespCorrs','TimeCorrs')
                               ) 
