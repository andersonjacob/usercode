import FWCore.ParameterSet.Config as cms

hcal_db_producer = cms.ESProducer("HcalDbProducer",
     dump = cms.untracked.vstring(''),
     file = cms.untracked.string('')
)
 
hcales_ascii = cms.ESSource(
    "HcalTextCalibrations",
    input = cms.VPSet(
        ## cms.PSet(
        ##     object = cms.string('ElectronicsMap'),
        ##     file = cms.FileInPath('andersj/HcalPlotter/data/tb2011_map_HBHO_HPD.txt')
        ##     ),
        ## cms.PSet(
        ##     object = cms.string('ElectronicsMap'),
        ##     file = cms.FileInPath('andersj/HcalPlotter/data/tb2011_map_HB_HPD_HO_SiPM.txt')
        ##     ),
##         cms.PSet(
##         object = cms.string('ElectronicsMap'),
##         file = cms.FileInPath('andersj/HcalPlotter/data/map_tb2011_HB_ODU_ODU_HO_SiPM.txt')
##         ),
        cms.PSet(
        object = cms.string('ElectronicsMap'),
        file = cms.FileInPath('andersj/HcalPlotter/data/map_tb2011_HB_HPD_ODU_HO_SiPM.txt')
        ),
        ## cms.PSet(
        ## object = cms.string('ElectronicsMap'),
        ## file = cms.FileInPath('andersj/HcalPlotter/data/map_tb2011_HB_EDU_ODU_HO_HPD_SiPM.txt')
        ## ),
        ## cms.PSet(
        ## object = cms.string('Pedestals'),
        ## file = cms.FileInPath('andersj/HcalPlotter/data/ped_tb2011_HBHO_HPD.txt')
        ## ),
        ## cms.PSet(
        ##     object = cms.string('Pedestals'),
        ##     file = cms.FileInPath('andersj/HcalPlotter/data/ped_tb2011_HB_HPD_HO_SiPM.txt')
        ##     ),
        ## cms.PSet(
        ##     object = cms.string('Pedestals'),
        ##     file = cms.FileInPath('andersj/HcalPlotter/data/ped_tb2011_HB_HPD_HO_HPD_SiPM.txt')
        ##     ),
##         cms.PSet(
##             object = cms.string('Pedestals'),
##             file = cms.FileInPath('andersj/HcalPlotter/data/ped_tb2011_HB_ODU_ODU_HO_SiPM.txt')
##             ),
##         cms.PSet(
##             object = cms.string('Pedestals'),
##             file = cms.FileInPath('andersj/HcalPlotter/data/ped_tb2011_HB_HamO_HamO_HO_HPD_SiPM.txt')
##             ),
        cms.PSet(
            object = cms.string('Pedestals'),
            file = cms.FileInPath('andersj/HcalPlotter/data/ped_tb2011_HB_HPD_ODU_HO_SiPM.txt')
            ),
        ## cms.PSet(
        ##     object = cms.string('Pedestals'),
        ##     file = cms.FileInPath('andersj/HcalPlotter/data/ped_tb2011_HB_HamE_HamO_HO_HPD_SiPM.txt')
        ##     ),
##         cms.PSet(
##             object = cms.string('Pedestals'),
##             file = cms.FileInPath('andersj/HcalPlotter/data/ped_tb2011_HB_FBK_HamO_HO_HPD_SiPM.txt')
##             ),
        ## cms.PSet(
        ##     object = cms.string('Gains'),
        ##     file = cms.FileInPath('andersj/HcalPlotter/data/gain_tb2011_HBHO_HPD.txt')
        ##     )
        ## cms.PSet(
        ## object = cms.string('Gains'),
        ## file = cms.FileInPath('andersj/HcalPlotter/data/gain_tb2011_HB_HPD_HO_HPD_SiPM.txt')
        ## )
##         cms.PSet(
##             object = cms.string('Gains'),
##             file = cms.FileInPath('andersj/HcalPlotter/data/gain_tb2011_HB_ODU_ODU_HO_SiPM.txt')
##             )
        cms.PSet(
            object = cms.string('Gains'),
            file = cms.FileInPath('andersj/HcalPlotter/data/gain_tb2011_HB_HPD_ODU_HO_SiPM.txt')
            )
        )
    )

hcales_hardcode = cms.ESSource(
    "HcalHardcodeCalibrations",
    toGet = cms.untracked.vstring(
        'Pedestals',
        'PedestalWidths', 'LutMetadata',
        'Gains',
        'GainWidths', 'LUTCorrs', 
        'PFCorrs', 'QIEData',
        'L1TriggerObjects','ZSThresholds','DcsValues',
        'ChannelQuality','RespCorrs','TimeCorrs',
        'RecoParams')
    )

hcalasciiprefer = cms.ESPrefer("HcalTextCalibrations", "hcales_ascii")
