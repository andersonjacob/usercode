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
#    'file:/uscms_data/d2/andersj/HO/SingleMuonE150_reco.root'
#    'file:../../muminE100Hot.root'
'/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/baf6984ab4a8bc85ba6a0d6aa70509bb/pimin1000_9.root',
'/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/baf6984ab4a8bc85ba6a0d6aa70509bb/pimin1000_8.root',
'/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/baf6984ab4a8bc85ba6a0d6aa70509bb/pimin1000_7.root',
'/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/baf6984ab4a8bc85ba6a0d6aa70509bb/pimin1000_6.root',
'/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/baf6984ab4a8bc85ba6a0d6aa70509bb/pimin1000_5.root',
'/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/baf6984ab4a8bc85ba6a0d6aa70509bb/pimin1000_4.root',
'/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/baf6984ab4a8bc85ba6a0d6aa70509bb/pimin1000_3.root',
'/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/baf6984ab4a8bc85ba6a0d6aa70509bb/pimin1000_25.root',
'/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/baf6984ab4a8bc85ba6a0d6aa70509bb/pimin1000_24.root',
'/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/baf6984ab4a8bc85ba6a0d6aa70509bb/pimin1000_23.root',
'/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/baf6984ab4a8bc85ba6a0d6aa70509bb/pimin1000_22.root',
'/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/baf6984ab4a8bc85ba6a0d6aa70509bb/pimin1000_21.root',
'/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/baf6984ab4a8bc85ba6a0d6aa70509bb/pimin1000_20.root',
'/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/baf6984ab4a8bc85ba6a0d6aa70509bb/pimin1000_2.root',
'/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/baf6984ab4a8bc85ba6a0d6aa70509bb/pimin1000_19.root',
'/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/baf6984ab4a8bc85ba6a0d6aa70509bb/pimin1000_18.root',
'/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/baf6984ab4a8bc85ba6a0d6aa70509bb/pimin1000_17.root',
'/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/baf6984ab4a8bc85ba6a0d6aa70509bb/pimin1000_16.root',
'/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/baf6984ab4a8bc85ba6a0d6aa70509bb/pimin1000_15.root',
'/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/baf6984ab4a8bc85ba6a0d6aa70509bb/pimin1000_14.root',
'/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/baf6984ab4a8bc85ba6a0d6aa70509bb/pimin1000_13.root',
'/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/baf6984ab4a8bc85ba6a0d6aa70509bb/pimin1000_12.root',
'/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/baf6984ab4a8bc85ba6a0d6aa70509bb/pimin1000_11.root',
'/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/baf6984ab4a8bc85ba6a0d6aa70509bb/pimin1000_10.root',
'/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta13/baf6984ab4a8bc85ba6a0d6aa70509bb/pimin1000_1.root'
    )
)

process.demo = cms.EDAnalyzer('HOSiPMAnalysis',
                 outfname = cms.untracked.string("SinglePion1000Eta13.root"),
                 #mipE = cms.untracked.double(4.45)
                 mipE = cms.untracked.double(1.0),
                 centralEta = cms.untracked.int32(13),
                 centralPhi = cms.untracked.int32(1)
                              
)


process.p = cms.Path(process.demo)
