import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load('Configuration/StandardSequences/Services_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration/StandardSequences/MixingNoPileUp_cff')
process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

process.GlobalTag.globaltag = 'MC_31X_V9::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(25000) )

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
process.source = cms.Source ("PoolSource",
                             duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                             fileNames = readFiles,
                             secondaryFileNames = secFiles)
readFiles.extend( [
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_1.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_10.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_11.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_12.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_13.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_14.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_15.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_16.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_17.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_18.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_19.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_2.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_20.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_21.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_22.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_23.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_24.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_25.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_26.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_27.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_28.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_29.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_3.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_30.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_31.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_32.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_33.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_34.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_35.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_36.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_37.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_38.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_39.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_4.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_40.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_41.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_42.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_43.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_44.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_45.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_46.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_47.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_48.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_49.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_5.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_50.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_6.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_7.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_8.root',
'/store/user/andersj/mumin_e_10_60_SiPM_CMSSW_3_3_0/mumin_e_10_60_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_9.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_1.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_10.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_11.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_12.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_13.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_14.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_15.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_16.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_17.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_18.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_19.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_2.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_20.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_21.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_22.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_23.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_24.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_25.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_26.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_27.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_28.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_29.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_3.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_30.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_31.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_32.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_33.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_34.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_35.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_36.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_37.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_38.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_39.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_4.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_40.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_41.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_42.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_43.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_44.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_45.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_46.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_47.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_48.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_49.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_5.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_50.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_6.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_7.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_8.root',
'/store/user/andersj/mumin_e_1_10_SiPM_CMSSW_3_3_0/mumin_e_1_10_SiPM_CMSSW_3_3_0_FullReReco/f522d24cc6a503c1256ec93603a25698/trainingSamples_9.root'
] )


secFiles.extend( [
    ] )

## process.source = cms.Source("PoolSource",
##     # replace 'myfile.root' with the source file you want to use
## #    fileNames = cms.untracked.vstring(
## #        'file:../../HO/SinglePionE300_reco_1.root'
## #    )
##     fileNames = cms.untracked.vstring(
##     #'file:/uscms_data/d2/andersj/HO/SingleMuonE150.root'
## #    'file:../../muminE100Hot.root'
##     'file:/uscms_data/d2/andersj/mu/muminE100Hot.root'
##     )
## )

process.demo = cms.EDAnalyzer('MuonTemplateProducer',
                 outfname = cms.untracked.string("MuonTemplates.root"),
                 #mipE = cms.untracked.double(4.45)
                 mipE = cms.untracked.double(1.0),
                 centralEta = cms.untracked.int32(8),
                 centralPhi = cms.untracked.int32(1),
                 doMuons = cms.untracked.bool(True),
                              
                 TrackAssociatorParameters = cms.PSet(
    muonMaxDistanceSigmaX = cms.double(0.0),
    muonMaxDistanceSigmaY = cms.double(0.0),
    CSCSegmentCollectionLabel = cms.InputTag("cscSegments"),
    dRHcal = cms.double(9999.0),
    dREcal = cms.double(9999.0),
    CaloTowerCollectionLabel = cms.InputTag("towerMaker"),
    useEcal = cms.bool(True),
    dREcalPreselection = cms.double(0.05),
    HORecHitCollectionLabel = cms.InputTag("horeco"),
    dRMuon = cms.double(9999.0),
    trajectoryUncertaintyTolerance = cms.double(-1.0),
    crossedEnergyType = cms.string('SinglePointAlongTrajectory'),
    muonMaxDistanceX = cms.double(5.0),
    muonMaxDistanceY = cms.double(5.0),
    useHO = cms.bool(True),
    accountForTrajectoryChangeCalo = cms.bool(False),
    DTRecSegment4DCollectionLabel = cms.InputTag("dt4DSegments"),
    EERecHitCollectionLabel = cms.InputTag("ecalRecHit",
                                           "EcalRecHitsEE"),
    dRHcalPreselection = cms.double(0.2),
    useMuon = cms.bool(True),
    useCalo = cms.bool(True),
    EBRecHitCollectionLabel = cms.InputTag("ecalRecHit",
                                           "EcalRecHitsEB"),
    dRMuonPreselection = cms.double(0.2),
    usePreshower = cms.bool(False),
    dRPreshowerPreselection = cms.double(0.2),
    truthMatch = cms.bool(False),
    HBHERecHitCollectionLabel = cms.InputTag("hbhereco"),
    useHcal = cms.bool(True)
    )

)


process.p = cms.Path(process.demo)
