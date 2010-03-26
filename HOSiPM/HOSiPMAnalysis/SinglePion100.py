import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

readFiles = cms.untracked.vstring(
       '/store/user/andersj/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/118ffc99144cb1d75baba1d815f01ac6/pimin100_9.root',
       '/store/user/andersj/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/118ffc99144cb1d75baba1d815f01ac6/pimin100_8.root',
       '/store/user/andersj/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/118ffc99144cb1d75baba1d815f01ac6/pimin100_7.root',
       '/store/user/andersj/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/118ffc99144cb1d75baba1d815f01ac6/pimin100_6.root',
       '/store/user/andersj/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/118ffc99144cb1d75baba1d815f01ac6/pimin100_5.root',
       '/store/user/andersj/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/118ffc99144cb1d75baba1d815f01ac6/pimin100_4.root',
       '/store/user/andersj/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/118ffc99144cb1d75baba1d815f01ac6/pimin100_3.root',
       '/store/user/andersj/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/118ffc99144cb1d75baba1d815f01ac6/pimin100_25.root',
       '/store/user/andersj/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/118ffc99144cb1d75baba1d815f01ac6/pimin100_24.root',
       '/store/user/andersj/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/118ffc99144cb1d75baba1d815f01ac6/pimin100_23.root',
       '/store/user/andersj/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/118ffc99144cb1d75baba1d815f01ac6/pimin100_22.root',
       '/store/user/andersj/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/118ffc99144cb1d75baba1d815f01ac6/pimin100_21.root',
       '/store/user/andersj/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/118ffc99144cb1d75baba1d815f01ac6/pimin100_20.root',
       '/store/user/andersj/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/118ffc99144cb1d75baba1d815f01ac6/pimin100_2.root',
       '/store/user/andersj/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/118ffc99144cb1d75baba1d815f01ac6/pimin100_19.root',
       '/store/user/andersj/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/118ffc99144cb1d75baba1d815f01ac6/pimin100_18.root',
       '/store/user/andersj/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/118ffc99144cb1d75baba1d815f01ac6/pimin100_17.root',
       '/store/user/andersj/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/118ffc99144cb1d75baba1d815f01ac6/pimin100_16.root',
       '/store/user/andersj/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/118ffc99144cb1d75baba1d815f01ac6/pimin100_15.root',
       '/store/user/andersj/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/118ffc99144cb1d75baba1d815f01ac6/pimin100_14.root',
       '/store/user/andersj/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/118ffc99144cb1d75baba1d815f01ac6/pimin100_13.root',
       '/store/user/andersj/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/118ffc99144cb1d75baba1d815f01ac6/pimin100_12.root',
       '/store/user/andersj/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/118ffc99144cb1d75baba1d815f01ac6/pimin100_11.root',
       '/store/user/andersj/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/118ffc99144cb1d75baba1d815f01ac6/pimin100_10.root',
       '/store/user/andersj/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/piMinus_100GeV_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0_Eta8/118ffc99144cb1d75baba1d815f01ac6/pimin100_1.root'
       );

process.source = cms.Source ("PoolSource",
                             duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                             fileNames = readFiles)

process.demo = cms.EDAnalyzer('HOSiPMAnalysis',
                 outfname = cms.untracked.string("SinglePion100Eta8.root"),
                 #mipE = cms.untracked.double(4.45)
                 mipE = cms.untracked.double(1.0),
                 centralEta = cms.untracked.int32(8),
                 centralPhi = cms.untracked.int32(1)
                              
)

process.p = cms.Path(process.demo)
