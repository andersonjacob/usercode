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
    #'file:/uscms_data/d2/andersj/HO/SingleMuonE150.root'
#    'file:../../muminE100Hot.root'
    'file:/uscms_data/d2/andersj/mu/mumin_e_10_60_hacked.root'
## '/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/6cd5d5c7f0880963f9ca2f4fdf44d555/pimin1000_9.root',
## '/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/6cd5d5c7f0880963f9ca2f4fdf44d555/pimin1000_8.root',
## '/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/6cd5d5c7f0880963f9ca2f4fdf44d555/pimin1000_7.root',
## '/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/6cd5d5c7f0880963f9ca2f4fdf44d555/pimin1000_6.root',
## '/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/6cd5d5c7f0880963f9ca2f4fdf44d555/pimin1000_5.root',
## '/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/6cd5d5c7f0880963f9ca2f4fdf44d555/pimin1000_4.root',
## '/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/6cd5d5c7f0880963f9ca2f4fdf44d555/pimin1000_3.root',
## '/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/6cd5d5c7f0880963f9ca2f4fdf44d555/pimin1000_25.root',
## '/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/6cd5d5c7f0880963f9ca2f4fdf44d555/pimin1000_24.root',
## '/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/6cd5d5c7f0880963f9ca2f4fdf44d555/pimin1000_23.root',
## '/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/6cd5d5c7f0880963f9ca2f4fdf44d555/pimin1000_22.root',
## '/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/6cd5d5c7f0880963f9ca2f4fdf44d555/pimin1000_21.root',
## '/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/6cd5d5c7f0880963f9ca2f4fdf44d555/pimin1000_20.root',
## '/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/6cd5d5c7f0880963f9ca2f4fdf44d555/pimin1000_2.root',
## '/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/6cd5d5c7f0880963f9ca2f4fdf44d555/pimin1000_19.root',
## '/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/6cd5d5c7f0880963f9ca2f4fdf44d555/pimin1000_18.root',
## '/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/6cd5d5c7f0880963f9ca2f4fdf44d555/pimin1000_17.root',
## '/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/6cd5d5c7f0880963f9ca2f4fdf44d555/pimin1000_16.root',
## '/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/6cd5d5c7f0880963f9ca2f4fdf44d555/pimin1000_15.root',
## '/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/6cd5d5c7f0880963f9ca2f4fdf44d555/pimin1000_14.root',
## '/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/6cd5d5c7f0880963f9ca2f4fdf44d555/pimin1000_13.root',
## '/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/6cd5d5c7f0880963f9ca2f4fdf44d555/pimin1000_12.root',
## '/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/6cd5d5c7f0880963f9ca2f4fdf44d555/pimin1000_11.root',
## '/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/6cd5d5c7f0880963f9ca2f4fdf44d555/pimin1000_10.root',
## '/store/user/andersj/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/TeVpiMinus_HOSiPM_LowGain_5fCMIP_CMSSW_3_3_0/6cd5d5c7f0880963f9ca2f4fdf44d555/pimin1000_1.root'

##     '/store/user/andersj/SingleJetPt1000_ZecotekSiPMHO/SingleJetPt1000_ZecotekSiPMHO_Reco/1a4de225bfc658a405ac459d9e14eead/SingleJetPt1000_includeHO_reco_1.root',
##     '/store/user/andersj/SingleJetPt1000_ZecotekSiPMHO/SingleJetPt1000_ZecotekSiPMHO_Reco/1a4de225bfc658a405ac459d9e14eead/SingleJetPt1000_includeHO_reco_2.root',
##     '/store/user/andersj/SingleJetPt1000_ZecotekSiPMHO/SingleJetPt1000_ZecotekSiPMHO_Reco/1a4de225bfc658a405ac459d9e14eead/SingleJetPt1000_includeHO_reco_3.root',
##     '/store/user/andersj/SingleJetPt1000_ZecotekSiPMHO/SingleJetPt1000_ZecotekSiPMHO_Reco/1a4de225bfc658a405ac459d9e14eead/SingleJetPt1000_includeHO_reco_4.root'

##    '/store/user/andersj/SinglePionE300_HamamatsuSiPMHO_HcalOnly/SinglePionE300_HamamatsuSiPMHO_HcalOnly_Reco/87d375cd16331b810d688abd73696634/SinglePionE300_reco_2.root',
##    '/store/user/andersj/SinglePionE300_HamamatsuSiPMHO_HcalOnly/SinglePionE300_HamamatsuSiPMHO_HcalOnly_Reco/87d375cd16331b810d688abd73696634/SinglePionE300_reco_4.root',
##    '/store/user/andersj/SinglePionE300_HamamatsuSiPMHO_HcalOnly/SinglePionE300_HamamatsuSiPMHO_HcalOnly_Reco/87d375cd16331b810d688abd73696634/SinglePionE300_reco_6.root',
##    '/store/user/andersj/SinglePionE300_HamamatsuSiPMHO_HcalOnly/SinglePionE300_HamamatsuSiPMHO_HcalOnly_Reco/87d375cd16331b810d688abd73696634/SinglePionE300_reco_8.root',
##     '/store/user/andersj/SinglePionE300_HamamatsuSiPMHO_HcalOnly/SinglePionE300_HamamatsuSiPMHO_HcalOnly_Reco/87d375cd16331b810d688abd73696634/SinglePionE300_reco_10.root',
##     '/store/user/andersj/SinglePionE300_HamamatsuSiPMHO_HcalOnly/SinglePionE300_HamamatsuSiPMHO_HcalOnly_Reco/87d375cd16331b810d688abd73696634/SinglePionE300_reco_12.root',
##     '/store/user/andersj/SinglePionE300_HamamatsuSiPMHO_HcalOnly/SinglePionE300_HamamatsuSiPMHO_HcalOnly_Reco/87d375cd16331b810d688abd73696634/SinglePionE300_reco_14.root',
##     '/store/user/andersj/SinglePionE300_HamamatsuSiPMHO_HcalOnly/SinglePionE300_HamamatsuSiPMHO_HcalOnly_Reco/87d375cd16331b810d688abd73696634/SinglePionE300_reco_16.root',
##     '/store/user/andersj/SinglePionE300_HamamatsuSiPMHO_HcalOnly/SinglePionE300_HamamatsuSiPMHO_HcalOnly_Reco/87d375cd16331b810d688abd73696634/SinglePionE300_reco_18.root',
##     '/store/user/andersj/SinglePionE300_HamamatsuSiPMHO_HcalOnly/SinglePionE300_HamamatsuSiPMHO_HcalOnly_Reco/87d375cd16331b810d688abd73696634/SinglePionE300_reco_20.root',
##     '/store/user/andersj/SinglePionE300_HamamatsuSiPMHO_HcalOnly/SinglePionE300_HamamatsuSiPMHO_HcalOnly_Reco/87d375cd16331b810d688abd73696634/SinglePionE300_reco_1.root',
##     '/store/user/andersj/SinglePionE300_HamamatsuSiPMHO_HcalOnly/SinglePionE300_HamamatsuSiPMHO_HcalOnly_Reco/87d375cd16331b810d688abd73696634/SinglePionE300_reco_3.root',
##     '/store/user/andersj/SinglePionE300_HamamatsuSiPMHO_HcalOnly/SinglePionE300_HamamatsuSiPMHO_HcalOnly_Reco/87d375cd16331b810d688abd73696634/SinglePionE300_reco_5.root',
##     '/store/user/andersj/SinglePionE300_HamamatsuSiPMHO_HcalOnly/SinglePionE300_HamamatsuSiPMHO_HcalOnly_Reco/87d375cd16331b810d688abd73696634/SinglePionE300_reco_7.root',
##     '/store/user/andersj/SinglePionE300_HamamatsuSiPMHO_HcalOnly/SinglePionE300_HamamatsuSiPMHO_HcalOnly_Reco/87d375cd16331b810d688abd73696634/SinglePionE300_reco_9.root',
##     '/store/user/andersj/SinglePionE300_HamamatsuSiPMHO_HcalOnly/SinglePionE300_HamamatsuSiPMHO_HcalOnly_Reco/87d375cd16331b810d688abd73696634/SinglePionE300_reco_11.root',
##     '/store/user/andersj/SinglePionE300_HamamatsuSiPMHO_HcalOnly/SinglePionE300_HamamatsuSiPMHO_HcalOnly_Reco/87d375cd16331b810d688abd73696634/SinglePionE300_reco_13.root',
##     '/store/user/andersj/SinglePionE300_HamamatsuSiPMHO_HcalOnly/SinglePionE300_HamamatsuSiPMHO_HcalOnly_Reco/87d375cd16331b810d688abd73696634/SinglePionE300_reco_15.root',
##     '/store/user/andersj/SinglePionE300_HamamatsuSiPMHO_HcalOnly/SinglePionE300_HamamatsuSiPMHO_HcalOnly_Reco/87d375cd16331b810d688abd73696634/SinglePionE300_reco_17.root',
##     '/store/user/andersj/SinglePionE300_HamamatsuSiPMHO_HcalOnly/SinglePionE300_HamamatsuSiPMHO_HcalOnly_Reco/87d375cd16331b810d688abd73696634/SinglePionE300_reco_19.root'
    )
)

process.demo = cms.EDAnalyzer('HOSiPMAnalysis',
                 outfname = cms.untracked.string("SinglePion1000.root"),
                 #mipE = cms.untracked.double(4.45)
                 mipE = cms.untracked.double(1.0),
                 centralEta = cms.untracked.int32(8),
                 centralPhi = cms.untracked.int32(2)
                              
)


process.p = cms.Path(process.demo)
