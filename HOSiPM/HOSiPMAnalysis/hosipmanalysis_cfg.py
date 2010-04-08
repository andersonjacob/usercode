import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
from TrackingTools.TrackAssociator.default_cfi import TrackAssociatorParameterBlock

process.GlobalTag.globaltag = 'MC_3XY_V21::All'

process.load("RecoMuon.MuonIdentification.links_cfi")


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 

process.source = cms.Source(
    "PoolSource", fileNames = readFiles,
    secondaryFileNames = secFiles,
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
    )

readFiles.extend( [
    '/store/user/andersj/SingleMum_e_10_60_GEN_SIM_RECO_50fCperMIP/SingleMum_e_10_60_GEN_SIM_RECO_50fCperMIP/f229d1db935df3676e6e05502ae8f563/SingleMum_e_10_60_GEN_SIM_RECO_1.root',
    '/store/user/andersj/SingleMum_e_10_60_GEN_SIM_RECO_50fCperMIP/SingleMum_e_10_60_GEN_SIM_RECO_50fCperMIP/f229d1db935df3676e6e05502ae8f563/SingleMum_e_10_60_GEN_SIM_RECO_2.root',
    '/store/user/andersj/SingleMum_e_10_60_GEN_SIM_RECO_50fCperMIP/SingleMum_e_10_60_GEN_SIM_RECO_50fCperMIP/f229d1db935df3676e6e05502ae8f563/SingleMum_e_10_60_GEN_SIM_RECO_3.root',
    '/store/user/andersj/SingleMum_e_10_60_GEN_SIM_RECO_50fCperMIP/SingleMum_e_10_60_GEN_SIM_RECO_50fCperMIP/f229d1db935df3676e6e05502ae8f563/SingleMum_e_10_60_GEN_SIM_RECO_4.root',
    '/store/user/andersj/SingleMum_e_10_60_GEN_SIM_RECO_50fCperMIP/SingleMum_e_10_60_GEN_SIM_RECO_50fCperMIP/f229d1db935df3676e6e05502ae8f563/SingleMum_e_10_60_GEN_SIM_RECO_5.root',
    '/store/user/andersj/SingleMum_e_10_60_GEN_SIM_RECO_50fCperMIP/SingleMum_e_10_60_GEN_SIM_RECO_50fCperMIP/f229d1db935df3676e6e05502ae8f563/SingleMum_e_10_60_GEN_SIM_RECO_6.root',
    '/store/user/andersj/SingleMum_e_10_60_GEN_SIM_RECO_50fCperMIP/SingleMum_e_10_60_GEN_SIM_RECO_50fCperMIP/f229d1db935df3676e6e05502ae8f563/SingleMum_e_10_60_GEN_SIM_RECO_7.root',
    '/store/user/andersj/SingleMum_e_10_60_GEN_SIM_RECO_50fCperMIP/SingleMum_e_10_60_GEN_SIM_RECO_50fCperMIP/f229d1db935df3676e6e05502ae8f563/SingleMum_e_10_60_GEN_SIM_RECO_8.root',
    '/store/user/andersj/SingleMum_e_10_60_GEN_SIM_RECO_50fCperMIP/SingleMum_e_10_60_GEN_SIM_RECO_50fCperMIP/f229d1db935df3676e6e05502ae8f563/SingleMum_e_10_60_GEN_SIM_RECO_9.root',
    '/store/user/andersj/SingleMum_e_10_60_GEN_SIM_RECO_50fCperMIP/SingleMum_e_10_60_GEN_SIM_RECO_50fCperMIP/f229d1db935df3676e6e05502ae8f563/SingleMum_e_10_60_GEN_SIM_RECO_10.root'
   ] )

secFiles.extend( [
    ] )

process.muons.inputCollectionLabels = cms.VInputTag(
    cms.InputTag("generalTracks"),
    cms.InputTag("globalMuonLinks"), 
    cms.InputTag("standAloneMuons","UpdatedAtVtx"))
process.muons.inputCollectionTypes = cms.vstring('inner tracks', 
						 'links', 
						 'outer tracks')

process.muons.minCaloCompatibility = 0.
process.muons.minPt = 0.
process.muons.minP = 0.
process.muons.minNumberOfMatches = 0
#process.muons.MuonCaloCompatibility.MuonTemplateFileName = 'RecoMuon/MuonIdentification/data/MuID_templates_muons_barrel_correctedHO_SiPMs.root';
#process.muons.MuonCaloCompatibility.PionTemplateFileName = 'RecoMuon/MuonIdentification/data/MuID_templates_pions_barrel_correctedHO_SiPMs.root';

process.MuonNumberingInitialization = cms.ESProducer("MuonNumberingInitialization")

process.demo = cms.EDAnalyzer(
    'HOSiPMAnalysis',
    outfname = cms.untracked.string("mum_e_10_60_test.root"),
    #mipE = cms.untracked.double(4.45)
    mipE = cms.untracked.double(1.0),
    delta_eta = cms.untracked.double(0.02),
    delta_phi = cms.untracked.double(0.02),
    doMuons = cms.untracked.bool(True),
    TrackAssociatorParameters = TrackAssociatorParameterBlock.TrackAssociatorParameters
)


process.p = cms.Path(process.globalMuonLinks*process.muons*process.demo)
