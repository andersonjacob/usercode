import FWCore.ParameterSet.Config as cms

ecalDetIdToBeRecovered = cms.EDProducer("EcalDetIdToBeRecoveredProducer",

        # SRP collections
        ebSrFlagCollection = cms.InputTag("ecalDigis"),
        eeSrFlagCollection = cms.InputTag("ecalDigis"),

        # Integrity for xtal data
        ebIntegrityGainErrors = cms.InputTag("ecalEBunpacker:EcalIntegrityGainErrors"),
        ebIntegrityGainSwitchErrors = cms.InputTag("ecalEBunpacker:EcalIntegrityGainSwitchErrors"),
        ebIntegrityChIdErrors = cms.InputTag("ecalEBunpacker:EcalIntegrityChIdErrors"),

        # Integrity for xtal data - EE specific (to be rivisited towards EB+EE common collection)
        eeIntegrityGainErrors = cms.InputTag("ecalEBunpacker:EcalIntegrityGainErrors"),
        eeIntegrityGainSwitchErrors = cms.InputTag("ecalEBunpacker:EcalIntegrityGainSwitchErrors"),
        eeIntegrityChIdErrors = cms.InputTag("ecalEBunpacker:EcalIntegrityChIdErrors"),

        # Integrity Errors
        integrityTTIdErrors = cms.InputTag("ecalEBunpacker:EcalIntegrityTTIdErrors"),
        integrityBlockSizeErrors = cms.InputTag("ecalEBunpacker:EcalIntegrityBlockSizeErrors"),

        # output collections
        ebDetIdToBeRecovered = cms.string("ebDetId"),
        eeDetIdToBeRecovered = cms.string("eeDetId"),
        ebFEToBeRecovered = cms.string("ebFE"),
        eeFEToBeRecovered = cms.string("eeFE")
)
