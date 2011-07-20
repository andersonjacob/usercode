import FWCore.ParameterSet.Config as cms

#from RecoLocalCalo.EcalRecAlgos.ecalCleaningAlgo import cleaningAlgoConfig # This is needed in 4_2_2

# rechit producer
ecalRecHit = cms.EDProducer("EcalRecHitProducer",
    EErechitCollection = cms.string('EcalRecHitsEE'),
    EEuncalibRecHitCollection = cms.InputTag("ecalGlobalUncalibRecHit","EcalUncalibRecHitsEE"),
    EBuncalibRecHitCollection = cms.InputTag("ecalGlobalUncalibRecHit","EcalUncalibRecHitsEB"),
    EBrechitCollection = cms.string('EcalRecHitsEB'),
    # channel flags to be exluded from reconstruction, e.g { 1, 2 }
    ChannelStatusToBeExcluded = cms.vint32(),
    # apply laser corrections
    laserCorrection = cms.bool(True),
    # avoid propagation of dead channels other than after recovery
    killDeadChannels = cms.bool(True),
    algo = cms.string("EcalRecHitWorkerSimple"),
    # reco flags association to DB flag
    # the vector index corresponds to the DB flag
    # the value correspond to the reco flag
    flagsMapDBReco = cms.vint32(
             0,   0,   0,  0, # standard reco
             4,               # faulty hardware (noisy)
            -1,  -1,  -1,     # not yet assigned
             4,   4,          # faulty hardware (fixed gain)
             7,   7,   7,     # dead channel with trigger
             8,               # dead FE
             9                # dead or recovery failed
            ),
    # for channel recovery
    algoRecover = cms.string("EcalRecHitWorkerRecover"),
    #algoRecover = cms.string(""),
    recoverEBIsolatedChannels = cms.bool(False),
    recoverEEIsolatedChannels = cms.bool(False),
    recoverEBVFE  = cms.bool(False),
    recoverEEVFE  = cms.bool(False),
    recoverEBFE = cms.bool(False),
    recoverEEFE = cms.bool(False),
    ## ebDetIdToBeRecovered = cms.InputTag("ecalDetIdToBeRecovered:ebDetId"),
    ## eeDetIdToBeRecovered = cms.InputTag("ecalDetIdToBeRecovered:eeDetId"),
    ## ebFEToBeRecovered = cms.InputTag("ecalDetIdToBeRecovered:ebFE"),
    ## eeFEToBeRecovered = cms.InputTag("ecalDetIdToBeRecovered:eeFE"),
    logWarningEtThreshold_EB_FE = cms.double(50),# in EB logWarningThreshold is actually in E (GeV)
    logWarningEtThreshold_EE_FE = cms.double(50),# in EE the energy should correspond to Et (GeV) but the recovered values of energies are not tested if make sense
    ebDetIdToBeRecovered = cms.InputTag(""),
    eeDetIdToBeRecovered = cms.InputTag(""),
    ebFEToBeRecovered = cms.InputTag(""),
    eeFEToBeRecovered = cms.InputTag(""),
    singleChannelRecoveryMethod = cms.string("NeuralNetworks"),
    singleChannelRecoveryThreshold = cms.double(8),
    ## cleaningConfig=cleaningAlgoConfig,  # This is needed in 4_2_2
    triggerPrimitiveDigiCollection = cms.InputTag("ecalEBunpacker:EcalTriggerPrimitives")
)
