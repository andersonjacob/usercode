import FWCore.ParameterSet.Config as cms

# producer of rechits starting from digis
ecalGlobalUncalibRecHit = cms.EDProducer("EcalUncalibRecHitProducer",
    EEdigiCollection = cms.InputTag("ecalEBunpacker","eeDigis"),
    EBdigiCollection = cms.InputTag("ecalEBunpacker","ebDigis"),
    EEhitCollection = cms.string("EcalUncalibRecHitsEE"),
    betaEB = cms.double(1.7),
    betaEE = cms.double(1.37),
    AlphaBetaFilename = cms.untracked.string("NOFILE"),
    MinAmplEndcap = cms.double(16.0),
    MinAmplBarrel = cms.double(12.0),
    UseDynamicPedestal = cms.bool(True),
    alphaEB = cms.double(1.2),
    alphaEE = cms.double(1.63),
    EBhitCollection = cms.string("EcalUncalibRecHitsEB"),
    algo = cms.string("EcalUncalibRecHitWorkerFixedAlphaBetaFit")
)
