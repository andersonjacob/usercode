import FWCore.ParameterSet.Config as cms
import sys

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = 500

process.load("andersj.HcalPlotter.hcal_tb10_cff")

# We are getting the runnumber from the command line now. No need to edit this file.
#print "Enter Runnumber: "
#rn = sys.stdin.readline()
RUNNUMBER = 2264
print "Running on: " + str(RUNNUMBER)

# Use TB source for HTB data fromat
#process.source = cms.Source("HcalTBSource",
#                            fileNames = cms.untracked.vstring( 
#"file:/data/spool/HTB_00" + str(RUNNUMBER) + ".root"                            
#)
#)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    "file:/tmp_mnt/data0/spool/EcalHcalCombined2010_0000" + str(RUNNUMBER) + ".0.root",
    )
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('test/tb2_0000' + str(RUNNUMBER) +'.root')
)
process.load("EventFilter.EcalRawToDigi.EcalUnpackerData_cfi")
process.load("EventFilter.EcalRawToDigi.EcalUnpackerMapping_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("Geometry.CaloEventSetup.EcalTrigTowerConstituents_cfi")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.load("Geometry.EcalMapping.EcalMapping_cfi")
process.load("Geometry.EcalMapping.EcalMappingRecord_cfi")
process.load("CalibCalorimetry.EcalLaserCorrection.ecalLaserCorrectionService_cfi")
process.load("CalibCalorimetry.EcalTrivialCondModules.EcalTrivialCondRetriever_cfi")

process.ecalEBunpacker.FEDs = cms.vint32(610)
process.ecalEBunpacker.memUnpacking = cms.bool(False)
process.ecalEBunpacker.srpUnpacking = cms.bool(False)

process.load("ecalLocalRecoSequence_cff")
## #uncalibrated rechits
## import RecoLocalCalo.EcalRecProducers.ecalFixedAlphaBetaFitUncalibRecHit_cfi
## process.ecalUncalibHit = RecoLocalCalo.EcalRecProducers.ecalFixedAlphaBetaFitUncalibRecHit_cfi.ecalFixedAlphaBetaFitUncalibRecHit.clone()

## #rechits
## process.load("RecoLocalCalo.EcalRecProducers.ecalRecHit_cfi")

## process.ecalUncalibHit.MinAmplBarrel = 12.
## process.ecalUncalibHit.MinAmplEndcap = 16.
## process.ecalUncalibHit.EBdigiCollection = 'ecalEBunpacker:ebDigis'
## process.ecalUncalibHit.EEdigiCollection = 'ecalEBunpacker:eeDigis'

## process.ecalRecHit.EBuncalibRecHitCollection = 'ecalUncalibHit:EcalUncalibRecHitsEB'
## process.ecalRecHit.EEuncalibRecHitCollection = 'ecalUncalibHit:EcalUncalibRecHitsEE'

## process.ecalDataSequence = cms.Sequence(process.ecalEBunpacker*process.ecalUncalibHit*process.ecalRecHit)

## process.ecalUncalibRecHit = cms.EDProducer(
##     "EcalUncalibRecHitProducer",
##     EBdigiCollection = cms.InputTag("ecalEBunpacker","ebDigis"),
##     EEhitCollection = cms.string('EcalUncalibRecHitsEE'),
##     EEdigiCollection = cms.InputTag("ecalEBunpacker","eeDigis"),
##     EBhitCollection = cms.string('EcalUncalibRecHitsEB'),
##     algo = cms.string("EcalUncalibRecHitWorkerWeights")
##     )

## process.ecalRecHit = cms.EDProducer(
##     "EcalRecHitProducer",
##     EErechitCollection = cms.string(''),
##     EEuncalibRecHitCollection = cms.InputTag(''),
##     EBuncalibRecHitCollection = cms.InputTag('ecalUncalibRecHit', 'EcalUncalibRecHitsEB'),
##     EBrechitCollection = cms.string('EcalRecHitsEB'),
##     algo = cms.string('EcalRecHitWorkerSimple')
##     )

process.load("EventFilter.HcalRawToDigi.HcalRawToDigi_cfi")
process.load("RecoLocalCalo.HcalRecProducers.HcalSimpleReconstructor_hbhe_cfi")
process.load("RecoLocalCalo.HcalRecProducers.HcalSimpleReconstructor_ho_cfi")
process.load("RecoLocalCalo.HcalRecProducers.HcalSimpleReconstructor_hf_cfi")
process.load("RecoTBCalo.HcalTBObjectUnpacker.HcalTBObjectUnpacker_Normal_cfi")

process.horeco.firstSample=2
process.horeco.samplesToAdd=6

process.hfreco.firstSample=2
process.hfreco.samplesToAdd=9

## process.plotanal=cms.EDAnalyzer(
##     "HcalHOTBPlotAnal",
##     hbheRHtag = cms.untracked.InputTag("hbhereco"),
##     hoRHtag   = cms.untracked.InputTag("horeco"),
##     hfRHtag   = cms.untracked.InputTag("hfreco"),
##     hcalDigiTag = cms.untracked.InputTag("hcalDigis"),
##     hcalTrigTag = cms.untracked.InputTag("tbunpack"),
##     calibFC2GeV = cms.untracked.double(0.01),
##     HistoParameters = cms.PSet(
##         pedGeVlo   = cms.double(-15),
##         pedGeVhi   = cms.double(15),
##         pedADClo   = cms.double(0),
##         pedADChi   = cms.double(49),
##         ledGeVlo   = cms.double(-5),
##         ledGeVhi   = cms.double(250),
##         laserGeVlo = cms.double(-5),
##         laserGeVhi = cms.double(350),
##         otherGeVlo = cms.double(-5),
##         otherGeVhi = cms.double(250),
##         beamGeVlo  = cms.double(-20),
##         beamGeVhi  = cms.double(50),
##         #beamGeVhi  = cms.double(80*200),
##         timeNSlo   = cms.double(50),
##         timeNShi   = cms.double(250)
##         )
##     )

process.dump = cms.EDAnalyzer("EventContentAnalyzer")

#process.p=cms.Path(process.ecalEBunpacker+process.hcalDigis+process.hbhereco+process.horeco+process.hfreco+process.tbunpack+process.dump)
process.p=cms.Path(process.ecalEBunpacker+process.ecalGlobalUncalibRecHit+process.ecalRecHit+process.hcalDigis+process.hbhereco+process.horeco+process.hfreco+process.tbunpack+process.dump)
#process.p=cms.Path(process.ecalEBunpacker+process.ecalLocalRecoSequence+process.hcalDigis+process.hbhereco+process.horeco+process.hfreco+process.tbunpack+process.dump)
