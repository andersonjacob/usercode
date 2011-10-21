import FWCore.ParameterSet.Config as cms
import sys
import os
import re

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = 500
process.MessageLogger.suppressWarning.extend(["ecalRecHit"])

process.load("andersj.HcalPlotter.hcal_tb11_cff")

# We are getting the runnumber from the command line now. No need to edit this file.
if len(sys.argv) > 2:
    #print sys.argv
    arg1 = sys.argv[2]
    RUNNUMBER = int(arg1)
    inputieta = int(sys.argv[3])
    inputiphi = int(sys.argv[4])
    #print arg1
else:
    RUNNUMBER = 0
if (RUNNUMBER < 1):
    rn = input("Enter Runnumber: ")
    RUNNUMBER = int(rn)
    tmp = input("Enter hb table ieta: ")
    inputieta = int(tmp)
    tmp = input("Enter hb table iphi: ")
    inputiphi = int(tmp)
print "Running on: {0:08d} table coordinates ({1},{2})".format(RUNNUMBER,
                                                               inputieta,
                                                               inputiphi)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

datadir = 'tb_data'

files = ['file:{0}/{1}'.format(datadir,x) for x in os.listdir(datadir) \
         if re.search('_{0:08d}\.\d+'.format(RUNNUMBER), x)]

#print files

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(files
    )
)

#process.source.fileNames.extend(files)

print process.source.fileNames

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('test/tb2011_{0:08d}.root'.format(RUNNUMBER))
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
process.ecalEBunpacker.syncCheck = cms.bool(False)
process.ecalEBunpacker.orderedFedList = cms.vint32(610)
process.ecalEBunpacker.silentMode = cms.untracked.bool(False)

import RecoLocalCalo.EcalRecProducers.ecalRatioUncalibRecHit_cfi
process.ecalUncalibHit = RecoLocalCalo.EcalRecProducers.ecalRatioUncalibRecHit_cfi.ecalRatioUncalibRecHit.clone()

## import RecoLocalCalo.EcalRecProducers.ecalFixedAlphaBetaFitUncalibRecHit_cfi

## process.ecalUncalibHit = RecoLocalCalo.EcalRecProducers.ecalFixedAlphaBetaFitUncalibRecHit_cfi.ecalFixedAlphaBetaFitUncalibRecHit.clone()

process.ecalUncalibHit.EBdigiCollection = 'ecalEBunpacker:ebDigis'
process.ecalUncalibHit.EEdigiCollection = 'ecalEBunpacker:eeDigis'

process.load("RecoLocalCalo.EcalRecProducers.ecalRecHit_cfi")
process.ecalRecHit.ebDetIdToBeRecovered = cms.InputTag("","")
process.ecalRecHit.eeDetIdToBeRecovered = cms.InputTag("","")
process.ecalRecHit.eeFEToBeRecovered = cms.InputTag("","")
process.ecalRecHit.ebFEToBeRecovered = cms.InputTag("","")
process.ecalRecHit.recoverEBFE = cms.bool(False)
process.ecalRecHit.recoverEEFE = cms.bool(False)
process.ecalRecHit.EBuncalibRecHitCollection = 'ecalUncalibHit:EcalUncalibRecHitsEB'
process.ecalRecHit.EEuncalibRecHitCollection = 'ecalUncalibHit:EcalUncalibRecHitsEE'

process.ecalDataSequence = cms.Sequence(process.ecalEBunpacker*process.ecalUncalibHit*process.ecalRecHit)

process.load("EventFilter.HcalRawToDigi.HcalRawToDigi_cfi")
process.load("RecoLocalCalo.HcalRecProducers.HcalSimpleReconstructor_hbhe_cfi")
process.load("RecoLocalCalo.HcalRecProducers.HcalSimpleReconstructor_ho_cfi")
#process.load("RecoLocalCalo.HcalRecProducers.HcalSimpleReconstructor_hf_cfi")
process.load("andersj.HcalPlotter.HcalTBObjectUnpacker_Normal_cfi")

## process.horeco.firstSample=2
## process.horeco.samplesToAdd=6
process.horeco.firstSample=3
process.horeco.samplesToAdd=4
process.horeco.correctForPhaseContainment = False
process.horeco.correctForTimeslew = False

#process.hfreco.firstSample=2
#process.hfreco.samplesToAdd=9

process.plotanal=cms.EDAnalyzer(
    "HcalHOTBPlotAnal",
    hbheRHtag = cms.untracked.InputTag("hbheprereco"),
    hoRHtag   = cms.untracked.InputTag("horeco"),
    hfRHtag   = cms.untracked.InputTag("hfreco"),
    hcalDigiTag = cms.untracked.InputTag("hcalDigis"),
    hcalTrigTag = cms.untracked.InputTag("tbunpack"),
    ebRHtag = cms.untracked.InputTag("ecalRecHit:EcalRecHitsEB"),
    doBeamCounters = cms.untracked.bool(True),
    calibFC2GeV = cms.untracked.double(0.01),
    HistoParameters = cms.PSet(
        pedGeVlo   = cms.double(-15),
        pedGeVhi   = cms.double(15),
        pedADClo   = cms.double(0),
        pedADChi   = cms.double(49),
        ledGeVlo   = cms.double(-5),
        ledGeVhi   = cms.double(250),
        laserGeVlo = cms.double(-5),
        laserGeVhi = cms.double(350),
        otherGeVlo = cms.double(-5),
        otherGeVhi = cms.double(250),
        beamGeVlo  = cms.double(-20),
        beamGeVhi  = cms.double(500),
        #beamGeVhi  = cms.double(80*200),
        timeNSlo   = cms.double(50),
        timeNShi   = cms.double(250),
        ieta = cms.int32(inputieta),
        iphi = cms.int32(inputiphi)
        )
    )

##process.dump = cms.EDAnalyzer("EventContentAnalyzer")

#process.p=cms.Path(process.ecalDataSequence+process.hcalDigis+process.hbheprereco+process.horeco+process.hfreco+process.tbunpack+process.plotanal)
process.p=cms.Path(process.ecalDataSequence+process.hcalDigis+process.hbheprereco+process.horeco+process.tbunpack+process.plotanal)
