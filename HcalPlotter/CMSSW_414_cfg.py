import FWCore.ParameterSet.Config as cms
import sys

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = 500
process.MessageLogger.suppressWarning.extend(["ecalRecHit"])

process.load("CaloTestBeam_2011.HcalPlotter.hcal_tb11_cff")

# We are getting the runnumber from the command line now. No need to edit this file.
print "Enter Runnumber: "
rn = sys.stdin.readline()
RUNNUMBER = int(rn.strip())
print "Running on: " + str(RUNNUMBER)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    "file:data/EcalHcalCombined2011_0000" + str(RUNNUMBER) + ".0.root",
    )
)

## process.source = cms.Source("HcalTBSource",
##     fileNames = cms.untracked.vstring(
##     "file:data/HTB_" + str(RUNNUMBER) + ".root",
##     )
## )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('test/tb2_' + str(RUNNUMBER) +'.root')
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

process.load("CaloTestBeam_2011.HcalPlotter.ecalLocalRecoSequence_cff")
process.ecalDataSequence = cms.Sequence(process.ecalEBunpacker*process.ecalGlobalUncalibRecHit*process.ecalRecHit)

process.load("EventFilter.HcalRawToDigi.HcalRawToDigi_cfi")
process.load("RecoLocalCalo.HcalRecProducers.HcalSimpleReconstructor_hbhe_cfi")
process.load("RecoLocalCalo.HcalRecProducers.HcalSimpleReconstructor_ho_cfi")
process.load("RecoLocalCalo.HcalRecProducers.HcalSimpleReconstructor_hf_cfi")
process.load("CaloTestBeam_2011.HcalPlotter.HcalTBObjectUnpacker_Normal_cfi")

process.horeco.firstSample=2
process.horeco.samplesToAdd=6

process.hfreco.firstSample=2
process.hfreco.samplesToAdd=9

#process.hbheprereco.tsFromDB = False  #This is needed for 4_2_2
#process.horeco.tsFromDB = False  #This is needed for 4_2_2
#process.hfreco.tsFromDB = False  #This is needed for 4_2_2

#process.hbheprereco.correctTiming             = False
#process.hbheprereco.setNoiseFlags             = False
#process.hbheprereco.setHSCPFlags              = False
#process.hbheprereco.setSaturationFlags        = False
#process.hbheprereco.setTimingShapedCutsFlags  = False
#process.hbheprereco.setPulseShapeFlags        = False


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
        timeNShi   = cms.double(250)
        )
    )

##process.dump = cms.EDAnalyzer("EventContentAnalyzer")

process.p=cms.Path(process.ecalDataSequence+process.hcalDigis+process.hbheprereco+process.horeco+process.hfreco+process.tbunpack+process.plotanal)