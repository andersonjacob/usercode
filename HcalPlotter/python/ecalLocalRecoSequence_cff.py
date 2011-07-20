import FWCore.ParameterSet.Config as cms

# Calo geometry service model
#ECAL conditions
#
# removed : this goes into CalibCalorimetry/Configuration/data/Ecal_FakeCalibrations.cff
#
#  include "CalibCalorimetry/EcalTrivialCondModules/data/EcalTrivialCondRetriever.cfi"
#
#ECAL reconstruction
from ecalUncalibRecHit_cfi import *
from ecalRecHit_cfi import *
from ecalDetIdToBeRecovered_cfi import *
from ecalFixedAlphaBetaFitUncalibRecHit_cfi import *
ecalLocalRecoSequence = cms.Sequence(ecalUncalibRecHit*ecalDetIdToBeRecovered*ecalRecHit)
