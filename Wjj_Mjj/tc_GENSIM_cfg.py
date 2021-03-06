# Auto generated configuration file
# using: 
# Revision: 1.303.2.7 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: Configuration/Generator/python/QCD_Pt_30_TuneZ2_7TeV_pythia6_cff.py -s GEN,SIM --conditions auto:mc --eventcontent RAWSIM --datatier GEN-SIM -n 10 --no_exec --mc --fileout tc_GENSIM.root --python_filename tc_GENSIM_cfg.py
import FWCore.ParameterSet.Config as cms

process = cms.Process('SIM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic7TeV2011Collision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(5)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.303.2.7 $'),
    annotation = cms.untracked.string('Configuration/Generator/python/QCD_Pt_30_TuneZ2_7TeV_pythia6_cff.py nevts:10'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('tc_GENSIM.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'MC_42_V12::All'

from Configuration.Generator.PythiaUEZ2Settings_cfi import *
process.generator = cms.EDFilter("Pythia6GeneratorFilter",
    pythiaPylistVerbosity = cms.untracked.int32(0),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(7000.0),
    maxEventsToPrint = cms.untracked.int32(1),
    PythiaParameters = cms.PSet(
	pythiaUESettingsBlock,
        processParameters = cms.vstring(
						 'MSEL=0',
						 'MSUB(362)=1',
						 'MSUB(368)=1',
						 'MSUB(371)=1',
						 'MSUB(376)=1',
						 'MSUB(366)=0',
						 'MSUB(372)=0',
						 'MSUB(375)=0',
						 'PMAS(329,1)=160.0d0',
						 'PMAS(330,1)=160.0d0',
						 'PMAS(333,1)=290.0d0',
						 'PMAS(334,1)=290.0d0',
						 'PMAS(335,1)=290.0d0',
						 'PMAS(368,1)=320.0d0',
						 'RTCM(12)=290.0d0',
						 'RTCM(13)=290.0d0',
						 'RTCM(48)=290.0d0',
						 'RTCM(49)=290.0d0',
						 'RTCM(50)=1000.0d0',
						 'RTCM(51)=1000.0d0',
						 'MDME(174,1)=0',
						 'MDME(175,1)=0',
						 'MDME(176,1)=0',
						 'MDME(177,1)=0',
						 'MDME(178,1)=0',
						 'MDME(179,1)=0',
						 'MDME(182,1)=0',
						 'MDME(184,1)=0',
						 'MDME(186,1)=0',
						 'MDME(190,1)=0',
						 'MDME(191,1)=0',
						 'MDME(192,1)=0',
						 'MDME(194,1)=0',
						 'MDME(195,1)=0',
						 'MDME(196,1)=0',
						 'MDME(198,1)=0',
						 'MDME(199,1)=0',
						 'MDME(200,1)=0',
						 'MDME(208,1)=0',
						 '24:ALLOFF',
						 '24:ONIFANY 11 13'),
        parameterSets = cms.vstring('pythiaUESettings', 
            'processParameters')
    )
)


process.ProductionFilterSequence = cms.Sequence(process.generator)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.RAWSIMoutput_step)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq 
