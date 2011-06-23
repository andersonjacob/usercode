import FWCore.ParameterSet.Config as cms

process = cms.Process('HLT')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mix_E7TeV_FlatDist10_2011EarlyData_50ns_PoissonOOT')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('HLTrigger.Configuration.HLT_GRun_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:WHZH_GENSIM.root')
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.303.2.7 $'),
    annotation = cms.untracked.string('GEN-Fragment nevts:10'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('WHZH_RAWSIM.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RAW')
    )
)

# Additional output definition

process.mix.input.fileNames = cms.untracked.vstring(
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0010/B2F2175B-B159-E011-96FA-0024E87681FD.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0010/7CB4CB0F-EA59-E011-B9DE-00A0D1EEAA00.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0009/FC4CB683-B959-E011-9E1C-00266CF2580C.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0009/F4BDF130-B859-E011-9A12-0015178C6630.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0009/EE932991-B959-E011-BBC1-00151796C1CC.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0009/E2A32C80-B959-E011-AB98-00266CF253C4.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0009/D639D3CC-A559-E011-A5BB-00266CF85EA8.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0009/D22D24DA-B659-E011-8E1E-0015178C4BE4.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0009/D09019E2-7B59-E011-8C63-00266CF834E0.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0009/BE360FA7-8759-E011-9CC3-00266CF833EC.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0009/A0533002-B759-E011-B1CA-00151796C1CC.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0009/9EF955E3-B659-E011-AC50-00151796D708.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0009/82C56689-B959-E011-8F4D-00A0D1EEF6B8.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0009/528E75E9-B159-E011-8212-00151796C0F0.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0009/521C3627-B259-E011-BC6F-00151796C140.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0009/486182FA-B659-E011-BB64-00151796C470.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0009/220C1ADB-9B59-E011-992C-0015178C4B28.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0008/F67D4F59-0858-E011-A223-00215A491956.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0008/EE85ECB8-4258-E011-9189-00266CF8377C.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0008/D653CC4C-4B58-E011-9B56-0024E850DF73.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0008/D05C8EE0-2E58-E011-89E4-0026B93F4ADF.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0008/D03AC72A-2858-E011-ADDA-0026B95A45EC.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0008/CCE20350-4A58-E011-A42C-0024E85A4EB3.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0008/CA1242FC-4A58-E011-9237-0022198273C0.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0008/C81C116C-3658-E011-88DC-001D09645B2D.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0008/B6A75789-2C58-E011-9B2A-00266CF83958.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0008/9EDA8825-1258-E011-9A20-000AE4889FB6.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0008/7E7892B7-4858-E011-BF8C-BC305B390AA7.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0008/7C0B1D56-2F58-E011-9B8C-00188B89EEFA.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0008/706173B0-4858-E011-8BC0-BC305B390A25.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0008/6AD0394D-4A58-E011-8E43-0024E85A4EB3.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0008/62D7FA7B-4159-E011-91E2-0024E8766422.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0008/4E077A9F-AA59-E011-809A-0024E876635F.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0008/36E1C33B-4A58-E011-B606-0024E85A4FEF.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0008/36371B4C-4B58-E011-A853-0024E850DF73.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0008/2ED9C9F9-4A58-E011-9EC6-0022198273F0.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0008/0C723277-8458-E011-B6BD-000AE488BBBA.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0008/04AE43EF-1C58-E011-82C3-0015C5E5B22E.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0007/FAC66105-6059-E011-AC34-00A0D1EEF69C.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0007/E82FAC88-E559-E011-A0E9-00151796C06C.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0007/E68DB9E5-AB59-E011-8F3C-0024E8769A01.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0007/E4C0DDB3-AB59-E011-ABD2-0024E8768C71.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0007/E43DA512-EA59-E011-B871-00266CF2580C.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0007/D6016472-E859-E011-A78F-0024E86E8DB4.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0007/D4FF2FB5-AB59-E011-A2EB-0024E86E8D25.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0007/D0B20346-AC59-E011-8A56-00151796C470.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0007/CE0E9F6B-4159-E011-A1F6-0024E8768BE2.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0007/CCA4766B-E959-E011-A3B4-00266CF25F34.root',
    '/store/mc/Summer11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START311_V2-v2/0007/CA913AF2-D456-E011-9CCF-BC305B390A0B.root'
    )
# Other statements
process.GlobalTag.globaltag = 'MC_42_V12::All'

# Path and EndPath definitions
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.digitisation_step,process.L1simulation_step,process.digi2raw_step)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.endjob_step,process.RAWSIMoutput_step])

# customisation of the process.

# Automatic addition of the customisation function from Configuration.StandardSequences.DigiToRecoNoPU
from Configuration.StandardSequences.DigiToRecoNoPU import customise 

#call to customisation function customise imported from Configuration.StandardSequences.DigiToRecoNoPU
process = customise(process)

# End of customisation functions
