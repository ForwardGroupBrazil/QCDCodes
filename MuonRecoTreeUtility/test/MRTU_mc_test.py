# Auto generated configuration file
# using: 
# Revision: 1.151 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: promptCollisionReco -s RAW2DIGI,L1Reco,RECO,DQM,ALCA:SiStripCalZeroBias --datatier RECO --eventcontent RECO --conditions GR09_H_V7OFF::All --scenario pp --no_exec --data --magField AutoFromDBCurrent -n 100
import FWCore.ParameterSet.Config as cms

process = cms.Process('MRTU')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/RawToDigi_cff')
process.load('Configuration/StandardSequences/L1Reco_cff') #
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('DQMOffline/Configuration/DQMOffline_cff')
#process.load('Configuration/StandardSequences/AlCaRecoStreams_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.151 $'),
    annotation = cms.untracked.string('promptCollisionReco nevts:100'),
    name = cms.untracked.string('PyReleaseValidation')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('OtherCMS',
                                    'StdException',
                                    'Unknown',
                                    'BadAlloc',
                                    'BadExceptionType',
                                    'ProductNotFound',
                                    'DictionaryNotFound',
                                    'InsertFailure',
                                    'Configuration',
                                    'LogicError',
                                    'UnimplementedFeature',
                                    'InvalidReference',
                                    'NullPointerError',
                                    'NoProductSpecified',
                                    'EventTimeout',
                                    'EventCorruption',
                                    'ScheduleExecutionFailure',
                                    'EventProcessorFailure',
                                    'FileInPathError',
                                    'FileOpenError',
                                    'FileReadError',
                                    'FatalRootError',
                                    'MismatchedInputFiles',
                                    'ProductDoesNotSupportViews',
                                    'ProductDoesNotSupportPtr',
                                    'NotFound')    
    )
# Input source

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    '/store/user/aeverett/aeverett//Commisioning//MC900//PurdueSkim0307_01//muonSkim_1.root',
    )
                            )

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
   splitLevel = cms.untracked.int32(0),
   outputCommands = process.FEVTEventContent.outputCommands,
   fileName = cms.untracked.string('reco_outputFileName.root'),
   dataset = cms.untracked.PSet(
       dataTier = cms.untracked.string('RAW-RECO'),
       filterName = cms.untracked.string('')
   )
)



# Other statements
process.GlobalTag.globaltag = 'START3X_V18::All'
process.ttrhbwr.ComputeCoarseLocalPositionFromDisk = True

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.muonreconstruction_step = cms.Path(process.muonrecoComplete)
process.dqmoffline_step = cms.Path(process.DQMOffline)
process.endjob_step = cms.Path(process.endOfProcess)
process.out_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.muonreconstruction_step,process.endjob_step,process.out_step)



def customise(process):
    from Workspace.MuonRecoTreeUtility.muonRecoTreeUtility_cff import insertMRTU
    insertMRTU(process)
    return (process)

process = customise(process)

process.recoMuonTreeMaker.outputFileName = "MRTU_outputFileName.root"



