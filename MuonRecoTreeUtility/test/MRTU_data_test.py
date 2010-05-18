# Auto generated configuration file
# using: 
# Revision: 1.151 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: promptCollisionReco -s RAW2DIGI,L1Reco,RECO,DQM,ALCA:SiStripCalZeroBias --datatier RECO --eventcontent RECO --conditions GR09_H_V7OFF::All --scenario pp --no_exec --data --magField AutoFromDBCurrent -n 100
import FWCore.ParameterSet.Config as cms

process = cms.Process('MRTU')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# Destinations
process.MessageLogger.destinations += ['FitterMessages',
                                       'DetailedMessages']

# Categories
process.MessageLogger.categories   += ['GlobalMuonTrackMatcher']

process.MessageLogger.categories   += ['GlobalMuonTrajectoryBuilder',
                                       'GlobalTrajectoryBuilderBase',
                                       'GlobalMuonProducer',
                                       'MuonTrajectoryCleaner']

process.MessageLogger.categories   += ['MuonTrackLoader','TrackFitters','MuonRecoTreeUtility']

# Modules to Debug.  Note that including globalMuons and
# globalMatchingAnalyser will cause the matching printout to be repeated
# since there is a call to the TrackMatcher within the matching
# analyzer.
process.MessageLogger.debugModules += ['globalMuons','recoMuonTreeMaker']

# This one gives us some top level information about the STA, TK, GLB,
# regional TK, matched TK, and the global refit.
process.MessageLogger.DetailedMessages = cms.untracked.PSet(
    threshold  = cms.untracked.string('DEBUG'),
    default    = cms.untracked.PSet(limit = cms.untracked.int32(0)),
    INFO = cms.untracked.PSet(limit = cms.untracked.int32(0)),
    GlobalMuonTrackMatcher = cms.untracked.PSet(limit = cms.untracked.int32(-1)),
    GlobalMuonTrajectoryBuilder = cms.untracked.PSet(limit = cms.untracked.int32(-1)),
    GlobalTrajectoryBuilderBase = cms.untracked.PSet(limit = cms.untracked.int32(-1)),
    GlobalMuonProducer = cms.untracked.PSet(limit = cms.untracked.int32(-1)),
    MuonTrajectoryCleaner = cms.untracked.PSet(limit = cms.untracked.int32(-1)),

    )

# This one gives us the status of the track loader as well as the
# possibility of a hit-by-hit status for the KFFitter and KFSmoother.
# The KFFitter and KFSmoother stuff is only really useful if you suspect
# a problem with the fit or the smooth.  Such a problem would manifest
# itself, for example, in the DetailedMessages as a failure for the
# sta-tk-pair to refit as a trajectory.
process.MessageLogger.FitterMessages  = cms.untracked.PSet(
     threshold  = cms.untracked.string('DEBUG'),
     default    = cms.untracked.PSet(limit = cms.untracked.int32(0)),
     INFO = cms.untracked.PSet(limit = cms.untracked.int32(0)),
     MuonTrackLoader = cms.untracked.PSet(limit = cms.untracked.int32(0)),
     TrackFitters = cms.untracked.PSet(limit = cms.untracked.int32(0)),
     MuonRecoTreeUtility = cms.untracked.PSet(limit = cms.untracked.int32(-1)),
     )



# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration/StandardSequences/RawToDigi_Data_cff')
process.load('Configuration/StandardSequences/L1Reco_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('DQMOffline/Configuration/DQMOffline_cff')
#process.load('Configuration/StandardSequences/AlCaRecoStreams_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    annotation = cms.untracked.string('promptCollisionReco nevts:100'),
    name = cms.untracked.string('PyReleaseValidation')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
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
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
    'file:/home/ba01/u112/aeverett/work/CMSSW_3_5_6/src/FastAnalysis/Skims/test/muonSkim.data.root'
    ##     '/store/data/Commissioning10/MinimumBias/RAW-RECO/Apr1Skim_GOODCOLL-v1/0140/F4EA88D2-C83E-DF11-AF55-00261894392C.root',
    ##     '/store/data/Commissioning10/MinimumBias/RAW-RECO/Apr1Skim_GOODCOLL-v1/0140/E27B88D1-8040-DF11-B3FC-00261894391B.root',
    ##     '/store/data/Commissioning10/MinimumBias/RAW-RECO/Apr1Skim_GOODCOLL-v1/0140/E246ADD2-C83E-DF11-A9BE-00261894390C.root',
    ##     '/store/data/Commissioning10/MinimumBias/RAW-RECO/Apr1Skim_GOODCOLL-v1/0140/E0F8D4D2-C83E-DF11-9249-002618943877.root',
    ##     '/store/data/Commissioning10/MinimumBias/RAW-RECO/Apr1Skim_GOODCOLL-v1/0140/E0B003DD-C83E-DF11-8A5D-002618943947.root',
    ),
    )

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = process.FEVTEventContent.outputCommands,
    fileName = cms.untracked.string('promptCollisionReco_RAW2DIGI_L1Reco_RECO_DQM_ALCA.root'),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('RAW-RECO'),
        filterName = cms.untracked.string('')
    )
)

# Other statements
process.GlobalTag.globaltag = 'GR_R_35X_V6::All'
process.ttrhbwr.ComputeCoarseLocalPositionFromDisk = True
process.ttrhbwor.ComputeCoarseLocalPositionFromDisk = True

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.muon_reco_step = cms.Path(process.muonrecoComplete)
process.dqmoffline_step = cms.Path(process.DQMOffline)
process.endjob_step = cms.Path(process.endOfProcess)
process.out_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.muon_reco_step,process.endjob_step,process.out_step)


def customise(process):
    from Workspace.MuonRecoTreeUtility.muonRecoTreeUtilityForData_cff import insertMRTU
    insertMRTU(process)
    return (process)

process = customise(process)

process.recoMuonTreeMaker.outputFileName = "MRTU_data.root"

