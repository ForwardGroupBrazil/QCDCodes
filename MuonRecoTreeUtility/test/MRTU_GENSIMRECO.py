# Auto generated configuration file
# using: 
# Revision: 1.123 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: MRTU -s L1,HLT:1E31 -n 100 --processName MHTU --conditions FrontierConditions_GlobalTag,MC_31X_V1::All --filein HLTDEBUG-file.root --no_output --customise Workspace/MuonRecoTreeUtility/muonRecoTreeUtility_customise.py --no_exec
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
process.MessageLogger.categories   += ['MuonTrackLoader','TrackFitters']

# Modules to Debug.  Note that including globalMuons and
# globalMatchingAnalyser will cause the matching printout to be repeated
# since there is a call to the TrackMatcher within the matching
# analyzer.
process.MessageLogger.debugModules += ['globalMuons']

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
     )



# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
#process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/MixingNoPileUp_cff')
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Sim_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration/StandardSequences/SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration/StandardSequences/RawToDigi_cff')
process.load('Configuration/StandardSequences/PostRecoGenerator_cff')
process.load('Configuration/StandardSequences/Validation_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    annotation = cms.untracked.string('MHTU nevts:100'),
    name = cms.untracked.string('PyReleaseValidation')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(50)
)
process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)
# Input source
process.source = cms.Source(
     "PoolSource",
     skipEvents = cms.untracked.uint32(6),
     fileNames = cms.untracked.vstring(
     '/store/mc/Summer09/ppMuX/GEN-SIM-RECO/MC_31X_V3_7TeV-v1/0012/626BBC75-C2A3-DE11-8C9F-00E0813006C6.root'
     ),
     secondaryFileNames = cms.untracked.vstring(
     )
     )

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    fileName = cms.untracked.string('MRTU_Reco.root'),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        filterName = cms.untracked.string('')
    )
)

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'START3X_V18::All'
process.ttrhbwor.ComputeCoarseLocalPositionFromDisk = True
process.ttrhbwr.ComputeCoarseLocalPositionFromDisk = True

# Path and EndPath definitions
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path(process.reconstruction)
process.muon_reco_step = cms.Path(process.muonrecoComplete)
process.validation_step = cms.Path(process.validation)
process.endjob_step = cms.Path(process.endOfProcess)
process.out_step = cms.EndPath(process.output)


# Schedule definition
process.schedule = cms.Schedule()
process.schedule.extend([process.simulation_step,process.digitisation_step,process.L1simulation_step])
process.schedule.extend([process.muon_reco_step])
#process.schedule.extend([process.validation_step])
process.schedule.extend([process.endjob_step,process.out_step])


## process.TFileService = cms.Service(
##     "TFileService",
##     fileName = cms.string("TFS_debug.root")
##     )


# Automatic addition of the customisation function
def customise(process):
     from Workspace.MuonRecoTreeUtility.muonRecoTreeUtility_cff import insertMRTU
     insertMRTU(process)
     return (process)


process = customise(process)

