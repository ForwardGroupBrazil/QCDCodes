# Import configurations
import FWCore.ParameterSet.Config as cms

# set up process
process = cms.Process("DiMuonIntoNtuples")

# initialize MessageLogger and output report
## process.load("FWCore.MessageLogger.MessageLogger_cfi")
## process.MessageLogger.cerr.threshold = 'INFO'
## process.MessageLogger.cerr.INFO = cms.untracked.PSet(
##     default          = cms.untracked.PSet( limit = cms.untracked.int32(0)  ),
## )


process.MessageLogger = cms.Service("MessageLogger",
    detailedInfo = cms.untracked.PSet(
        threshold = cms.untracked.string('INFO'),
        INFO = cms.untracked.PSet(
            limit = cms.untracked.int32(-1)
        ),
        extension = cms.untracked.string('.txt')
    ),
    debug = cms.untracked.PSet(
        threshold = cms.untracked.string('DEBUG'),
        INFO = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        NTUPLE = cms.untracked.PSet(
            limit = cms.untracked.int32(-1)
        ),
        TrackTransformer = cms.untracked.PSet(
            limit = cms.untracked.int32(-1)
        ),
#        TracksToTrajecories = cms.untracked.PSet(
#            limit = cms.untracked.int32(-1)
#        ),
        GlobalMuonRefitter = cms.untracked.PSet(
            limit = cms.untracked.int32(-1)
        ),
        DEBUG = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        extension = cms.untracked.string('.txt')
    ),
    debugModules  = cms.untracked.vstring('aodmuontontuple',),
    categories = cms.untracked.vstring(
        'NTUPLE',"TrackTransformer","GlobalMuonRefitter"
        ),
    destinations = cms.untracked.vstring('detailedInfo','debug')
)


process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# Load geometry
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
## # aaa why this next one needed?
## process.load('Configuration/StandardSequences/Services_cff')
## #process.load('FWCore/MessageService/MessageLogger_cfi')
## process.load('Configuration/StandardSequences/MixingNoPileUp_cff')
## process.load('Configuration/StandardSequences/GeometryPilot2_cff')
## #process.load('Configuration/StandardSequences/MagneticField_38T_cff')
## process.load('Configuration/StandardSequences/RawToDigi_cff')
from Configuration.StandardSequences.Reconstruction_cff import*
process.load('Configuration/StandardSequences/Reconstruction_cff')

## process.raw2digi_step = cms.Path(process.RawToDigi)
## process.reconstruction_step = cms.Path(process.reconstruction)

process.GlobalTag.globaltag = cms.string('MC_31X_V1::All')
process.load("Configuration.StandardSequences.MagneticField_38T_cff")

# this defines the input files
from testInput_cfi import *

# this inputs the input files from the previous function
process.source = testInput()

# set the number of events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

# input pat sequences
#process.load("PhysicsTools.PatAlgos.patLayer0_cff")
#process.load("PhysicsTools.PatAlgos.patLayer1_cff")
process.load('PhysicsTools.PatAlgos.patSequences_cff')

# input AOD analyzer sequence
process.load("UserCode.MuonToNtuple.aodmuontontuple_cff")

# talk to TFileService for output histograms
#process.TFileService = cms.Service("TFileService",
#    fileName = cms.string('test.root')
#)

# request a summary at the end of the file
process.options = cms.untracked.PSet(
     wantSummary = cms.untracked.bool(True)
)

process.hltL1gtTrigReport = cms.EDAnalyzer( "L1GtTrigReport",
    UseL1GlobalTriggerRecord = cms.bool( False ),
    L1GtRecordInputTag = cms.InputTag( "hltGtDigis" )
)
process.hltTrigReport = cms.EDAnalyzer( "HLTrigReport",
    HLTriggerResults = cms.InputTag( 'TriggerResults','','HLT' )
)
#process.HLTAnalyzerEndpath = cms.EndPath( process.hltL1gtTrigReport + process.hltTrigReport )

#process.load("PhysicsTools.PatAlgos.patLayer1_EventContent_cff")
# setup event content: keep everything before PAT
#process.patEventContent = cms.PSet(
#    outputCommands = cms.untracked.vstring('keep *')
#)
# extend event content to include PAT objects
#process.patEventContent.outputCommands.extend(process.patLayer1EventContent.outputCommands)

## #process.load("RecoTracker.TrackProducer.RefitterWithMaterial_cff")
## process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
## process.TrackRefitter.src = "globalMuons"
## process.TrackRefitter.TTRHBuilder = "WithoutRefit"
## process.TrackRefitter.TrajectoryInEvent = False
## process.load("RecoTracker.TransientTrackingRecHit.TransientTrackingRecHitBuilderWithoutRefit_cfi")


#process.aodmuontontuple.RefitterParameters.TrackerRecHitBuilder= 'WithoutRefit'

process.p2 =cms.Path(process.siPixelRecHits+process.siStripMatchedRecHits+process.ckftracks) #+process.muonrecowith_TeVRefinemen)

process.p = cms.Path(
    process.siPixelRecHits+process.siStripMatchedRecHits+
    process.patDefaultSequence
#    process.patLayer0
#    *process.patLayer1
    *process.aodmuontontuple
    #+process.mcTruthForDimuons
)

process.schedule = cms.Schedule(process.p2,process.p)

