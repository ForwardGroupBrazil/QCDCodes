# Import configurations
import FWCore.ParameterSet.Config as cms

# set up process
process = cms.Process("DiMuonIntoNtuples")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('PATLayer0Summary')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    default          = cms.untracked.PSet( limit = cms.untracked.int32(0)  ),
    PATLayer0Summary = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
)


## process.MessageLogger = cms.Service("MessageLogger",
##     detailedInfo = cms.untracked.PSet(
##         threshold = cms.untracked.string('INFO'),
##         INFO = cms.untracked.PSet(
##             limit = cms.untracked.int32(-1)
##         ),
##         extension = cms.untracked.string('.txt')
##     ),
##     debug = cms.untracked.PSet(
##         threshold = cms.untracked.string('DEBUG'),
##         INFO = cms.untracked.PSet(
##             limit = cms.untracked.int32(0)
##         ),
##         NTUPLE = cms.untracked.PSet(
##             limit = cms.untracked.int32(-1)
##         ),
##         DEBUG = cms.untracked.PSet(
##             limit = cms.untracked.int32(0)
##         ),
##         extension = cms.untracked.string('.txt')
##     ),
##     debugModules  = cms.untracked.vstring('muontontuple',),
##     categories = cms.untracked.vstring(
##         'NTUPLE',
##         ),
##     destinations = cms.untracked.vstring('detailedInfo','debug')
## )


process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# Load geometry
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('IDEAL_V11::All')
process.load("Configuration.StandardSequences.MagneticField_38T_cff")

# this defines the input files
from testInput_cfi import *

# this inputs the input files from the previous function
process.source = testInput()

# set the number of events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)


# input pat sequences
process.load("PhysicsTools.PatAlgos.patLayer0_cff")
process.load("PhysicsTools.PatAlgos.patLayer1_cff")

# input pat analyzer sequence
process.load("UserCode.MuonToNtuple.muontontuple_cfi")

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
process.HLTAnalyzerEndpath = cms.EndPath( process.hltL1gtTrigReport + process.hltTrigReport )


process.load("PhysicsTools.PatAlgos.patLayer1_EventContent_cff")
# setup event content: keep everything before PAT
#aaa process.patEventContent = cms.PSet(
#aaa    outputCommands = cms.untracked.vstring('keep *')
#aaa )
# extend event content to include PAT objects
#aaa process.patEventContent.outputCommands.extend(process.patLayer1EventContent.outputCommands)
#process.patLayer1EventContent.outputCommands.extend(['keep *_offlinePrimaryVerticesFromCTFTracks_*_*'])

#process.load("ElectroWeakAnalysis.ZReco.mcTruthForDimuons_cff")
#process.load("ElectroWeakAnalysis.ZReco.dimuonsHLTFilter_cfi")
#process.load("ElectroWeakAnalysis.ZReco.dimuonsFilter_cfi")
#process.load("ElectroWeakAnalysis.ZReco.dimuonsOneTrackFilter_cfi")

process.p = cms.Path(
    process.patLayer0
    *process.patLayer1
    *process.muontontuple
    #+process.mcTruthForDimuons
)

from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('inclusiveMuSkim.root'),
    # save only events passing the full path
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    # save PAT Layer 1 output
    outputCommands = cms.untracked.vstring(
#        'drop *',
        'keep *_genMuons_*_*',
        'keep *_selectedLayer1Muons_*_*',
        'keep *_slimLayer1Jets_*_*',
        'keep *_goodTracks_*_*',
        'keep *_layer1METs_*_*',
        'keep *_layer1PFMETs_*_*',
        'keep l1extraL1MuonParticles_hltL1extraParticles_*_*',
        'keep *_patHLTMu3_*_*',
        'keep *_patHLTMu9_*_*',
        'keep *_patHLTDiMu3_*_*'
    ),
    dropMetaDataForDroppedData = cms.untracked.bool(True),
)
#process.outpath = cms.EndPath(process.out)