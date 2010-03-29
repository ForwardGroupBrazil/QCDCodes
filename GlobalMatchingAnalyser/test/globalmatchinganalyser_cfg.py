import FWCore.ParameterSet.Config as cms

process = cms.Process("DebugMuon")
# Messages
##process.load("RecoMuon.Configuration.MessageLogger_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.MessageLogger.destinations += ['AnalyzerMessages',
                                       'DetailedMessages',
                                       'FitterMessages']
process.MessageLogger.categories   += ['MatchAnalyzer',
                                       'GlobalMuonTrackMatcher',
                                       'MuonTrackLoader']
process.MessageLogger.categories   += ['GlobalMuonTrajectoryBuilder',
                                       'GlobalTrajectoryBuilderBase',
                                       'GlobalMuonProducer',
                                       'MuonTrajectoryCleaner']
process.MessageLogger.categories   += ['MuonTrackLoader','TrackFitters']
process.MessageLogger.debugModules += ['globalMuons','globalMatchingAnalyser']


process.MessageLogger.DetailedMessages = cms.untracked.PSet(
    threshold  = cms.untracked.string('DEBUG'),
    default    = cms.untracked.PSet(limit = cms.untracked.int32(0)),
    INFO = cms.untracked.PSet(limit = cms.untracked.int32(0)),
    GlobalMuonTrackMatcher = cms.untracked.PSet(limit = cms.untracked.int32(-1)),
    GlobalMuonTrajectoryBuilder = cms.untracked.PSet(limit = cms.untracked.int32(-1)),
    GlobalTrajectoryBuilderBase = cms.untracked.PSet(limit = cms.untracked.int32(-1)),
    GlobalMuonProducer = cms.untracked.PSet(limit = cms.untracked.int32(-1)),
    MuonTrajectoryCleaner = cms.untracked.PSet(limit = cms.untracked.int32(-1)),
    #MuonTrackLoader = cms.untracked.PSet(limit = cms.untracked.int32(-1)),
    )

process.MessageLogger.AnalyzerMessages = cms.untracked.PSet(
    threshold  = cms.untracked.string('DEBUG'),
    default    = cms.untracked.PSet(limit = cms.untracked.int32(0)),
    MatchAnalyzer = cms.untracked.PSet(limit = cms.untracked.int32(-1))
    )

process.MessageLogger.FitterMessages  = cms.untracked.PSet(
    threshold  = cms.untracked.string('DEBUG'),
    default    = cms.untracked.PSet(limit = cms.untracked.int32(0)),
    INFO = cms.untracked.PSet(limit = cms.untracked.int32(0)),
    MuonTrackLoader = cms.untracked.PSet(limit = cms.untracked.int32(-1)),
    TrackFitters = cms.untracked.PSet(limit = cms.untracked.int32(-1)),
    )


# Muon Reco
process.load("RecoLocalMuon.Configuration.RecoLocalMuon_cff")

process.load("RecoMuon.Configuration.RecoMuon_cff")

process.load("Configuration.StandardSequences.Services_cff")

process.load("Configuration.StandardSequences.Geometry_cff")

process.load("Configuration.StandardSequences.MagneticField_38T_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.load("Configuration.StandardSequences.RawToDigi_cff")

process.load("Configuration.StandardSequences.Reconstruction_cff")


process.source = cms.Source(
    "PoolSource",
##     eventsToProcess = cms.untracked.VEventRange(
## ##     '124020:19131248',
## ##     '124020:28258367',
## ##     '124022:20851004',
## ##     '124022:28309645',
## ##     '124023:13014532',
## ##     '124023:21857034',
## ##     '124230:15994834',
## #
##     '124020:5100126',
##     '124020:28258367',
##     '124022:13322446',
##     '124022:28309645',
##     '124023:10537672',
##     '124024:5006249',
##     '124027:10742619',
##     '124230:11859102',
##     '124230:18631437',
    
##     ),
    fileNames = cms.untracked.vstring(
    
    # 'file:/scratch/scratch96/a/aeverett/FastAnalysis/split-mc-matched.root',

    # 'file:/scratch/scratch96/a/aeverett/FastAnalysis/data/badGlb/fastReco-badGlobalMuons_10/myFastReco.root',
    'file:/scratch/scratch96/a/aeverett/FastAnalysis/data/goodGlb/fastReco-goodGlobalMuons_21/myFastReco.root',
    
    )
    )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(2)
    )

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('myDebugMuonReco.zero.allBad.root')
                               )

process.p = cms.Path(process.muonrecoComplete)
#process.p = cms.Path(process.reconstruction)

process.load("ISpy/Analyzers/ISpy_Producer_cff")
process.add_(
    cms.Service("ISpyService",
    outputFileName = cms.untracked.string('myJustMuonRecoSpy.ig'),
    outputMaxEvents = cms.untracked.int32(100),
    )
)
process.p1 = cms.Path(process.iSpy_sequence)

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("TFS_allBad.root")
    )

process.load("UserCode.GlobalMatchingAnalyser.globalmatchinganalyser_cfi")

process.globalMatchingAnalyser.GlobalMuonTrackMatcher.UseTFileService = cms.untracked.bool(True)

process.analyser_step = cms.Path(process.globalMatchingAnalyser) 

process.this_is_the_end = cms.EndPath(process.out)

#
# Additional tags and parameter changes
#

process.GlobalTag.globaltag = 'GR09_R_34X_V2::All'

process.ttrhbwr.ComputeCoarseLocalPositionFromDisk = True

process.globalMuons.TrackLoaderParameters.UseTFileService = cms.untracked.bool(True)
process.standAloneMuons.TrackLoaderParameters.UseTFileService = cms.untracked.bool(True)
process.globalMuons.GLBTrajBuilderParameters.GlobalMuonTrackMatcher.UseTFileService = cms.untracked.bool(True)
process.globalMuons.GLBTrajBuilderParameters.UseTFileService = cms.untracked.bool(True)
process.globalMuons.UseTFileService = cms.untracked.bool(True)

#process.globalMuons.GLBTrajBuilderParameters.PtCut = 0.5

#process.globalMuons.GLBTrajBuilderParameters.GlobalMuonTrackMatcher.DeltaRCut_2 = 0.3
#process.globalMatchingAnalyser.GlobalMuonTrackMatcher.DeltaRCut_2 = 0.3

#process.globalMuons.MuonCollectionLabel = cms.InputTag("standAloneMuons","")

def customise(process):
    from Workspace.MuonRecoTreeUtility.muonRecoTreeUtility_cff import insertMRTU
    insertMRTU(process)
    return (process)

process = customise(process)
