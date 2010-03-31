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
process.MessageLogger.categories   += ['TrackFitters']
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
    # eventsToProcess = cms.untracked.VEventRange( '123977:11474609', ),
    fileNames = cms.untracked.vstring(
    '/store/user/aeverett/MinimumBias//CMSSW_354_PUSkim_123978//679e85c2781449dc0ef0f390bc4dc871//muonSkim_1_1.root',
    '/store/user/aeverett/MinimumBias//CMSSW_354_PUSkim_123818//679e85c2781449dc0ef0f390bc4dc871//muonSkim_1_1.root',
    '/store/user/aeverett/MinimumBias//CMSSW_354_PUSkim_123976//679e85c2781449dc0ef0f390bc4dc871//muonSkim_1_1.root',
    '/store/user/aeverett/MinimumBias//CMSSW_354_PUSkim_124030//679e85c2781449dc0ef0f390bc4dc871//muonSkim_1_1.root',
    '/store/user/aeverett/MinimumBias//CMSSW_354_PUSkim_123591//679e85c2781449dc0ef0f390bc4dc871//muonSkim_1_1.root',
    '/store/user/aeverett/MinimumBias//CMSSW_354_PUSkim_123592//679e85c2781449dc0ef0f390bc4dc871//muonSkim_1_1.root',
    #'/store/user/aeverett/CMSSW_3_4_1/SingleMuPt0_500/aeverett/SingleMuPt0_500_CMSSW_3_4_1_step1//SingleMuPt0_500_CMSSW_3_4_1_step2//bd7cda0b32f61291da1aab637c754773//step2_23.root'
    
    )
    )
    
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('myJustMuonReco.root')
    )

process.p = cms.Path(process.muonrecoComplete)
#process.p = cms.Path(process.reconstruction)

## process.load("ISpy/Analyzers/ISpy_Producer_cff")
## process.add_(
##     cms.Service("ISpyService",
##     outputFileName = cms.untracked.string('myJustMuonRecoSpy.ig'),
##     outputMaxEvents = cms.untracked.int32(100),
##     )
## )
## process.p1 = cms.Path(process.iSpy_sequence)

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("TFS_debug.root")
    )

process.load("UserCode.GlobalMatchingAnalyser.globalmatchinganalyser_cfi")

process.analyser_step = cms.Path(process.globalMatchingAnalyser) 

process.this_is_the_end = cms.EndPath(process.out)

process.schedule = cms.Schedule(process.p,process.analyser_step,process.this_is_the_end)

#
# Additional tags and parameter changes
#

process.GlobalTag.globaltag = 'GR10_P_V4::All'

#process.globalMuons.GLBTrajBuilderParameters.PtCut = 0.5

#process.globalMuons.GLBTrajBuilderParameters.GlobalMuonTrackMatcher.DeltaRCut_2 = 0.3
#process.globalMatchingAnalyser.GlobalMuonTrackMatcher.DeltaRCut_2 = 0.3

process.ttrhbwor.ComputeCoarseLocalPositionFromDisk = True
process.ttrhbwr.ComputeCoarseLocalPositionFromDisk = True

process.globalMatchingAnalyser.GlobalMuonTrackMatcher.UseTFileService = cms.untracked.bool(True)

process.globalMuons.UseTFileService = cms.untracked.bool(True)
process.globalMuons.TrackLoaderParameters.UseTFileService = cms.untracked.bool(True)
process.globalMuons.GLBTrajBuilderParameters.UseTFileService = cms.untracked.bool(True)
process.globalMuons.GLBTrajBuilderParameters.GlobalMuonTrackMatcher.UseTFileService = cms.untracked.bool(True)

process.standAloneMuons.TrackLoaderParameters.UseTFileService = cms.untracked.bool(True)

#process.globalMuons.MuonCollectionLabel = cms.InputTag("standAloneMuons","")

def customise(process):
    from Workspace.MuonRecoTreeUtility.muonRecoTreeUtilityForData_cff import insertMRTU
    insertMRTU(process)
    return (process)

process = customise(process)
