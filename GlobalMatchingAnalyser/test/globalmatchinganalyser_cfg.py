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
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring(
    #'/store/user/aeverett/ExoMu100824/Inclusive15/aeverett/InclusiveMu15/Inc15/e6de73f96a386852504504b09e609805/pat_579_0_UOm.root',
    #'/store/user/aeverett/CMSSW_3_4_1/TTbar_Tauola//aeverett//TTbar_Tauola_CMSSW_3_4_1_step1//TTbar_Tauola_CMSSW_3_4_1_step2//bd7cda0b32f61291da1aab637c754773//step2_13.root',
    #'/store/user/aeverett/CMSSW_3_4_1//SingleMuPt0_500//aeverett//SingleMuPt0_500_CMSSW_3_4_1_step1//SingleMuPt0_500_CMSSW_3_4_1_step2//bd7cda0b32f61291da1aab637c754773//step2_23.root',
    #'/store//mc//Summer10//QCD_Pt-15_InclusiveMu5_7TeV-pythia6//GEN-SIM-RECO//START36_V10-v1//0010//00B6A6AE-B47C-DF11-AA57-002481CFE888.root'
    #$inputFileNames
    '/store/user/aeverett/CMSSW_3_8_2//SingleMuPt0_500//aeverett//SingleMuPt0_500_CMSSW_3_8_2_step1//SingleMuPt0_500_CMSSW_3_8_2_step2//0766a8d077887d01b4cb494d974baf3d//step2_RAW2DIGI_L1Reco_RECO_VALIDATION_DQM_21_1_UhQ.root'
    )
    )
    
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(15)
    )

process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('myJustMuonReco.root'),
    SelectEvents   = cms.untracked.PSet(SelectEvents = cms.vstring('pFilter')),
    )

#process.load("RecoMuon.Configuration.RecoMuonPPonly_cff")
#process.load("RecoMuon.GlobalTrackingTools.GlobalTrackQuality_cfi") # import *
#process.muons.fillGlobalTrackQuality = True

process.p = cms.Path(process.muonrecoComplete)
#process.p = cms.Path(process.muonreco)
#process.p = cms.Path(process.globalMuons*process.glbTrackQual*process.muons)
#process.p = cms.Path(process.reconstruction)

process.tagCands = cms.EDFilter(
    "MuonRefSelector",
    src = cms.InputTag("muons"),
    cut = cms.string('isGlobalMuon > 0 & combinedQuality.localDistance > 5'),
    filter = cms.bool(True),
)
process.pFilter = cms.Path(process.tagCands)

process.load("ISpy/Analyzers/ISpy_Producer_cff")
process.add_(
    cms.Service("ISpyService",
    outputFileName = cms.untracked.string('myJustMuonRecoSpy.ig'),
    outputMaxEvents = cms.untracked.int32(100),
    )
)
process.p1 = cms.Path(process.tagCands + process.iSpy_sequence)

process.TFileService = cms.Service(
    "TFileService",
    #fileName = cms.string("OTB_$outputFileName")
    fileName = cms.string("OTB_TFS.root")
    )

process.load("UserCode.GlobalMatchingAnalyser.globalmatchinganalyser_cfi")

process.globalMatchingAnalyserGood = process.globalMatchingAnalyser.clone(useAll=2)
process.globalMatchingAnalyserBad = process.globalMatchingAnalyser.clone(useAll=3)

process.analyser_step = cms.Path(process.globalMatchingAnalyser * process.globalMatchingAnalyserGood * process.globalMatchingAnalyserBad) 

process.this_is_the_end = cms.EndPath(process.out)

process.schedule = cms.Schedule(process.p,process.pFilter,process.p1,process.analyser_step,process.this_is_the_end)

#
# Additional tags and parameter changes
#

process.GlobalTag.globaltag = 'START38_V9::All'

process.globalMuons.GLBTrajBuilderParameters.GlbRefitterParameters.MuonHitsOption = 1
#process.globalMuons.GLBTrajBuilderParameters.GlbRefitterParameters.SkipStation = 1

#process.globalMuons.GLBTrajBuilderParameters.PtCut = 0.5

#process.globalMuons.GLBTrajBuilderParameters.GlobalMuonTrackMatcher.DeltaRCut_2 = 0.3
#process.globalMatchingAnalyser.GlobalMuonTrackMatcher.DeltaRCut_2 = 0.3
#process.globalMuons.GLBTrajBuilderParameters.GlobalMuonTrackMatcher.Pt_threshold1 = 10.0
#process.globalMuons.GLBTrajBuilderParameters.GlobalMuonTrackMatcher.DeltaDCut_2 = 25.0

process.ttrhbwor.ComputeCoarseLocalPositionFromDisk = True
process.ttrhbwr.ComputeCoarseLocalPositionFromDisk = True

process.globalMatchingAnalyser.GlobalMuonTrackMatcher.UseTFileService = cms.untracked.bool(True)
process.globalMatchingAnalyserGood.GlobalMuonTrackMatcher.UseTFileService = cms.untracked.bool(True)
process.globalMatchingAnalyserBad.GlobalMuonTrackMatcher.UseTFileService = cms.untracked.bool(True)

process.globalMuons.UseTFileService = cms.untracked.bool(True)
process.globalMuons.TrackLoaderParameters.UseTFileService = cms.untracked.bool(True)
process.globalMuons.GLBTrajBuilderParameters.UseTFileService = cms.untracked.bool(True)
process.globalMuons.GLBTrajBuilderParameters.GlobalMuonTrackMatcher.UseTFileService = cms.untracked.bool(True)

process.standAloneMuons.TrackLoaderParameters.UseTFileService = cms.untracked.bool(True)

#process.globalMuons.MuonCollectionLabel = cms.InputTag("standAloneMuons","")

## def customise(process):
##     from Workspace.MuonRecoTreeUtility.muonRecoTreeUtilityForData_cff import insertMRTU
##     insertMRTU(process)
##     return (process)

## process = customise(process)
