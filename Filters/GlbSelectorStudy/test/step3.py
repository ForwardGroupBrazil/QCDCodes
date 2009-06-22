# Auto generated configuration file
# using: 
# Revision: 1.99.2.3 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: step2 -s RAW2DIGI,RECO,POSTRECO,ALCA:MuAlCalIsolatedMu+RpcCalHLT,VALIDATION --relval 25000,100 --datatier GEN-SIM-RECO --eventcontent RECOSIM --conditions FrontierConditions_GlobalTag,IDEAL_V9::All --filein file:XXX --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('VAL3a')

## process.MessageLogger = cms.Service("MessageLogger",
##     detailedInfo = cms.untracked.PSet(
##         threshold = cms.untracked.string('INFO'),
##         INFO = cms.untracked.PSet(
##             limit = cms.untracked.int32(0)
##         ),
##         extension = cms.untracked.string('.txt')
##     ),
##     debug = cms.untracked.PSet(
##         threshold = cms.untracked.string('DEBUG'),
##         INFO = cms.untracked.PSet(
##             limit = cms.untracked.int32(0)
##         ),
##         SimLevelParent = cms.untracked.PSet(
##             limit = cms.untracked.int32(1000000)
##         ),
##         DEBUG = cms.untracked.PSet(
##             limit = cms.untracked.int32(0)
##         ),
##         MotherSearch = cms.untracked.PSet(
##             limit = cms.untracked.int32(1000000)
##         ),
##         GlbSelectorStudy = cms.untracked.PSet(
##             limit = cms.untracked.int32(1000000)
##         )
##     ),
##     debugModules  = cms.untracked.vstring('glbSelStudy',),
##     categories = cms.untracked.vstring('SimHitsAnlzrImproved', 
##         'GlbSelectorStudy', 
##         'MotherSearch', 
##         'SimLevelParent'),
##     destinations = cms.untracked.vstring('detailedInfo3')
## )

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.categories = ['MotherSearch', 'GlbSelectorStudy','MuonIdentification']
process.MessageLogger.debugModules = ['glbSelStudy','muons']
process.MessageLogger.cout = cms.untracked.PSet(
    threshold = cms.untracked.string('DEBUG'),
    default = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
    ),
    MotherSearch = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
    ),
    GlbSelectorStudy = cms.untracked.PSet(
        limit = cms.untracked.int32(-1)
    ),
        MuonIdentification = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
        ),
)
process.MessageLogger.cerr = cms.untracked.PSet(
    placeholder = cms.untracked.bool(True)
)



# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
#process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/MixingNoPileUp_cff')
process.load('Configuration/StandardSequences/GeometryPilot2_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/RawToDigi_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/StandardSequences/PostRecoGenerator_cff')
#process.load('Configuration/StandardSequences/AlCaRecoStreams_cff')
process.load('Configuration/StandardSequences/Validation_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.7 $'),
    annotation = cms.untracked.string('step2 nevts:1'),
    name = cms.untracked.string('PyReleaseValidation')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)
process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)
# Input source
process.source = cms.Source("PoolSource",
                            skipEvents = cms.untracked.uint32(0),         
#                                fileNames = cms.untracked.vstring('file:step2K_RAW2DIGI_RECO_POSTRECO_ALCA_VALIDATION.root')
#                            fileNames = cms.untracked.vstring('file:/home/ba01/u112/aeverett/scratch_rcac/fullOutput.QCDpt800.root')
                            #   fileNames = cms.untracked.vstring('file:/home/ba01/u112/aeverett/scratch_rcac/step3.root')
#                                    fileNames = cms.untracked.vstring('/store/user/aeverett/SingleKPt2_200_CMSSW_2_2_5_IDEAL_step1//SingleKPt2_200_CMSSW_2_2_5_IDEAL_step1//abc92f2ca79c035bec8931c1df74b704//SingleKPt2_200_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_162.root',)
                            fileNames = cms.untracked.vstring('/store/user/aeverett//CMSSW_2_2_5//SingleKPt2_200//aeverett//SingleKPt2_200_CMSSW_2_2_5_IDEAL_step1//SingleKPt2_200_CMSSW_2_2_5_IDEAL_step2//49d2e03eccbedc0e7ba634f37fb81980//step2_RAW2DIGI_RECO_295.root',),
                            
                            secondaryFileNames = cms.untracked.vstring('/store/user/aeverett/SingleKPt2_200_CMSSW_2_2_5_IDEAL_step1/SingleKPt2_200_CMSSW_2_2_5_IDEAL_step1/abc92f2ca79c035bec8931c1df74b704/SingleKPt2_200_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_329.root',),

                            
)

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    outputCommands = process.RECOSIMEventContent.outputCommands,
    fileName = cms.untracked.string('step3.root'),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    )
)


# Additional output definition
process.output.outputCommands.append("keep *")
process.output.outputCommands.append("keep *_MEtoEDMConverter_*_*")
# Other statements
process.TFileService = cms.Service("TFileService", fileName = cms.string("histo3.root"))

process.GlobalTag.globaltag = 'IDEAL_V11::All'
process.validation = cms.Sequence(process.mix+process.globaldigisanalyze*process.globalhitsanalyze*process.globalrechitsanalyze*process.globalValidation)

#process.load("TrackingTools.TrackRefitter.TracksToTrajectories_cff")

#####
process.muonCand = cms.EDFilter(
    "MuonRefSelector",
    src = cms.InputTag("muons"),
    cut = cms.string('isGood("GMGoldAll") > 0 && isGood("GlobalMuonPromptTight")')
    )
process.tpToCandTrackAssociation = cms.EDProducer(
    'TrackAssociatorEDProducer',
    associator = cms.string('TrackAssociatorByDeltaR'),
    label_tp = cms.InputTag('mergedtruth', 'MergedTrackTruth'),
    label_tr = cms.InputTag('globalMuons')
    #    label_tr = cms.InputTag('muonGlb')
    )
import Validation.RecoMuon.MultiTrackValidator_cfi
process.glbCandTrackVTrackAssoc = Validation.RecoMuon.MultiTrackValidator_cfi.multiTrackValidator.clone()
 
process.glbCandTrackVTrackAssoc.associatormap = 'tpToCandTrackAssociation'
process.glbCandTrackVTrackAssoc.associators = ('TrackAssociatorByDeltaR',)
process.glbCandTrackVTrackAssoc.label = ('globalMuons',)

import Validation.RecoMuon.RecoMuonValidator_cfi
process.candMuonVTrackAssoc = Validation.RecoMuon.RecoMuonValidator_cfi.recoMuonValidator.clone()
 
process.candMuonVTrackAssoc.subDir = 'RecoMuonV/CandMuon_TrackAssoc'

process.candMuonVTrackAssoc.trkMuLabel = 'generalTracks'
process.candMuonVTrackAssoc.staMuLabel = 'standAloneMuons:UpdatedAtVtx'
process.candMuonVTrackAssoc.glbMuLabel = 'globalMuons'
process.candMuonVTrackAssoc.muonLabel = 'muonCand'
 
process.candMuonVTrackAssoc.trkMuAssocLabel = 'tpToTkmuTrackAssociation'
process.candMuonVTrackAssoc.staMuAssocLabel = 'tpToStaTrackAssociation'
process.candMuonVTrackAssoc.glbMuAssocLabel = 'tpToCandTrackAssociation'

import UserCode.GlbSelectorStudy.glbselectorstudy_cff
process.goldenSelStudy = UserCode.GlbSelectorStudy.glbselectorstudy_cfi.glbSelStudy.clone()
process.goldenSelStudy.doAssoc = False
process.goldenSelStudy.trkMuAssocLabel = "tpToTkmuTrackAssociation"
process.goldenSelStudy.staMuAssocLabel = "tpToStaTrackAssociation"
process.goldenSelStudy.glbMuAssocLabel = "tpToCandTrackAssociation"
process.goldenSelStudy.glbMuLabel = "globalMuons"
process.goldenSelStudy.muonLabel = "muonCand"
#process.glbSelStudy.tpSelector.tip = 10000
#process.glbSelStudy.tpSelector.lip = 10000

#####
process.load("UserCode.GlbSelectorStudy.glbselectorstudy_cff")
process.glbSelStudy.doAssoc = False
process.glbSelStudy.trkMuAssocLabel = "tpToTkmuTrackAssociation"
process.glbSelStudy.staMuAssocLabel = "tpToStaTrackAssociation"
process.glbSelStudy.glbMuAssocLabel = "tpToGlbTrackAssociation"
#process.glbSelStudy.tpSelector.tip = 10000
#process.glbSelStudy.tpSelector.lip = 10000
process.p = cms.Path(process.muonCand+process.muonAssociation_seq*process.tpToCandTrackAssociation+(process.glbCandTrackVTrackAssoc*process.candMuonVTrackAssoc)+process.glbSelStudy*process.goldenSelStudy)


# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path(process.reconstruction)
process.postreco_step = cms.Path(process.postreco_generator)

process.validation_step = cms.Path(process.validation)
process.endjob_step = cms.Path(process.endOfProcess)
process.out_step = cms.EndPath(process.output)

process.load("SimGeneral.MixingModule.mixNoPU_cfi")
process.load("SimGeneral.TrackingAnalysis.trackingParticles_cfi")
#process.mergedtruth.vertexDistanceCut = 2000
process.mergedtruth.volumeRadius = 1200.0 # 1200.0 4K
process.mergedtruth.volumeZ = 3000.0 # 3000.0 6K

process.reTP_step = cms.Path(process.mix*process.trackingParticles)


process.load("SimTracker.TrackHistory.TrackClassifier_cff")

process.trackHistoryAnalyzer = cms.EDAnalyzer("TrackHistoryAnalyzer",
    process.trackClassifier
)
#process.trackHistoryAnalyzer.trackProducer = 'globalMuon'
process.trackHistoryAnalyzer.trackAssociator = 'TrackAssociatorByHits'
process.history = cms.Path(process.trackHistoryAnalyzer)

process.add_( 
  cms.Service("TFileService",
      fileName = cms.string("test.root")
  )
)
process.trackCategoriesAnalyzer = cms.EDFilter("TrackCategoriesAnalyzer",
    process.trackClassifier,
    minimumNumberOfHits = cms.untracked.int32(8),
    minimumTransverseMomentum = cms.untracked.double(1.),
    minimumNumberOfPixelHits = cms.untracked.int32(2),
    maximumChiSquared = cms.untracked.double(5.),
    trackQualityClass = cms.untracked.string('loose')
)

process.p2 = cms.Path(process.trackCategoriesAnalyzer)



# Schedule definition
#process.schedule = cms.Schedule(process.reTP_step,process.raw2digi_step,process.reconstruction_step,process.p) #,process.postreco_step,process.validation_step,process.endjob_step,process.out_step)

#process.schedule = cms.Schedule(process.validation_step,process.p)

process.schedule = cms.Schedule(process.reTP_step,process.raw2digi_step,process.reconstruction_step,process.postreco_step,process.validation_step,process.p,process.endjob_step,process.out_step)
#process.schedule = cms.Schedule(process.p)

process.glbSelStudy.trackProducer = "globalMuons"
process.glbSelStudy.trackAssociator = "TrackAssociatorByPosition"


