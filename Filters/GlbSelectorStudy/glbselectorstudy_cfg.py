import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

#process.load("FWCore.MessageService.MessageLogger_cfi")

MessageLogger = cms.Service("MessageLogger",
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
        SimLevelParent = cms.untracked.PSet(
            limit = cms.untracked.int32(1000000)
        ),
        DEBUG = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        MotherSearch = cms.untracked.PSet(
            limit = cms.untracked.int32(1000000)
        )
    ),
    categories = cms.untracked.vstring('SimHitsAnlzrImproved', 
        'SinglePiFilter', 
        'MotherSearch', 
        'SimLevelParent'),
    destinations = cms.untracked.vstring('detailedInfo')
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/user/aeverett/CMSSW_2_2_5/SingleMuPt100//aeverett//SingleMuPt100_CMSSW_2_2_5_IDEAL_step1//SingleMuPt100_CMSSW_2_2_5_IDEAL_step2//49d2e03eccbedc0e7ba634f37fb81980//step2_RAW2DIGI_RECO_8.root'
    )
)

from Validation.RecoMuon.associators_cff import *
process.load("UserCode.GlbSelectorStudy.glbselectorstudy_cfi")


process.p = cms.Path(process.glbSelStudy)
