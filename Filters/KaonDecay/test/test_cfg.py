import FWCore.ParameterSet.Config as cms

process = cms.Process("R2")
process.load("FWCore.MessageLogger.MessageLogger_cfi")

#process.MessageLogger = cms.Service("MessageLogger",
#               debugModules = cms.untracked.vstring('kaonDecayTracker','KaonDecalCal'),
#               categories = cms.untracked.vstring('KaonDecay'),
#               destinations = cms.untracked.vstring('cout'),
#          cout = cms.untracked.PSet(
#    threshold = cms.untracked.string('DEBUG'),
#    noTimeStamps = cms.untracked.bool(True),
#    noLineBreaks = cms.untracked.bool(True),
#    default = cms.untracked.PSet(
#    limit = cms.untracked.int32(0)
#    ),
#    FwkReport = cms.untracked.PSet(
#    limit = cms.untracked.int32(0)
#    ),        
#    preEventProcessing = cms.untracked.PSet(
#    limit = cms.untracked.int32(0)
#    ),
#    FwkJob = cms.untracked.PSet(
#    limit = cms.untracked.int32(0)
#    ),        
#    KaonDecay = cms.untracked.PSet(
#    limit = cms.untracked.int32(-1)
#    )
#    )
#                                    )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring('')
)

#process.PoolSource.fileNames = ['file:/scratch/scratch96/a/aeverett/biasStudy/newFilter2012/SingleKMinusPt2_200-newStep3-robot_cfg/robot_cfg-newStep3-SingleKMinusPt2_200-0062/robot_cfg-newStep3-SingleKMinusPt2_200-0062.root']

process.PoolSource.fileNames = ['/store/user/hyxu/SingleKMinusPt2_200-newStep3/newStep3-SingleKMinusPt2_200-0062.root',]

#process.load("UserCode.MissingHitsFilter.MissingHitsFilter_cfi")
#process.pFilter = cms.Path(process.missingHitsFilter)

process.kaonDecayTracker  = cms.EDFilter("KaonDecay",
                                  inMuonBound = cms.untracked.bool(False),
                                  inTrackerBound = cms.untracked.bool(True),
                                  inCal = cms.untracked.bool(False)
                                  )
process.pFilterTracker = cms.Path(process.kaonDecayTracker)

process.kaonDecayMuonBound  = cms.EDFilter("KaonDecay",
                                  inMuonBound = cms.untracked.bool(True),
                                  inTrackerBound = cms.untracked.bool(False),
                                  inCal = cms.untracked.bool(False)
                                  )
process.pFilterMuonBound = cms.Path(process.kaonDecayMuonBound)

process.kaonDecayCal  = cms.EDFilter("KaonDecay",
                                  inMuonBound = cms.untracked.bool(False),
                                  inTrackerBound = cms.untracked.bool(False),
                                  inCal = cms.untracked.bool(True)
                                  )
process.pFilterCal = cms.Path(process.kaonDecayCal)

process.kaonDecayMuon  = cms.EDFilter("KaonDecay",
                                  inMuonBound = cms.untracked.bool(False),
                                  inTrackerBound = cms.untracked.bool(False),
                                  inCal = cms.untracked.bool(False),
                                  simMuon = cms.untracked.bool(True)
                                  )
process.pFilterMuon = cms.Path(process.kaonDecayMuon)

#process.missingHitsFilter.hitCut = -1

process.ROBOT_Tk = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('pFilterTracker')
    ),
    fileName = cms.untracked.string('testTracker.root')
)

process.ROBOT_Cal = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('pFilterCal')
    ),
    fileName = cms.untracked.string('testCal.root')
)

process.ROBOT_MuonBound = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('pFilterMuonBound')
    ),
    fileName = cms.untracked.string('testMuonBound.root')
)

process.ROBOT_Other = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('!pFilterCal', '!pFilterTracker', '!pFilterMuonBound')
    ),
    fileName = cms.untracked.string('testOther.root')
)

process.ROBOT_Muon = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('pFilterMuon')
    ),
    fileName = cms.untracked.string('testMuon.root')
)

process.ROBOT_NoMuon = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('!pFilterMuon')
    ),
    fileName = cms.untracked.string('testNoMuon.root')
)

process.e = cms.EndPath(process.ROBOT_Tk+process.ROBOT_Cal+process.ROBOT_Other+process.ROBOT_Muon+process.ROBOT_NoMuon+process.ROBOT_MuonBound)
