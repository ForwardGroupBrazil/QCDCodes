import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

#process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger = cms.Service("MessageLogger",
               debugModules = cms.untracked.vstring('l2l3analyzer',),
               categories = cms.untracked.vstring('L2L3PtAnalyzer'),
               destinations = cms.untracked.vstring('cout'),
          cout = cms.untracked.PSet(
    threshold = cms.untracked.string('DEBUG'),
    noTimeStamps = cms.untracked.bool(True),
    noLineBreaks = cms.untracked.bool(True),
    default = cms.untracked.PSet(
    limit = cms.untracked.int32(0)
    ),
    FwkReport = cms.untracked.PSet(
    limit = cms.untracked.int32(0)
        ),        
    preEventProcessing = cms.untracked.PSet(
    limit = cms.untracked.int32(0)
    ),
    FwkJob = cms.untracked.PSet(
    limit = cms.untracked.int32(0)
    ),        
    L2L3PtAnalyzer = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
    )
    )
                                    )


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/home/ba01/u112/aeverett/scratch_rcac/muPt500_2112b.root'
    )
)

process.load("UserCode.L2L3PtAnalyzer.l2l3ptanalyzer_cfi")


process.p = cms.Path(process.l2l3analyzer)
