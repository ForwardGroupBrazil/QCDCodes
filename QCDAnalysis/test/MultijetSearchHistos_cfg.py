import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.TFileService=cms.Service("TFileService",fileName=cms.string('MultijetSearchHistos_superTight_HT450.root'))

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(1)
        )
process.source = cms.Source("EmptySource")

##-------------------- User analyzer  --------------------------------
process.multijet    = cms.EDAnalyzer('MultijetSearchHistos',
    filename        = cms.string('/uscms_data/d2/kkousour/7TeV/2011/Jets/data/July8th/InclusiveJetsTree_HT_data.root'),
    treename        = cms.string('ProcessedTree'),
    dirname         = cms.string('ak5'),
    triggers        = cms.vstring('HLT_HT450_v1','HLT_HT450_v2','HLT_HT450_v3','HLT_HT450_v4',
                                  'HLT_HT450_v5','HLT_HT450_v6','HLT_HT450_v7'),
    minPt           = cms.vdouble(160,120,90,70,50,40,30,30),
    maxEta          = cms.double(2.5),
    rank            = cms.int32(8),
    jetID           = cms.int32(3),
    hcalNoiseFilter = cms.int32(1),
    nEvents         = cms.int32(-1)
)

process.p = cms.Path(process.multijet)

