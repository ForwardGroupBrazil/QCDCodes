import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.TFileService=cms.Service("TFileService",fileName=cms.string('MultijetSearchTree_HT550.root'))

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(1)
        )
process.source = cms.Source("EmptySource")

##-------------------- User analyzer  --------------------------------
process.multijet    = cms.EDAnalyzer('MultijetSearchTree',
    filename        = cms.string('/uscms_data/d2/kkousour/7TeV/2011/Jets/data/Aug22nd/MultiJetsTree_data.root'),
    treename        = cms.string('ProcessedTree'),
    dirname         = cms.string('ak5'),
    triggers        = cms.vstring('HLT_HT550_v1','HLT_HT550_v2','HLT_HT550_v3','HLT_HT550_v4',
                                  'HLT_HT550_v5','HLT_HT550_v6','HLT_HT550_v7','HLT_HT550_v8'),
    jetID           = cms.int32(2),
    hcalNoiseFilter = cms.int32(1),
    etaMAX          = cms.double(2.5),
    ptMIN           = cms.double(30),
    nEvents         = cms.int32(-1)
)

process.p = cms.Path(process.multijet)

