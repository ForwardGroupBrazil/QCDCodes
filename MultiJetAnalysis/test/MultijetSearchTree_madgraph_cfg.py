import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.TFileService=cms.Service("TFileService",fileName=cms.string('/uscms_data/d2/kkousour/MultijetSearchTree_madgraph.root'))

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(1)
        )
process.source = cms.Source("EmptySource")

##-------------------- User analyzer  --------------------------------
process.HT500    = cms.EDAnalyzer('MultijetSearchTree',
    filenames       = cms.vstring('/uscms_data/d2/kkousour/7TeV/2011/Jets/mc/MadgraphSummer11_HT-500-1000_Multijets.root'),
    treename        = cms.string('ProcessedTree'),
    dirname         = cms.string('ak5'),
    triggers        = cms.vstring(),
    jetID           = cms.int32(2),
    hcalNoiseFilter = cms.int32(0),
    etaMAX          = cms.double(2.5),
    ptMIN           = cms.double(30),
    puTagMIN        = cms.double(0.0),
    nEvents         = cms.int32(-1),
    isMC            = cms.untracked.bool(True),
)

process.HT1000 = process.HT500.clone(filenames = cms.vstring('/uscms_data/d2/kkousour/7TeV/2011/Jets/mc/MadgraphSummer11_HT-1000_Multijets.root'))

process.p = cms.Path(process.HT500 * process.HT1000)

