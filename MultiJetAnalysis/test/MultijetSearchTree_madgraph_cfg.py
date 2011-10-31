import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.TFileService=cms.Service("TFileService",fileName=cms.string('MultijetSearchTree_madgraph.root'))

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(1)
        )
process.source = cms.Source("EmptySource")

##-------------------- User analyzer  --------------------------------
process.multijet    = cms.EDAnalyzer('MultijetSearchTree',
    filenames       = cms.vstring('/uscms_data/d2/kkousour/7TeV/2011/Jets/mc/MadgraphSummer11_Multijets.root'),
    treename        = cms.string('ProcessedTree'),
    dirname         = cms.string('ak5'),
    triggers        = cms.vstring(),
    jetID           = cms.int32(2),
    hcalNoiseFilter = cms.int32(0),
    etaMAX          = cms.double(2.5),
    ptMIN           = cms.double(30),
    puTagMIN        = cms.double(0.2),
    nEvents         = cms.int32(-1),
    isMC            = cms.untracked.bool(True),
)

process.p = cms.Path(process.multijet)

