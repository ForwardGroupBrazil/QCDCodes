import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.TFileService=cms.Service("TFileService",fileName=cms.string('ResponseHistos.root'))

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(1)
        )
process.source = cms.Source("EmptySource")

##-------------------- User analyzer  --------------------------------
process.response   = cms.EDAnalyzer('ResponseHistos',
    filename        = cms.string('/uscms_data/d2/kkousour/7TeV/2011/Jets/mc/Summer11Flat_ProcessedTree_mc.root'),
    treename        = cms.string('ProcessedTree'),
    dirname         = cms.string('ak5'),
    etaBnd          = cms.vdouble(-5.0,-3.0,-1.3,1.3,3.0,5.0),
    ptBnd           = cms.vdouble(30,50,100,200,500,1000,2000,3500),
    maxDR           = cms.double(0.25),
    nEvents         = cms.int32(-1)
)

process.p = cms.Path(process.response)

