import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.TFileService=cms.Service("TFileService",fileName=cms.string('ResponseHistos_ak7.root'))

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(1)
        )
process.source = cms.Source("EmptySource")

##-------------------- User analyzer  --------------------------------
process.response2J   = cms.EDAnalyzer('ResponseHistos',
    filename        = cms.string('/uscms_data/d2/kkousour/7TeV/2011/Jets/mc/Summer11PythiaZ2Flat_InclusiveJetsTree_mc.root'),
    treename        = cms.string('ProcessedTree'),
    dirname         = cms.string('ak7'),
    PFBiasCorName   = cms.string('YBiasCorrection.txt'),
    yBnd            = cms.vdouble(0.0,0.5,1.0,1.5,2.0,2.05,3.0),
    yFineBnd        = cms.vdouble(0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,4.0,5.0),
    ptBnd           = cms.vdouble(30,50,100,200,500,1000,2000,3500),
    maxDR           = cms.double(0.25),
    nJets           = cms.int32(2),
    nEvents         = cms.int32(-1)
)
process.responseAll = process.response2J.clone(nJets = 100)
process.p = cms.Path(process.responseAll)

