import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.TFileService=cms.Service("TFileService",fileName=cms.string('InclusiveHistos.root'))

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(1)
        )
process.source = cms.Source("EmptySource")

##-------------------- User analyzer  --------------------------------
process.hlt110v1  = cms.EDAnalyzer('InclusiveHistos',
    filename     = cms.string('/data/kkousour/7TeV/2011/Jets/ProcessedTree_data_hlt110v1.root'),
    treename     = cms.string('ak7/ProcessedTree'),
    yBnd         = cms.vdouble(0.0,0.5,1.0,1.5,2.0,2.5,3.0),
    ptBnd        = cms.vdouble(30, 40, 53, 67, 81, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3500),
    minPFPt      = cms.double(30),
    minCaloPt    = cms.double(30)     
)

process.p = cms.Path(process.hlt110v1)
