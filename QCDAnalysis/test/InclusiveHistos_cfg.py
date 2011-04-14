import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.TFileService=cms.Service("TFileService",fileName=cms.string('InclusiveHistos.root'))

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(1)
        )
process.source = cms.Source("EmptySource")

##-------------------- User analyzer  --------------------------------
process.hlt30v1  = cms.EDAnalyzer('InclusiveHistos',
    filename     = cms.string('/uscms_data/d2/kkousour/7TeV/2011/Jets/ProcessedTree_data_hlt30v1.root'),
    treename     = cms.string('ak7/ProcessedTree'),
    yBnd         = cms.vdouble(0.0,0.5,1.0,1.5,2.0,2.5,3.0),
    ptBnd        = cms.vdouble(30, 40, 53, 67, 81, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3500),
    minPt        = cms.vdouble(30,30,30,30,30,30),
    maxPt        = cms.vdouble(3500,3500,3500,3500,3500,3500),
    usePF        = cms.bool(False)     
)
process.hlt60v1  = process.hlt30v1.clone(filename = '/uscms_data/d2/kkousour/7TeV/2011/Jets/ProcessedTree_data_hlt60v1.root')
process.hlt80v1  = process.hlt30v1.clone(filename = '/uscms_data/d2/kkousour/7TeV/2011/Jets/ProcessedTree_data_hlt80v1.root')
process.hlt110v1 = process.hlt30v1.clone(filename = '/uscms_data/d2/kkousour/7TeV/2011/Jets/ProcessedTree_data_hlt110v1.root')
process.hlt150v1 = process.hlt30v1.clone(filename = '/uscms_data/d2/kkousour/7TeV/2011/Jets/ProcessedTree_data_hlt150v1.root')
process.hlt190v1 = process.hlt30v1.clone(filename = '/uscms_data/d2/kkousour/7TeV/2011/Jets/ProcessedTree_data_hlt190v1.root')
process.hlt240v1 = process.hlt30v1.clone(filename = '/uscms_data/d2/kkousour/7TeV/2011/Jets/ProcessedTree_data_hlt240v1.root')
process.hlt370v1 = process.hlt30v1.clone(filename = '/uscms_data/d2/kkousour/7TeV/2011/Jets/ProcessedTree_data_hlt370v1.root')
process.p = cms.Path(
  process.hlt30v1 * 
  process.hlt60v1 * 
  process.hlt80v1 * 
  process.hlt110v1 * 
  process.hlt150v1 * 
  process.hlt190v1 * 
  process.hlt240v1 *
  process.hlt370v1
)
