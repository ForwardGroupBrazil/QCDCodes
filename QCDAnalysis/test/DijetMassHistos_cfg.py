import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.TFileService=cms.Service("TFileService",fileName=cms.string('DijetMassHistos.root'))

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(1)
        )
process.source = cms.Source("EmptySource")

##-------------------- User analyzer  --------------------------------
process.hlt30v1  = cms.EDAnalyzer('DijetMassHistos',
    filename     = cms.string('/uscms_data/d2/kkousour/7TeV/2011/Jets/ProcessedTree_data_hlt30v1.root'),
    treename     = cms.string('ak7/ProcessedTree'),
    yBnd         = cms.vdouble(0.0,0.5,1.0,1.5,2.0,2.5,3.0),
    massBnd        = cms.vdouble(103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7000),
    minMass      = cms.vdouble(100,100,100,100,100,100),
    maxMass      = cms.vdouble(7000,7000,7000,7000,7000,7000),
    minPt1       = cms.double(60),
    minPt2       = cms.double(30),
    usePF        = cms.bool(True)     
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
