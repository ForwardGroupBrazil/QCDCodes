import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.TFileService=cms.Service("TFileService",fileName=cms.string('InclusiveHistos_pythia_puWeighted.root'))

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(1)
        )
process.source = cms.Source("EmptySource")

##-------------------- User analyzer  --------------------------------
process.minPT114   = cms.EDAnalyzer('InclusiveHistos',
    filename        = cms.string('/uscms_data/d2/kkousour/7TeV/2011/Jets/mc/Summer11PythiaZ2Flat_InclusiveJetsTree_mc.root'),
    treename        = cms.string('ProcessedTree'),
    dirname         = cms.string('ak7'),
    puFile          = cms.untracked.string('PileUpWeight.root'),
    puHisto         = cms.untracked.string('Weight'),
    yBnd            = cms.vdouble(0.0,0.5,1.0,1.5,2.0,2.5,3.0),
    ptBnd           = cms.vdouble(40,53,67,81,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,905,967,1032,1101,1172,1248,1327,1410,1497,1588,1684,1784,1890,2000,2116,2238,2366,2500,2640,2787,2941,3103,3273,3500),
    triggers        = cms.vstring(),
    minPt           = cms.vdouble(114),
    isMC            = cms.bool(True),    
    jetID           = cms.int32(2),
    hcalNoiseFilter = cms.int32(0),
    nEvents         = cms.int32(-1)
)

process.minPT196 = process.minPT114.clone(minPt = cms.vdouble(196))
process.minPT300 = process.minPT114.clone(minPt = cms.vdouble(300))
process.minPT507 = process.minPT114.clone(minPt = cms.vdouble(507))

process.p = cms.Path(process.minPT114 * process.minPT300)
#process.p = cms.Path(process.minPT114)

