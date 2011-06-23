import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.TFileService=cms.Service("TFileService",fileName=cms.string('InclusiveHistos_data.root'))

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(1)
        )
process.source = cms.Source("EmptySource")

##-------------------- User analyzer  --------------------------------
process.inclusive   = cms.EDAnalyzer('InclusiveHistos',
    filename        = cms.string('/uscms_data/d2/kkousour/7TeV/2011/Jets/data/June21st/InclusiveJetsTree_data.root'),
    treename        = cms.string('ProcessedTree'),
    dirname         = cms.string('ak7'),
    yBnd            = cms.vdouble(0.0,0.5,1.0,1.5,2.0,2.5,3.0),
    ptBnd           = cms.vdouble(40,53,67,81,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,905,967,1032,1101,1172,1248,1327,1410,1497,1588,1684,1784,1890,2000,2116,2238,2366,2500,2640,2787,2941,3103,3273,3500),
    triggers        = cms.vstring('HLT_Jet30_v1','HLT_Jet30_v2','HLT_Jet30_v3','HLT_Jet30_v4',
                                  'HLT_Jet60_v1','HLT_Jet60_v2','HLT_Jet60_v3','HLT_Jet60_v4',
                                  'HLT_Jet80_v1','HLT_Jet80_v2','HLT_Jet80_v3','HLT_Jet80_v4',
                                  'HLT_Jet110_v1','HLT_Jet110_v2','HLT_Jet110_v3','HLT_Jet110_v4',
                                  'HLT_Jet150_v1','HLT_Jet150_v2','HLT_Jet150_v3','HLT_Jet150_v4',
                                  'HLT_Jet190_v1','HLT_Jet190_v2','HLT_Jet190_v3','HLT_Jet190_v4',
                                  'HLT_Jet240_v1','HLT_Jet240_v2','HLT_Jet240_v3','HLT_Jet240_v4',
                                  'HLT_Jet300_v1', 'HLT_Jet300_v2','HLT_Jet300_v3',
                                  'HLT_Jet370_v1','HLT_Jet370_v2','HLT_Jet370_v3','HLT_Jet370_v4'),
    minPt           = cms.vdouble(53,53,53,53,81,81,81,81,114,114,114,114,153,153,153,153,220,220,220,220,
                                  272,272,272,272,330,330,330,330,362,362,362,468,468,468,468), 
    #minPt           = cms.vdouble(50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50),
    isMC            = cms.bool(False),    
    jetID           = cms.int32(2),
    hcalNoiseFilter = cms.int32(1),
    nEvents         = cms.int32(-1)
)

process.p = cms.Path(process.inclusive)

