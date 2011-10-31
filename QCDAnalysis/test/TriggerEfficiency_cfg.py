import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.TFileService=cms.Service("TFileService",fileName=cms.string('TriggerEfficiency.root'))

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(1)
        )
process.source = cms.Source("EmptySource")

##-------------------- User analyzer  --------------------------------
process.efficiency  = cms.EDAnalyzer('TriggerEfficiency',
    filenames       = cms.vstring(
                                   '/uscms_data/d2/kkousour/7TeV/2011/Jets/data/Oct11th/Jet_Run2011A_May10_ProcessedTree_data.root',
                                   '/uscms_data/d2/kkousour/7TeV/2011/Jets/data/Oct11th/Jet_Run2011A_PromptV4_ProcessedTree_data.root',
                                   '/uscms_data/d2/kkousour/7TeV/2011/Jets/data/Oct11th/Jet_Run2011A_Aug05_ProcessedTree_data.root',
                                   '/uscms_data/d2/kkousour/7TeV/2011/Jets/data/Oct11th/Jet_Run2011A_PromptV6_ProcessedTree_data.root',
                                   '/uscms_data/d2/kkousour/7TeV/2011/Jets/data/Oct11th/Jet_Run2011B_PromptV1_ProcessedTree_data.root'
    ),
    treename        = cms.string('ProcessedTree'),
    dirname         = cms.string('ak7'),
    refTrigger      = cms.vstring('HLT_Jet30_v1','HLT_Jet30_v2','HLT_Jet30_v3','HLT_Jet30_v4','HLT_Jet30_v5','HLT_Jet30_v6',
                                  'HLT_Jet60_v1','HLT_Jet60_v2','HLT_Jet60_v3','HLT_Jet60_v4','HLT_Jet60_v5','HLT_Jet60_v6',
                                  'HLT_Jet110_v1','HLT_Jet110_v2','HLT_Jet110_v3','HLT_Jet110_v4','HLT_Jet110_v5','HLT_Jet110_v6',
                                  'HLT_Jet190_v1','HLT_Jet190_v2','HLT_Jet190_v3','HLT_Jet190_v4','HLT_Jet190_v5','HLT_Jet190_v6',
                                  'HLT_Jet240_v1','HLT_Jet240_v2','HLT_Jet240_v3','HLT_Jet240_v4','HLT_Jet240_v5','HLT_Jet240_v6'
                                 ),    
    yBnd            = cms.vdouble(0.0,0.5,1.0,1.5,2.0,2.5,3.0),
    ptBnd           = cms.vdouble(20,40,53,67,81,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,905,967,1032,1101,1172,1248,1327,1410,1497,1588,1684,1784,1890,2000,2116,2238,2366,2500,2640,2787,2941,3103,3273,3500),
    minPt           = cms.double(20),
    L1Pt            = cms.vdouble(36,36,36,36,36,36,68,68,68,68,68,68,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92),
    HLTPt           = cms.vdouble(60,60,60,60,60,60,110,110,110,110,110,110,190,190,190,190,190,190,240,240,240,240,240,240,370,370,370,370,370,370),
    jetID           = cms.int32(2),
    hcalNoiseFilter = cms.int32(1),
    nEvents         = cms.int32(-1)
)

process.p = cms.Path(process.efficiency)

