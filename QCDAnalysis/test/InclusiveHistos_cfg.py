import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.TFileService=cms.Service("TFileService",fileName=cms.string('InclusiveHistos_MinPt50.root'))

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(1)
        )
process.source = cms.Source("EmptySource")

##-------------------- User analyzer  --------------------------------
process.inclusive   = cms.EDAnalyzer('InclusiveHistos',
    filename        = cms.string('/uscms_data/d2/kkousour/7TeV/2011/Jets/May17th/ProcessedTree_data.root'),
    treename        = cms.string('ProcessedTree'),
    dirname         = cms.string('ak7'),
    yBnd            = cms.vdouble(0.0,0.5,1.0,1.5,2.0,2.5,3.0),
    ptBnd           = cms.vdouble(30, 40, 55, 70, 85, 100, 115, 135, 155, 175, 195, 220, 245, 270, 300, 330, 360, 390, 430, 470, 510, 550, 600, 640, 690, 740, 790, 850, 910, 970, 1030, 1100, 1170, 1250, 1330, 1410, 1500, 1590, 1680, 1780, 1890, 2000, 2120, 2240, 2370, 2500, 2650, 2800, 2950, 3100, 3270, 3500),
    triggers        = cms.vstring('HLT_Jet30_v1','HLT_Jet30_v2',
                                'HLT_Jet60_v1','HLT_Jet60_v2',
                                'HLT_Jet80_v1','HLT_Jet80_v2',
                                'HLT_Jet110_v1','HLT_Jet110_v2',
                                'HLT_Jet150_v1','HLT_Jet150_v2',
                                'HLT_Jet190_v1','HLT_Jet190_v2',
                                'HLT_Jet240_v1','HLT_Jet240_v2',
                                'HLT_Jet300_v1',
                                'HLT_Jet370_v1','HLT_Jet370_v2'),
    #minPt           = cms.vdouble(50,50,85,85,115,115,155,155,220,220,270,270,330,330,360,470,470), 
    minPt           = cms.vdouble(50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50),
    isMC            = cms.bool(False),    
    jetID           = cms.int32(2),
    hcalNoiseFilter = cms.int32(1),
    nEvents         = cms.int32(-1)
)

process.p = cms.Path(process.inclusive)

