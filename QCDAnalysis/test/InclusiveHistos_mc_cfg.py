import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.TFileService=cms.Service("TFileService",fileName=cms.string('InclusiveHistos_MinPt100_mc.root'))

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(1)
        )
process.source = cms.Source("EmptySource")

##-------------------- User analyzer  --------------------------------
process.inclusive = cms.EDAnalyzer('InclusiveHistos',
    filename      = cms.string('/uscms_data/d2/kkousour/7TeV/Spring11/ProcessedTree_mc.root'),
    treename      = cms.string('ProcessedTree'),
    dirname       = cms.string('ak7'),
    yBnd          = cms.vdouble(0.0,0.5,1.0,1.5,2.0,2.5,3.0),
    ptBnd         = cms.vdouble(30, 40, 55, 70, 85, 100, 115, 135, 155, 175, 195, 220, 245, 270, 300, 330, 360, 390, 430, 470, 510, 550, 600, 640, 690, 740, 
790, 850, 910, 970, 1030, 1100, 1170, 1250, 1330, 1410, 1500, 1590, 1680, 1780, 1890, 2000, 2120, 2240, 2370, 2500, 2650, 2800, 2950, 3100, 3270, 3500),
    minPt         = cms.double(100),
    triggers      = cms.vstring(),
    isMC          = cms.bool(True)     
)

process.p = cms.Path(process.inclusive)

