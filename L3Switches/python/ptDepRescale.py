import FWCore.ParameterSet.Config as cms

ptRange = cms.vdouble(0 , 13, 30, 1000)
etaRange = cms.vdouble(0 ,10)

diagTerm = cms.PSet( values = cms.vdouble( 3, 5, 10 ),
action = cms.string("scale")
)
offDiagTerm = cms.PSet( values = cms.vdouble( 1, 1, 1 ),
action = cms.string("scale")
)

