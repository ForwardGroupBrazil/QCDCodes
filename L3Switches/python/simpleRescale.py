import FWCore.ParameterSet.Config as cms

ptRange = cms.vdouble(0 , 1000)
etaRange = cms.vdouble(0 ,10)
diagTerm = cms.PSet( values = cms.vdouble( 3 ))
offDiagTerm = cms.PSet( values = cms.vdouble( 1 ))
