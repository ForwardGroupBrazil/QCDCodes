import FWCore.ParameterSet.Config as cms

ptRange = cms.vdouble(0 , 13, 30, 70, 1000)
etaRange = cms.vdouble(0 , 1.0, 1.4, 10)

diagTerm = cms.PSet( values = cms.vdouble( 3, 5, 10, 10,
                                           3, 4,  7, 10,
                                           3, 5, 10, 10
                                           ))
offDiagTerm = cms.PSet( values = cms.vdouble( 1, 1, 1, 1,
                                              1, 1, 1, 1,
                                              1, 1, 1, 1
                                              ))

