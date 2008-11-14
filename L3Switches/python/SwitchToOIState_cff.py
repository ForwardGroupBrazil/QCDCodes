import FWCore.ParameterSet.Config as cms

from regularL3_cff import hltL3Muons
from HLTrigger.Configuration.HLT_2E30_cff import hltL3TrajectorySeed

#switch off a few things
hltL3TrajectorySeed.TrackerSeedCleaner = cms.PSet()
hltL3TrajectorySeed.TSGForRoadSearchIOpxl = cms.PSet()
hltL3TrajectorySeed.TSGFromPropagation = cms.PSet()
hltL3TrajectorySeed.MuonTrackingRegionBuilder = cms.PSet()
hltL3TrajectorySeed.TSGFromMixedPairs = cms.PSet()
hltL3TrajectorySeed.TSGFromPixelPairs = cms.PSet()
hltL3TrajectorySeed.TSGFromPixelTriplets = cms.PSet()
hltL3TrajectorySeed.TSGFromCombinedHits = cms.PSet()

#switch on the OIstate
hltL3TrajectorySeed.tkSeedGenerator = "TSGForRoadSearchOI"

ptRange = cms.vdouble(0 , 1000)
etaRange = cms.vdouble(0 ,10)

#factor3 = cms.PSet( values = cms.vdouble( 3 ))
#factor1 = cms.PSet( values = cms.vdouble( 1 ))

factor3 = cms.PSet( values = cms.vdouble( 3 ),
	action = cms.string("scale"))
factor1 = cms.PSet( values = cms.vdouble( 1 ),
	action = cms.string("scale"))

hltL3TrajectorySeed.TSGForRoadSearchOI.errorMatrixPset.action = "use"

hltL3TrajectorySeed.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.xAxis = ptRange
hltL3TrajectorySeed.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.yAxis = etaRange

hltL3TrajectorySeed.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V11 = factor3
hltL3TrajectorySeed.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V22 = factor3
hltL3TrajectorySeed.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V33 = factor3
hltL3TrajectorySeed.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V44 = factor3
hltL3TrajectorySeed.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V55 = factor3

hltL3TrajectorySeed.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V12 = factor1
hltL3TrajectorySeed.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V13 = factor1
hltL3TrajectorySeed.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V14 = factor1
hltL3TrajectorySeed.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V15 = factor1
hltL3TrajectorySeed.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V23 = factor1
hltL3TrajectorySeed.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V24 = factor1
hltL3TrajectorySeed.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V25 = factor1
hltL3TrajectorySeed.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V34 = factor1
hltL3TrajectorySeed.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V35 = factor1
hltL3TrajectorySeed.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V45 = factor1

from HLTrigger.Configuration.HLT_2E30_cff import MuonCkfTrajectoryBuilder

MuonCkfTrajectoryBuilder.useSeedLayer = True
#MuonCkfTrajectoryBuilder.rescaleErrorIfFail = 3.0

