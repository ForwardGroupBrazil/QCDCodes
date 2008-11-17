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

#from simpleRescale import *
#from ptDepRescale import *
from ptEtaRescale import *

hltL3TrajectorySeed.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.xAxis = ptRange
hltL3TrajectorySeed.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.yAxis = etaRange

hltL3TrajectorySeed.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V11 = diagTerm
hltL3TrajectorySeed.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V22 = diagTerm
hltL3TrajectorySeed.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V33 = diagTerm
hltL3TrajectorySeed.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V44 = diagTerm
hltL3TrajectorySeed.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V55 = diagTerm

hltL3TrajectorySeed.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V12 = offDiagTerm
hltL3TrajectorySeed.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V13 = offDiagTerm
hltL3TrajectorySeed.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V14 = offDiagTerm
hltL3TrajectorySeed.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V15 = offDiagTerm
hltL3TrajectorySeed.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V23 = offDiagTerm
hltL3TrajectorySeed.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V24 = offDiagTerm
hltL3TrajectorySeed.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V25 = offDiagTerm
hltL3TrajectorySeed.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V34 = offDiagTerm
hltL3TrajectorySeed.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V35 = offDiagTerm
hltL3TrajectorySeed.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V45 = offDiagTerm

from HLTrigger.Configuration.HLT_2E30_cff import MuonCkfTrajectoryBuilder

MuonCkfTrajectoryBuilder.useSeedLayer = True
#MuonCkfTrajectoryBuilder.rescaleErrorIfFail = 3.0

