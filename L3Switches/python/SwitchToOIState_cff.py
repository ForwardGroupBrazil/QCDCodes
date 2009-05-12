import FWCore.ParameterSet.Config as cms

from regularL3_cff import hltL3Muons
from HLTrigger.Configuration.HLT_2E30_cff import hltL3TrajectorySeedOrig

#switch off a few things
hltL3TrajectorySeedOrig.TrackerSeedCleaner = cms.PSet()
hltL3TrajectorySeedOrig.TSGForRoadSearchIOpxl = cms.PSet()
hltL3TrajectorySeedOrig.TSGFromPropagation = cms.PSet()
hltL3TrajectorySeedOrig.MuonTrackingRegionBuilder = cms.PSet()
hltL3TrajectorySeedOrig.TSGFromMixedPairs = cms.PSet()
hltL3TrajectorySeedOrig.TSGFromPixelPairs = cms.PSet()
hltL3TrajectorySeedOrig.TSGFromPixelTriplets = cms.PSet()
hltL3TrajectorySeedOrig.TSGFromCombinedHits = cms.PSet()

#switch on the OIstate
hltL3TrajectorySeedOrig.tkSeedGenerator = "TSGForRoadSearchOI"

#from simpleRescale import *
#from ptDepRescale import *
from ptEtaRescale import *

hltL3TrajectorySeedOrig.TSGForRoadSearchOI.errorMatrixPset.action = "use"

hltL3TrajectorySeedOrig.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.xAxis = ptRange
hltL3TrajectorySeedOrig.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.yAxis = etaRange

hltL3TrajectorySeedOrig.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V11 = diagTerm
hltL3TrajectorySeedOrig.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V22 = diagTerm
hltL3TrajectorySeedOrig.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V33 = diagTerm
hltL3TrajectorySeedOrig.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V44 = diagTerm
hltL3TrajectorySeedOrig.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V55 = diagTerm

hltL3TrajectorySeedOrig.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V12 = offDiagTerm
hltL3TrajectorySeedOrig.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V13 = offDiagTerm
hltL3TrajectorySeedOrig.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V14 = offDiagTerm
hltL3TrajectorySeedOrig.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V15 = offDiagTerm
hltL3TrajectorySeedOrig.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V23 = offDiagTerm
hltL3TrajectorySeedOrig.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V24 = offDiagTerm
hltL3TrajectorySeedOrig.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V25 = offDiagTerm
hltL3TrajectorySeedOrig.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V34 = offDiagTerm
hltL3TrajectorySeedOrig.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V35 = offDiagTerm
hltL3TrajectorySeedOrig.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V45 = offDiagTerm

hltL3TrajectorySeed = cms.EDFilter("L3SeedCombiner",
    seedCollections = cms.VInputTag( 
        cms.InputTag("hltL3TrajectorySeedOrig"),
    )
)

comboSeeds_seq = cms.Sequence(
    hltL3TrajectorySeedOrig
    )

from HLTrigger.Configuration.HLT_2E30_cff import MuonCkfTrajectoryBuilder

#MuonCkfTrajectoryBuilder.useSeedLayer = True
#MuonCkfTrajectoryBuilder.rescaleErrorIfFail = 3.0

