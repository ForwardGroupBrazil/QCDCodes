import FWCore.ParameterSet.Config as cms

from regularL3_cff import hltL3Muons
from HLTrigger.Configuration.HLT_2E30_cff import hltL3TrajectorySeed

hltL3TrajectorySeed.TSGForRoadSearchIOpxl = cms.PSet()
hltL3TrajectorySeed.TSGForRoadSearchOI = cms.PSet()
hltL3TrajectorySeed.MuonTrackingRegionBuilder = cms.PSet()
hltL3TrajectorySeed.TSGFromMixedPairs = cms.PSet()
hltL3TrajectorySeed.TSGFromPixelPairs = cms.PSet()
hltL3TrajectorySeed.TSGFromPixelTriplets = cms.PSet()
hltL3TrajectorySeed.TSGFromCombinedHits = cms.PSet()

hltL3TrajectorySeed.tkSeedGenerator = "TSGFromPropagation"

hltL3TrajectorySeed.TSGFromPropagation = cms.PSet(
        ErrorRescaling = cms.double(3.0),
        ComponentName = cms.string('TSGFromPropagation'),
        errorMatrixPset = cms.PSet(),
        UpdateState = cms.bool(True),
        UseSecondMeasurements = cms.bool(False),
        SelectState = cms.bool(True),
        MaxChi2 = cms.double(15.0),
        UseVertexState = cms.bool(True),
        Propagator = cms.string('SmartPropagatorAnyOpposite'),
        ResetRescaling = cms.bool(True),
        SigmaZ = cms.double(25.0),       
 )
