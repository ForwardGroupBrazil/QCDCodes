import FWCore.ParameterSet.Config as cms

from regularL3_cff import hltL3Muons
from HLTrigger.Configuration.HLT_2E30_cff import hltL3TrajectorySeedOrig

hltL3TrajectorySeedOrig.TSGForRoadSearchIOpxl = cms.PSet()
hltL3TrajectorySeedOrig.TSGForRoadSearchOI = cms.PSet()
hltL3TrajectorySeedOrig.MuonTrackingRegionBuilder = cms.PSet()
hltL3TrajectorySeedOrig.TSGFromMixedPairs = cms.PSet()
hltL3TrajectorySeedOrig.TSGFromPixelPairs = cms.PSet()
hltL3TrajectorySeedOrig.TSGFromPixelTriplets = cms.PSet()
hltL3TrajectorySeedOrig.TSGFromCombinedHits = cms.PSet()

hltL3TrajectorySeedOrig.tkSeedGenerator = "TSGFromPropagation"

hltL3TrajectorySeedOrig.TSGFromPropagation = cms.PSet(
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

hltL3TrajectorySeed = cms.EDFilter("SeedCombiner",
    seedCollections = cms.VInputTag( 
        cms.InputTag("hltL3TrajectorySeedOrig"),
    )
)

comboSeeds_seq = cms.Sequence(
    hltL3TrajectorySeedOrig
    )
