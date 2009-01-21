import FWCore.ParameterSet.Config as cms

from regularL3_cff import hltL3Muons
from HLTrigger.Configuration.HLT_2E30_cff import hltL3TrajectorySeedOrig

import HLTrigger.Configuration.HLT_2E30_cff

hltTrajSeedIOHit   = hltL3TrajectorySeedOrig.clone()
hltTrajSeedOIState = hltL3TrajectorySeedOrig.clone()
hltTrajSeedOIHit   = hltL3TrajectorySeedOrig.clone()

#hltL3TrajectorySeed = cmsPSet()

hltTrajSeedOIHit.TSGForRoadSearchIOpxl = cms.PSet()
hltTrajSeedOIHit.TSGForRoadSearchOI = cms.PSet()
hltTrajSeedOIHit.MuonTrackingRegionBuilder = cms.PSet()
hltTrajSeedOIHit.TSGFromMixedPairs = cms.PSet()
hltTrajSeedOIHit.TSGFromPixelPairs = cms.PSet()
hltTrajSeedOIHit.TSGFromPixelTriplets = cms.PSet()
hltTrajSeedOIHit.TSGFromCombinedHits = cms.PSet()

hltTrajSeedOIHit.tkSeedGenerator = "TSGFromPropagation"

hltTrajSeedOIHit.TSGFromPropagation = cms.PSet(
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
#import FWCore.ParameterSet.Config as cms
#
#from regularL3_cff import hltL3Muons
#from HLTrigger.Configuration.HLT_2E30_cff import hltL3TrajectorySeed

#switch off a few things
hltTrajSeedOIState.TrackerSeedCleaner = cms.PSet()
hltTrajSeedOIState.TSGForRoadSearchIOpxl = cms.PSet()
hltTrajSeedOIState.TSGFromPropagation = cms.PSet()
hltTrajSeedOIState.MuonTrackingRegionBuilder = cms.PSet()
hltTrajSeedOIState.TSGFromMixedPairs = cms.PSet()
hltTrajSeedOIState.TSGFromPixelPairs = cms.PSet()
hltTrajSeedOIState.TSGFromPixelTriplets = cms.PSet()
hltTrajSeedOIState.TSGFromCombinedHits = cms.PSet()

#switch on the OIstate
hltTrajSeedOIState.tkSeedGenerator = "TSGForRoadSearchOI"

#from simpleRescale import *
#from ptDepRescale import *
from ptEtaRescale import *

hltTrajSeedOIState.TSGForRoadSearchOI.errorMatrixPset.action = "use"

hltTrajSeedOIState.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.xAxis = ptRange
hltTrajSeedOIState.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.yAxis = etaRange

hltTrajSeedOIState.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V11 = diagTerm
hltTrajSeedOIState.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V22 = diagTerm
hltTrajSeedOIState.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V33 = diagTerm
hltTrajSeedOIState.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V44 = diagTerm
hltTrajSeedOIState.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V55 = diagTerm

hltTrajSeedOIState.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V12 = offDiagTerm
hltTrajSeedOIState.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V13 = offDiagTerm
hltTrajSeedOIState.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V14 = offDiagTerm
hltTrajSeedOIState.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V15 = offDiagTerm
hltTrajSeedOIState.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V23 = offDiagTerm
hltTrajSeedOIState.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V24 = offDiagTerm
hltTrajSeedOIState.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V25 = offDiagTerm
hltTrajSeedOIState.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V34 = offDiagTerm
hltTrajSeedOIState.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V35 = offDiagTerm
hltTrajSeedOIState.TSGForRoadSearchOI.errorMatrixPset.errorMatrixValuesPSet.pf3_V45 = offDiagTerm

hltL3TrajectorySeed = cms.EDFilter("SeedCombiner",
    seedCollections = cms.VInputTag( 
        cms.InputTag("hltTrajSeedIOHit"),
        cms.InputTag("hltTrajSeedOIState"),
        cms.InputTag("hltTrajSeedOIHit"),
    )
)

comboSeeds_seq = cms.Sequence(
    hltTrajSeedIOHit
    +hltTrajSeedOIState
    +hltTrajSeedOIHit
    )

from HLTrigger.Configuration.HLT_2E30_cff import MuonCkfTrajectoryBuilder

#MuonCkfTrajectoryBuilder.useSeedLayer = True
#MuonCkfTrajectoryBuilder.rescaleErrorIfFail = 3.0

