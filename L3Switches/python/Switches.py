import FWCore.ParameterSet.Config as cms

def regionBuilder():
    return cms.PSet(EtaR_UpperLimit_Par1 = cms.double( 0.25 ),
                    Eta_fixed = cms.double( 0.2 ),
                    beamSpot = cms.InputTag( "hltOfflineBeamSpot" ),
                    OnDemand = cms.double( -1.0 ),
                    Rescale_Dz = cms.double( 3.0 ),
                    Eta_min = cms.double( 0.1 ),
                    Rescale_phi = cms.double( 3.0 ),
                    PhiR_UpperLimit_Par1 = cms.double( 0.6 ),
                    DeltaZ_Region = cms.double( 15.9 ),
                    Phi_min = cms.double( 0.1 ),
                    PhiR_UpperLimit_Par2 = cms.double( 0.2 ),
                    vertexCollection = cms.InputTag( "pixelVertices" ),
                    Phi_fixed = cms.double( 0.2 ),
                    DeltaR = cms.double( 0.2 ),
                    EtaR_UpperLimit_Par2 = cms.double( 0.15 ),
                    UseFixedRegion = cms.bool( False ),
                    Rescale_eta = cms.double( 3.0 ),
                    UseVertex = cms.bool( False ),
                    EscapePt = cms.double( 1.5 )
                    )

def seedCleaner():
    return cms.PSet(cleanerFromSharedHits = cms.bool( True ),
                    ptCleaner = cms.bool( True ),
                    TTRHBuilder = cms.string( "WithTrackAngle" ),
                    beamSpot = cms.InputTag( "hltOfflineBeamSpot" ),
                    directionCleaner = cms.bool( True )
                    )

############# baseline pixel pair and pixel triplet seeding #####
        
def SwitchToBaseline(process):
    from RecoMuon.TrackerSeedGenerator.TSGs_cff import TSGsBlock
    process.hltL3TrajectorySeed.TkSeedGenerator = TSGsBlock.TSGFromCombinedHits
    process.hltL3TrajectorySeed.ServiceParameters.Propagators = cms.untracked.vstring()
    process.hltL3TrajectorySeed.MuonTrackingRegionBuilder = regionBuilder()
    process.hltL3TrajectorySeed.TrackerSeedCleaner = seedCleaner()
    
############# baseline plus mixed seeding ######
    
def SwitchToBaselinePP(process):
    from RecoMuon.TrackerSeedGenerator.TSGs_cff import TSGsBlock
    process.hltL3TrajectorySeed.TkSeedGenerator = TSGsBlock.TSGFromCombinedHits
    process.hltL3TrajectorySeed.TkSeedGenerator.PSetNames.append('thirdTSG')
    process.hltL3TrajectorySeed.ServiceParameters.Propagators = cms.untracked.vstring()
    process.hltL3TrajectorySeed.MuonTrackingRegionBuilder = regionBuilder()
    process.hltL3TrajectorySeed.TrackerSeedCleaner = seedCleaner()


############### OI-state based #############
    
def OIStatePropagators(process,pset):
    if (not hasattr(process.hltL3TrajectorySeed.ServiceParameters,"Propagators")):
        process.hltL3TrajectorySeed.ServiceParameters.Propagators = cms.untracked.vstring()
    process.hltL3TrajectorySeed.ServiceParameters.Propagators.append(pset.propagatorCompatibleName.value())
    process.hltL3TrajectorySeed.ServiceParameters.Propagators.append(pset.propagatorName.value())

    
def SwitchToOIState(process):
    from RecoMuon.TrackerSeedGenerator.TSGs_cff import TSGsBlock
    process.hltL3TrajectorySeed.MuonTrackingRegionBuilder = cms.PSet()
    process.hltL3TrajectorySeed.TrackerSeedCleaner = cms.PSet()
    process.hltL3TrajectorySeed.TkSeedGenerator = TSGsBlock.TSGForRoadSearchOI
    process.hltL3TrajectorySeed.ServiceParameters.Propagators = cms.untracked.vstring()
    OIStatePropagators(process,process.hltL3TrajectorySeed.TkSeedGenerator)

########## OI hit-based ###########

def OIHitPropagators(process,pset):
    if (not hasattr(process.hltL3TrajectorySeed.ServiceParameters,"Propagators")):
        process.hltL3TrajectorySeed.ServiceParameters.Propagators = cms.untracked.vstring()
    process.hltL3TrajectorySeed.ServiceParameters.Propagators.append('PropagatorWithMaterial')
    process.hltL3TrajectorySeed.ServiceParameters.Propagators.append(pset.Propagator.value())

def SwitchToOIHit(process):
    from RecoMuon.TrackerSeedGenerator.TSGs_cff import TSGsBlock
    process.hltL3TrajectorySeed.MuonTrackingRegionBuilder = cms.PSet()
    process.hltL3TrajectorySeed.TrackerSeedCleaner = seedCleaner()
    process.hltL3TrajectorySeed.TkSeedGenerator = TSGsBlock.TSGFromPropagation
    process.hltL3TrajectorySeed.ServiceParameters.Propagators = cms.untracked.vstring()
    OIHitPropagators(process,process.hltL3TrajectorySeed.TkSeedGenerator)

################## all OI combined in one module

def SwitchToOICombined(process):
    from RecoMuon.TrackerSeedGenerator.TSGs_cff import TSGsBlock
    process.hltL3TrajectorySeed.TkSeedGenerator = cms.PSet(
        ComponentName = cms.string("CombinedTSG"),
        PSetNames = cms.vstring('oiState','oiHit'),
        oiState= TSGsBlock.TSGForRoadSearchOI,
        oiHit= TSGsBlock.TSGFromPropagation
        )

    OIStatePropagators(process,process.hltL3TrajectorySeed.TkSeedGenerator.oiState)
    OIHitPropagators(process,process.hltL3TrajectorySeed.TkSeedGenerator.oiHit)

############### all combined from different modules
    
def SwitchToComboSeeds(process):
    from RecoMuon.TrackerSeedGenerator.TSGs_cff import TSGsBlock
    process.l3SeedCombination = cms.EDProducer(
        "L3MuonTrajectorySeedCombiner",
        labels = cms.VInputTag(
           cms.InputTag("hltL3TrajSeedIOHit"),
           cms.InputTag("hltL3TrajSeedOIState"),
           cms.InputTag("hltL3TrajSeedOIHit")
          )
        )
    SwitchToOIState(process)
    process.hltL3TrajSeedOIState = process.hltL3TrajectorySeed.clone()
    SwitchToOIHit(process)
    process.hltL3TrajSeedOIHit = process.hltL3TrajectorySeed.clone()
    process.hltL3TrajSeedOIHit.TkSeedGenerator = makeDualByIterativeOIHit()
    SwitchToBaseline(process)
    process.hltL3TrajSeedIOHit = process.hltL3TrajectorySeed.clone()
    process.hltL3TrajSeedIOHit.TkSeedGenerator = makeDualByIterative()
    
    process.hltL3TrajectorySeed = process.l3SeedCombination


################ iterative tracking version #########

def makeDualByIterative():
##    return makeBaseline()
    from RecoMuon.TrackerSeedGenerator.TSGs_cff import TSGsBlock
    return cms.PSet (
        ComponentName = cms.string('DualByL2TSG'),
        PSetNames = cms.vstring('skipTSG','iterativeTSG'),
        skipTSG = cms.PSet(    ),
        iterativeTSG = TSGsBlock.TSGFromCombinedHits,
        L3TkCollectionA = cms.InputTag('l3TkFromL2OICombination'),
        )

def makeDualByIterativeOIHit():
##    return makeBaseline()
    from RecoMuon.TrackerSeedGenerator.TSGs_cff import TSGsBlock
    return cms.PSet (
        ComponentName = cms.string('DualByL2TSG'),
        PSetNames = cms.vstring('skipTSG','iterativeTSG'),
        skipTSG = cms.PSet(    ),
        iterativeTSG = TSGsBlock.TSGFromPropagation,
        #aaa   L3TkCollectionA = cms.InputTag('hltL3TkTracksFromL2OIState'),
        L3TkCollectionA = cms.InputTag('hltL3MuonsOIState'),
        )



def SwitchToIterative3(process):
    SwitchToComboSeeds(process)
        
    process.hltL3TrackCandidateFromL2IOHit = process.hltL3TrackCandidateFromL2.clone()
    process.hltL3TrackCandidateFromL2IOHit.src = "hltL3TrajSeedIOHit"

    process.hltL3TrackCandidateFromL2OIHit = process.hltL3TrackCandidateFromL2.clone()
    process.hltL3TrackCandidateFromL2OIHit.src = "hltL3TrajSeedOIHit"

    process.hltL3TrackCandidateFromL2OIState = process.hltL3TrackCandidateFromL2.clone()
    process.hltL3TrackCandidateFromL2OIState.src = "hltL3TrajSeedOIState"

    process.l3TkCandFromL2Combination = cms.EDProducer(
        "L3TrackCandCombiner",
        labels = cms.VInputTag(
        cms.InputTag("hltL3TrackCandidateFromL2IOHit"),
        cms.InputTag("hltL3TrackCandidateFromL2OIHit"),
        cms.InputTag("hltL3TrackCandidateFromL2OIState"),
        )
        )
    
    process.hltL3TkTracksFromL2IOHit = process.hltL3TkTracksFromL2.clone()
    process.hltL3TkTracksFromL2IOHit.src = "hltL3TrackCandidateFromL2IOHit"

    process.hltL3TkTracksFromL2OIHit = process.hltL3TkTracksFromL2.clone()
    process.hltL3TkTracksFromL2OIHit.src = "hltL3TrackCandidateFromL2OIHit"
    
    process.hltL3TkTracksFromL2OIState = process.hltL3TkTracksFromL2.clone()
    process.hltL3TkTracksFromL2OIState.src = "hltL3TrackCandidateFromL2OIState"

    process.hltL3MuonsOIState = process.hltL3Muons.clone()
    process.hltL3MuonsOIState.L3TrajBuilderParameters.tkTrajLabel = "hltL3TkTracksFromL2OIState"

    process.hltL3MuonsOIHit = process.hltL3Muons.clone()
    process.hltL3MuonsOIHit.L3TrajBuilderParameters.tkTrajLabel = "hltL3TkTracksFromL2OIHit"

    process.hltL3MuonsIOHit = process.hltL3Muons.clone()
    process.hltL3MuonsIOHit.L3TrajBuilderParameters.tkTrajLabel = "hltL3TkTracksFromL2IOHit"

    process.l3TkFromL2OICombination = cms.EDProducer(
        "L3TrackCombiner",
        labels = cms.VInputTag(
        #aaa cms.InputTag("hltL3TkTracksFromL2OIHit"),
        #aaa cms.InputTag("hltL3TkTracksFromL2OIState")
        cms.InputTag("hltL3MuonsOIState"),
        cms.InputTag("hltL3MuonsOIHit"),
        )
        )
    
    process.l3TkFromL2Combination = cms.EDProducer(
        "L3TrackCombiner",
        labels = cms.VInputTag(
        cms.InputTag("hltL3TkTracksFromL2IOHit"),
        cms.InputTag("hltL3TkTracksFromL2OIHit"),
        cms.InputTag("hltL3TkTracksFromL2OIState")
        )
        )

    process.l3MuonsCombination = cms.EDProducer(
        "L3TrackCombiner",
        labels = cms.VInputTag(
        cms.InputTag("hltL3MuonsOIState"),
        cms.InputTag("hltL3MuonsOIHit"),
        cms.InputTag("hltL3MuonsIOHit")
        )
        )

    process.hltL3TkTracksFromL2 = process.l3TkFromL2Combination
    process.hltL3TrackCandidateFromL2 = process.l3TkCandFromL2Combination
    process.hltL3Muons = process.l3MuonsCombination

    process.HLTL3muonTkCandidateSequence = cms.Sequence(
         process.HLTDoLocalPixelSequence +
         process.HLTDoLocalStripSequence +
         #
         process.hltL3TrajSeedOIState +
         process.hltL3TrackCandidateFromL2OIState +
         process.hltL3TkTracksFromL2OIState +
         #
         process.hltL3MuonsOIState +
         #
         process.hltL3TrajSeedOIHit +
         process.hltL3TrackCandidateFromL2OIHit +
         process.hltL3TkTracksFromL2OIHit +        
         #
         process.hltL3MuonsOIHit +
         #
         process.l3TkFromL2OICombination + ######
         #
         process.hltL3TrajSeedIOHit +
         process.hltL3TrackCandidateFromL2IOHit +
         process.hltL3TkTracksFromL2IOHit +
         #
         process.hltL3MuonsIOHit +    #????
         #
         process.hltL3TrajectorySeed + ######
         process.hltL3TrackCandidateFromL2 ######
         )

    process.HLTL3muonrecoNocandSequence = cms.Sequence(
        process.HLTL3muonTkCandidateSequence +
        process.hltL3TkTracksFromL2 +  ######
        process.hltL3Muons ######
        )
    
