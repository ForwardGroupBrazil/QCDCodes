import FWCore.ParameterSet.Config as cms

def PCut(process):
    process.hltL3TrajectorySeed.PCut = cms.double(2.5)

def cleanTSG(m):
    m.TSGForRoadSearchIOpxl = cms.PSet()
    m.TSGForRoadSearchOI = cms.PSet()
    m.TSGFromPropagation = cms.PSet()
    m.TSGFromMixedPairs = cms.PSet()
    m.TSGFromPixelPairs = cms.PSet()
    m.TSGFromPixelTriplets = cms.PSet()
    m.TSGFromCombinedHits = cms.PSet()
    

############ baseline ##############

def makeBaseline():
    return cms.PSet(
        ComponentName = cms.string('CombinedTSG'),
        PSetNames = cms.vstring('firstTSG','secondTSG'),

        firstTSG = cms.PSet(
           ComponentName = cms.string('TSGFromOrderedHits'),
           OrderedHitsFactoryPSet = cms.PSet(
             ComponentName = cms.string('StandardHitTripletGenerator'),
             SeedingLayers = cms.string('PixelLayerTriplets'),
             GeneratorPSet = cms.PSet(
              useBending = cms.bool(True),
              useFixedPreFiltering = cms.bool(False),
              phiPreFiltering = cms.double(0.3),
              extraHitRPhitolerance = cms.double(0.06),
              useMultScattering = cms.bool(True),
              ComponentName = cms.string('PixelTripletHLTGenerator'),
              extraHitRZtolerance = cms.double(0.06)
              )
             ),
           TTRHBuilder = cms.string('WithTrackAngle')
           ),

        secondTSG = cms.PSet(
            ComponentName = cms.string('TSGFromOrderedHits'),
            OrderedHitsFactoryPSet = cms.PSet(
             ComponentName = cms.string('StandardHitPairGenerator'),
             SeedingLayers = cms.string('PixelLayerPairs')
             ),
            TTRHBuilder = cms.string('WithTrackAngle')
            )
        )


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
        
def SwitchToBaseline(process):
    PCut(process)
    cleanTSG(process.hltL3TrajectorySeed)
    process.hltL3TrajectorySeed.TSGFromCombinedHits = makeBaseline()
    process.hltL3TrajectorySeed.tkSeedGenerator = "TSGFromCombinedHits"
    process.hltL3TrajectorySeed.ServiceParameters.Propagators = cms.untracked.vstring()
    process.hltL3TrajectorySeed.MuonTrackingRegionBuilder = regionBuilder()
    process.hltL3TrajectorySeed.TrackerSeedCleaner = seedCleaner()
    

def makeBaselinePP():
    pset=makeBaseline()
    pset.PSetNames.append('thirdTSG')
    pset.thirdTSG = cms.PSet(
           ComponentName = cms.string('DualByEtaTSG'),
           PSetNames = cms.vstring('endcapTSG','barrelTSG'),
           barrelTSG = cms.PSet(    ),
           endcapTSG = cms.PSet(
            ComponentName = cms.string('TSGFromOrderedHits'),
            OrderedHitsFactoryPSet = cms.PSet(
              ComponentName = cms.string('StandardHitPairGenerator'),
              SeedingLayers = cms.string('MixedLayerPairs')
              ),
            TTRHBuilder = cms.string('WithTrackAngle')
            ),
           etaSeparation = cms.double(2.0)
          )
    return pset


############# baseline plus mixed seeding ######
    
def SwitchToBaselinePP(process):
    PCut(process)
    cleanTSG(process.hltL3TrajectorySeed)
    process.hltL3TrajectorySeed.TSGFromCombinedHits = makeBaselinePP()
    process.hltL3TrajectorySeed.tkSeedGenerator = "TSGFromCombinedHits"
    process.hltL3TrajectorySeed.ServiceParameters.Propagators = cms.untracked.vstring()
    process.hltL3TrajectorySeed.MuonTrackingRegionBuilder = regionBuilder()
    process.hltL3TrajectorySeed.TrackerSeedCleaner = seedCleaner()


############### OI-state based #############
    
def makeOIState():
    from UserCode.L3Switches.ptEtaRescale import ptRange,etaRange,diagTerm,offDiagTerm

    return  cms.PSet(
        ComponentName = cms.string('TSGForRoadSearch'),
        option = cms.uint32(3),
        propagatorCompatibleName = cms.string('SteppingHelixPropagatorOpposite'),
        propagatorName = cms.string('SteppingHelixPropagatorAlong'),
        maxChi2 = cms.double(40.0),
        manySeeds = cms.bool(False),
        copyMuonRecHit = cms.bool(False),
        errorMatrixPset = cms.PSet(
                   action = cms.string('use'),
                   atIP = cms.bool(True),
                   errorMatrixValuesPSet = cms.PSet(
                       xAxis = ptRange,
                       yAxis = etaRange,
                       zAxis = cms.vdouble(-3.14159, 3.14159),
                       
                       pf3_V11 = diagTerm,
                       pf3_V22 = diagTerm,
                       pf3_V33 = diagTerm,
                       pf3_V44 = diagTerm,
                       pf3_V55 = diagTerm,

                       pf3_V12 = offDiagTerm,
                       pf3_V13 = offDiagTerm,
                       pf3_V14 = offDiagTerm,
                       pf3_V15 = offDiagTerm,
                       pf3_V23 = offDiagTerm,
                       pf3_V24 = offDiagTerm,
                       pf3_V25 = offDiagTerm,
                       pf3_V34 = offDiagTerm,
                       pf3_V35 = offDiagTerm,
                       pf3_V45 = offDiagTerm
                       )
                   )
        )

def OIStatePropagators(process,pset):
    if (not hasattr(process.hltL3TrajectorySeed.ServiceParameters,"Propagators")):
        process.hltL3TrajectorySeed.ServiceParameters.Propagators = cms.untracked.vstring()
    process.hltL3TrajectorySeed.ServiceParameters.Propagators.append(pset.propagatorCompatibleName.value())
    process.hltL3TrajectorySeed.ServiceParameters.Propagators.append(pset.propagatorName.value())

    
def SwitchToOIState(process):
    PCut(process)
    #switch off a few things
    cleanTSG(process.hltL3TrajectorySeed)

    process.hltL3TrajectorySeed.MuonTrackingRegionBuilder = cms.PSet()
    process.hltL3TrajectorySeed.TrackerSeedCleaner = cms.PSet()

    #    process.hltL3TrajectorySeed
    
    #switch on the OIstate
    process.hltL3TrajectorySeed.tkSeedGenerator = "TSGForRoadSearchOI"
    process.hltL3TrajectorySeed.TSGForRoadSearchOI = makeOIState()
    process.hltL3TrajectorySeed.ServiceParameters.Propagators = cms.untracked.vstring()
    OIStatePropagators(process,process.hltL3TrajectorySeed.TSGForRoadSearchOI)



########## OI hit-based ###########

def makeOIHit():
    from UserCode.L3Switches.ptEtaRescale import ptRange,etaRange,diagTerm,offDiagTerm
    
    return cms.PSet(
        ResetMethod = cms.string('matrix'),
        ErrorRescaling = cms.double(3.0),
        ComponentName = cms.string('TSGFromPropagation'),
        errorMatrixPset = cms.PSet(
                   action = cms.string('use'),
                   atIP = cms.bool(True),
                   errorMatrixValuesPSet = cms.PSet(
                       xAxis = ptRange,
                       yAxis = etaRange,
                       zAxis = cms.vdouble(-3.14159, 3.14159),
                       
                       pf3_V11 = diagTerm,
                       pf3_V22 = diagTerm,
                       pf3_V33 = diagTerm,
                       pf3_V44 = diagTerm,
                       pf3_V55 = diagTerm,
                       
                       pf3_V12 = offDiagTerm,
                       pf3_V13 = offDiagTerm,
                       pf3_V14 = offDiagTerm,
                       pf3_V15 = offDiagTerm,
                       pf3_V23 = offDiagTerm,
                       pf3_V24 = offDiagTerm,
                       pf3_V25 = offDiagTerm,
                       pf3_V34 = offDiagTerm,
                       pf3_V35 = offDiagTerm,
                       pf3_V45 = offDiagTerm
                       )
                   ),
        UpdateState = cms.bool(True),
        SelectState = cms.bool(False),
        MaxChi2 = cms.double(40.0),
        UseVertexState = cms.bool(True),
        Propagator = cms.string('SmartPropagatorAnyOpposite'),
        SigmaZ = cms.double(25.0),
        )

def OIHitPropagators(process,pset):
    if (not hasattr(process.hltL3TrajectorySeed.ServiceParameters,"Propagators")):
        process.hltL3TrajectorySeed.ServiceParameters.Propagators = cms.untracked.vstring()
    process.hltL3TrajectorySeed.ServiceParameters.Propagators.append('PropagatorWithMaterial')
    process.hltL3TrajectorySeed.ServiceParameters.Propagators.append(pset.Propagator.value())

def SwitchToOIHit(process):
    PCut(process)
    cleanTSG(process.hltL3TrajectorySeed)
    process.hltL3TrajectorySeed.MuonTrackingRegionBuilder = cms.PSet()

    process.hltL3TrajectorySeed.tkSeedGenerator = "TSGFromPropagation"

    process.hltL3TrajectorySeed.TSGFromPropagation = makeOIHit()
    process.hltL3TrajectorySeed.ServiceParameters.Propagators = cms.untracked.vstring()
    OIHitPropagators(process,process.hltL3TrajectorySeed.TSGFromPropagation)

################## all combined in one module


def SwitchToAllCombined(process):
    cleanTSG(process.hltL3TrajectorySeed)
    PCut(process)
    
    process.hltL3TrajectorySeed.TSGFromCombinedSeeds = cms.PSet(
        ComponentName = cms.string("CombinedTSG"),
        PSetNames = cms.vstring('baselinePP','oiState','oiHit'),
        baselinePP = makeBaselinePP(),
        oiState= makeOIState(),
        oiHit= makeOIHit()
        )
    process.hltL3TrajectorySeed.tkSeedGenerator = "TSGFromCombinedSeeds"

    OIStatePropagators(process,process.hltL3TrajectorySeed.TSGFromCombinedSeeds.oiState)
    OIHitPropagators(process,process.hltL3TrajectorySeed.TSGFromCombinedSeeds.oiHit)

def SwitchToOICombined(process):
    cleanTSG(process.hltL3TrajectorySeed)
    PCut(process)
    
    process.hltL3TrajectorySeed.TSGFromCombinedSeeds = cms.PSet(
        ComponentName = cms.string("CombinedTSG"),
        PSetNames = cms.vstring('oiState','oiHit'),
        oiState= makeOIState(),
        oiHit= makeOIHit()
        )
    process.hltL3TrajectorySeed.tkSeedGenerator = "TSGFromCombinedSeeds"

    OIStatePropagators(process,process.hltL3TrajectorySeed.TSGFromCombinedSeeds.oiState)
    OIHitPropagators(process,process.hltL3TrajectorySeed.TSGFromCombinedSeeds.oiHit)

############### all combined from different modules
    
def SwitchToComboSeeds(process):
    cleanTSG(process.hltL3TrajectorySeed)
    PCut(process)
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
    process.hltL3TrajSeedOIHit.TSGFromPropagation = makeDualByIterativeOIHit()
    process.hltL3TrajSeedOIHit.tkSeedGenerator = 'TSGFromPropagation'
    SwitchToBaseline(process)
    process.hltL3TrajSeedIOHit = process.hltL3TrajectorySeed.clone()
    process.hltL3TrajSeedIOHit.TSGFromCombinedHits = makeDualByIterative()
    process.hltL3TrajSeedIOHit.tkSeedGenerator = 'TSGFromCombinedHits'
    
    process.hltL3TrajectorySeed = process.l3SeedCombination
#    process.HLTL3muonrecoSequence.replace(process.hltL3TrajectorySeed, process.hltTrajSeedIOHit + process.hltTrajSeedOIState + process.hltTrajSeedOIHit + process.hltL3TrajectorySeed)

################ iterative tracking version #########

def makeDualByIterative():
##    return makeBaseline()
    return cms.PSet (
        ComponentName = cms.string('DualByL2TSG'),
        PSetNames = cms.vstring('skipTSG','iterativeTSG'),
        skipTSG = cms.PSet(    ),
        iterativeTSG = makeBaseline(),
        L3TkCollectionA = cms.InputTag('l3TkFromL2OICombination'),
        )

def makeDualByIterativeOIHit():
##    return makeBaseline()
    return cms.PSet (
        ComponentName = cms.string('DualByL2TSG'),
        PSetNames = cms.vstring('skipTSG','iterativeTSG'),
        skipTSG = cms.PSet(    ),
        iterativeTSG = makeOIHit(),
        #aaa   L3TkCollectionA = cms.InputTag('hltL3TkTracksFromL2OIState'),
        L3TkCollectionA = cms.InputTag('hltL3MuonsOIState'),
        )

def SwitchToIterative(process):
    cleanTSG(process.hltL3TrajectorySeed)
    PCut(process)
    SwitchToOICombined(process)
        
    process.hltL3TrajSeedIOHit = process.hltL3TrajectorySeed.clone()
    cleanTSG(process.hltL3TrajSeedIOHit)
    process.hltL3TrajSeedIOHit.TSGFromCombinedSeeds = cms.PSet ( )
    process.hltL3TrajSeedIOHit.TSGFromCombinedHits = makeDualByIterative()
    #process.hltL3TrajSeedIOHit.TSGFromCombinedHits = makeBaseline()
    process.hltL3TrajSeedIOHit.tkSeedGenerator = "TSGFromCombinedHits"
    process.hltL3TrajSeedIOHit.ServiceParameters.Propagators = cms.untracked.vstring()
    process.hltL3TrajSeedIOHit.MuonTrackingRegionBuilder = regionBuilder()
    process.hltL3TrajSeedIOHit.TrackerSeedCleaner = seedCleaner()
    
    process.hltL3TrackCandidateFromL2IO = process.hltL3TrackCandidateFromL2.clone()
    process.hltL3TrackCandidateFromL2IO.src = "hltL3TrajSeedIOHit"
    
    process.hltL3TkTracksFromL2OI = process.hltL3TkTracksFromL2.clone()
    process.hltL3TkTracksFromL2OI.src = "hltL3TrackCandidateFromL2"
    process.hltL3TkTracksFromL2IO = process.hltL3TkTracksFromL2.clone()
    process.hltL3TkTracksFromL2IO.src = "hltL3TrackCandidateFromL2IO"

    process.l3TkFromL2OICombination = cms.EDProducer(
        "L3TrackCombiner",
        labels = cms.VInputTag(
        cms.InputTag("hltL3TkTracksFromL2OI"),
        )
        )


    process.l3TkFromL2Combination = cms.EDProducer(
        "L3TrackCombiner",
        labels = cms.VInputTag(
        cms.InputTag("hltL3TkTracksFromL2OI"),
        cms.InputTag("hltL3TkTracksFromL2IO")
        )
        )

    process.hltL3TkTracksFromL2 = process.l3TkFromL2Combination

    process.HLTL3muonrecoNocandSequence.replace(process.hltL3TkTracksFromL2,process.hltL3TkTracksFromL2OI * process.l3TkFromL2OICombination * process.hltL3TrajSeedIOHit + process.hltL3TrackCandidateFromL2IO + process.hltL3TkTracksFromL2IO + process.hltL3TkTracksFromL2)


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
    
