import FWCore.ParameterSet.Config as cms

process = cms.Process("MuonPlots")

# Messages
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.GlobalTag.globaltag = 'START3X_V25B::All'

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/data/gpetrucc/7TeV/inclusiveMuons/inclusiveMuons_MuHLT_MinBiasMC357_v3/tupleMC_1_1.root',
        'file:/data/gpetrucc/7TeV/inclusiveMuons/inclusiveMuons_MuHLT_MinBiasMC357_v3/tupleMC_2_1.root',
        'file:/data/gpetrucc/7TeV/inclusiveMuons/inclusiveMuons_MuHLT_MinBiasMC357_v3/tupleMC_3_1.root',
        'file:/data/gpetrucc/7TeV/inclusiveMuons/inclusiveMuons_MuHLT_MinBiasMC357_v3/tupleMC_4_1.root',
        'file:/data/gpetrucc/7TeV/inclusiveMuons/inclusiveMuons_MuHLT_MinBiasMC357_v3/tupleMC_5_1.root',
        'file:/data/gpetrucc/7TeV/inclusiveMuons/inclusiveMuons_MuHLT_MinBiasMC357_v3/tupleMC_6_1.root',
        'file:/data/gpetrucc/7TeV/inclusiveMuons/inclusiveMuons_MuHLT_MinBiasMC357_v3/tupleMC_7_1.root',
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.TFileService = cms.Service('TFileService',
    fileName=cms.string('inclusiveMuonPlots_MC.class.root')
)

from MuonAnalysis.Examples.inclusiveMuonPlots_cfi import makeInclusiveMuonPlots;
commonInputs = cms.PSet(
    muons     = cms.InputTag('patMuonsWithTrigger'),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
)
process.trackerMuons = cms.EDAnalyzer("InclusiveMuonPlots",
    makeInclusiveMuonPlots(),
    commonInputs,
    selection = cms.string("isTrackerMuon && muonID('TMLastStationAngTight')"),
)
process.globalMuons = process.trackerMuons.clone(
    selection = "isGlobalMuon"
)
process.standAloneMuons = process.trackerMuons.clone(
    selection = "isStandAloneMuon"
)


from HLTrigger.HLTfilters.hltHighLevelDev_cfi import hltHighLevelDev
reskim  = hltHighLevelDev.clone(TriggerResultsTag = cms.InputTag('TriggerResults','',''))
#process.bscMinBias = hltHighLevelDev.clone(HLTPaths = ['HLT_MinBiasBSC'], HLTPathsPrescales = [1])
process.bscMinBias = hltHighLevelDev.clone(HLTPaths = ['HLT_L1_BscMinBiasOR_BptxPlusORMinus'], HLTPathsPrescales = [1])
process.bit40    = reskim.clone(HLTPaths = ['Flag_Bit40'],    HLTPathsPrescales = [1])
process.haloVeto = reskim.clone(HLTPaths = ['Flag_HaloVeto'], HLTPathsPrescales = [1])
process.preFilter = cms.Sequence(process.bit40 + process.haloVeto + process.bscMinBias)

process.p = cms.Path(process.preFilter * process.trackerMuons * process.globalMuons * process.standAloneMuons)

def addClassAliases(process, label, classifier):
    module = getattr(process,label)
    sel = module.selection.value()
    modules = []
    classes = (("Ghost", "userInt('%s') <= -2"), ("Punch", "abs(userInt('%s')) <= 1"), ("Light", "userInt('%s') == 2"), ("Heavy", "userInt('%s') == 3")) 
    for X,C in classes:
        setattr(process, label+X, module.clone(selection = sel + ' && ' + (C % classifier)))
        modules.append( getattr(process, label+X) )
    setattr(process, label+"_Classes", cms.Path(sum(modules,process.preFilter)))

if True:
    addClassAliases(process, "trackerMuons",    "classByHitsTM")
    addClassAliases(process, "globalMuons",     "classByHitsGlb")
    addClassAliases(process, "standAloneMuons", "classByHitsSta")

if True:
    process.standAloneMuonsValidHits = process.standAloneMuons.clone(selection = "isStandAloneMuon && outerTrack.numberOfValidHits > 0")
    process.standAloneMuonsValidHits_Path = cms.Path(process.preFilter + process.standAloneMuonsValidHits)
    addClassAliases(process, "standAloneMuonsValidHits", "classByHitsSta")
    if True:
        for M in ('standAloneMuons', 'standAloneMuonsValidHits'):
            for X,C in (("B"," && abs(eta) < 1.1"), ("E", "&& abs(eta) > 1.1")):
                base = getattr(process,M)
                setattr(process, M+X, base.clone(selection = base.selection.value() + C))
                setattr(process, M+X+"_Path", cms.Path(process.preFilter + getattr(process,M+X)))
                addClassAliases(process, M+X, "classByHitsSta")
                    
if True:
    for X,C in (("Barrel", "abs(eta) <= 0.8"), ("Transition", "0.8 < abs(eta) <= 1.1"), ("Forward", "1.1 < abs(eta) <= 1.6"), ("VeryForward", "1.6 < abs(eta)")):
        setattr(process, "trackerMuons"+X, process.trackerMuons.clone(
            selection = "isTrackerMuon && muonID('TMLastStationAngTight') && "+C
        ))
        setattr(process, "pTM"+X, cms.Path(process.preFilter + getattr(process,"trackerMuons"+X)))
        if True:
            addClassAliases(process, "trackerMuons"+X, "classByHitsTM")
if True:
    type1 = " isGlobalMuon &&  muonID('GlobalMuonPromptTight')"
    type2 = type1 + "&& isTrackerMuon && track.numberOfValidHits > 10 & track.hitPattern.pixelLayersWithMeasurement >= 1"
    type3 = type2 + "&&  (userInt('muonStations:csc') >= 2 || userInt('muonStations:dt') >= 2)"
    type4 = type3 + "&&  (track.trackerExpectedHitsInner.numberOfLostHits <=1 && track.trackerExpectedHitsOuter.numberOfLostHits <= 1)"
    process.globalMuonsType1 = process.globalMuons.clone(selection = type1)
    process.globalMuonsType1B = process.globalMuonsType1.clone(selection = type1+" && abs(eta) < 1.0")
    process.globalMuonsType1E = process.globalMuonsType1.clone(selection = type1+" && abs(eta) > 1.0")
    process.globalMuonsType2 = process.globalMuons.clone(selection = type2)
    process.globalMuonsType2B = process.globalMuonsType2.clone(selection = type2+" && abs(eta) < 1.0")
    process.globalMuonsType2E = process.globalMuonsType2.clone(selection = type2+" && abs(eta) > 1.0")
    process.globalMuonsType3 = process.globalMuons.clone(selection = type3)
    process.globalMuonsType3B = process.globalMuonsType3.clone(selection = type3+" && abs(eta) < 1.0")
    process.globalMuonsType3E = process.globalMuonsType3.clone(selection = type3+" && abs(eta) > 1.0")
    process.globalMuonsType4 = process.globalMuons.clone(selection = type4)
    process.globalMuonsType4B = process.globalMuonsType4.clone(selection = type4+" && abs(eta) < 1.0")
    process.globalMuonsType4E = process.globalMuonsType4.clone(selection = type4+" && abs(eta) > 1.0")
    process.typedGlobalMuons = cms.Sequence(
        process.globalMuonsType1 + process.globalMuonsType1B + process.globalMuonsType1E +
        process.globalMuonsType2 + process.globalMuonsType2B + process.globalMuonsType2E +
        process.globalMuonsType3 + process.globalMuonsType3B + process.globalMuonsType3E +
        process.globalMuonsType4 + process.globalMuonsType4B + process.globalMuonsType4E 
    )
    process.p2 = cms.Path(process.preFilter + process.typedGlobalMuons)
    if True:
        for I in [1,2,3,4]:
            for P in ["", "B", "E"]:
                addClassAliases(process, "globalMuonsType%d%s" % (I,P), "classByHitsGlb")

if True:
    process.globalMuonsB = process.globalMuons.clone(selection = "isGlobalMuon && abs(eta) < 1.0") 
    process.globalMuonsE = process.globalMuons.clone(selection = "isGlobalMuon && abs(eta) > 1.0") 
    process.pBE = cms.Path(process.preFilter + process.globalMuonsB + process.globalMuonsE)
    if True:
        addClassAliases(process, "globalMuonsB", "classByHitsGlb")
        addClassAliases(process, "globalMuonsE", "classByHitsGlb")

if True:
    process.globalMuonsPt5 = cms.EDAnalyzer("InclusiveMuonPlots",
        makeInclusiveMuonPlots(2),
        commonInputs,
        selection = cms.string("isGlobalMuon && pt > 5"),
    )
    process.globalMuonsPt5B = process.globalMuonsPt5.clone(selection = "isGlobalMuon && pt > 5 && abs(eta) < 1.0") 
    process.globalMuonsPt5E = process.globalMuonsPt5.clone(selection = "isGlobalMuon && pt > 5 && abs(eta) > 1.0") 
    process.p5  = cms.Path(process.preFilter + process.globalMuonsPt5 + process.globalMuonsPt5B + process.globalMuonsPt5E)
    if True:
        for P in ["", "B", "E"]: addClassAliases(process, "globalMuonsPt5%s" % P, "classByHitsGlb")
if False:
    process.goodTrackGlobalMuonsPt5 = process.globalMuonsPt5.clone(
        selection = "isGlobalMuon && pt > 5 && track.trackerExpectedHitsOuter.numberOfLostHits <= 3 && track.hitPattern.pixelLayersWithMeasurement >= 2 && track.normalizedChi2 < 5",
    )
    process.goodMatchGlobalMuonsPt5 = process.globalMuonsPt5.clone(
        selection = "isGlobalMuon && pt > 5 && muonID('GlobalMuonPromptTight') && muonID('TMLastStationLoose')",
    )
    process.isolatedGlobalMuonsPt5  = process.globalMuonsPt5.clone(
        selection = process.globalMuonsPt5.selection.value() + " && ( isolationR03.hadEt + isolationR03.emEt + isolationR03.sumPt ) / pt < 0.3"
    )
    process.p.replace(process.globalMuons, process.globalMuons + 
        process.globalMuonsPt5 + process.globalMuonsPt5B + process.globalMuonsPt5E +
        process.goodTrackGlobalMuonsPt5 + process.goodMatchGlobalMuonsPt5 + process.isolatedGlobalMuonsPt5
    )

if True:
    process.TFileService.fileName = 'inclusiveMuonPlots_MC.Mu3.root'
    process.hltMu3 = hltHighLevelDev.clone(HLTPaths = ['HLT_Mu3'], HLTPathsPrescales = [1])
    process.preFilter.replace(process.bscMinBias, process.hltMu3)
    for N in process.analyzerNames().split():
        A = getattr(process,N)
        if A.type_() == 'InclusiveMuonPlots':
            A.selection = A.selection.value()+' && !triggerObjectMatchesByFilter("hltSingleMu3L3Filtered3").empty()';



