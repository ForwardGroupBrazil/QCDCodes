import FWCore.ParameterSet.Config as cms

process = cms.Process("MuonPlots")

# Messages
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#process.MessageLogger.cerr.FwkReport.reportEvery = 1
#process.MessageLogger.destinations += ['AnalyzerMessages']
#process.MessageLogger.categories   += ['IMPGP']
#process.MessageLogger.debugModules += ['globalMuons']
#process.MessageLogger.AnalyzerMessages = cms.untracked.PSet(
#    threshold  = cms.untracked.string('DEBUG'),
#    default    = cms.untracked.PSet(limit = cms.untracked.int32(0)),
#    IMPGP = cms.untracked.PSet(limit = cms.untracked.int32(-1))
#    )


process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = 'START38_V13::All'

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        ## Produced with MuonAnalysis/Examples/test/patSkims/patMuons_mc_cfg.py
    #'file:tupleMC.root'
    #'/store/user/tucker/InclusiveMu15/datamc_inclmu15/acbc16136637b3a76f8dd1607df2c067/pat_99_1_CJS.root'
    #'file:copy.root'
    #'/store/user/tucker/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia/datamc_zmumu/b4341788d83565203f0d6250b5475e6e/pat_9_1_etd.root',
    #'/store/user/tucker/TTJets_TuneZ2_7TeV-madgraph-tauola/ttbar_merge/46b5d51805a79be16a6eb0a24eeb03e8/merged_2_1_Ggi.root',
    '/store/user/tucker/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/datamc_inclmu15/6cdfcbb6530d498bbc64d1fe931c3bd6/pat_9_1_uOl.root',
    '/store/user/tucker/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/datamc_inclmu15/6cdfcbb6530d498bbc64d1fe931c3bd6/pat_99_1_6bA.root',
    '/store/user/tucker/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/datamc_inclmu15/6cdfcbb6530d498bbc64d1fe931c3bd6/pat_999_1_1nj.root',
    '/store/user/tucker/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/datamc_inclmu15/6cdfcbb6530d498bbc64d1fe931c3bd6/pat_998_1_V8c.root',
    '/store/user/tucker/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/datamc_inclmu15/6cdfcbb6530d498bbc64d1fe931c3bd6/pat_945_1_d8i.root',
    '/store/user/tucker/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/datamc_inclmu15/6cdfcbb6530d498bbc64d1fe931c3bd6/pat_944_1_OCU.root',
    '/store/user/tucker/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/datamc_inclmu15/6cdfcbb6530d498bbc64d1fe931c3bd6/pat_943_1_COD.root',
    '/store/user/tucker/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/datamc_inclmu15/6cdfcbb6530d498bbc64d1fe931c3bd6/pat_942_1_oAq.root',
    '/store/user/tucker/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/datamc_inclmu15/6cdfcbb6530d498bbc64d1fe931c3bd6/pat_941_1_oML.root',
    '/store/user/tucker/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/datamc_inclmu15/6cdfcbb6530d498bbc64d1fe931c3bd6/pat_940_1_lAh.root',
    '/store/user/tucker/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/datamc_inclmu15/6cdfcbb6530d498bbc64d1fe931c3bd6/pat_93_1_th8.root',
    '/store/user/tucker/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/datamc_inclmu15/6cdfcbb6530d498bbc64d1fe931c3bd6/pat_939_1_9D0.root',
    '/store/user/tucker/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/datamc_inclmu15/6cdfcbb6530d498bbc64d1fe931c3bd6/pat_938_1_xut.root',
    '/store/user/tucker/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/datamc_inclmu15/6cdfcbb6530d498bbc64d1fe931c3bd6/pat_937_1_BA2.root',
    )
                            )

#process.load("sample_2010A")
# qcd 0.12 * 29504866. / 10000. * 10.

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.TFileService = cms.Service('TFileService', fileName=cms.string('plots_simple.root') )

from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
reskim  = hltHighLevel.clone(TriggerResultsTag = cms.InputTag('TriggerResults','',''))
process.recoMu     = reskim.clone(HLTPaths = ['skim_RecoMu'],)
process.bscMinBias = hltHighLevel.clone(HLTPaths = ['HLT_L1_BscMinBiasOR_BptxPlusORMinus'], )

## Common plots
loose_cut = 'isGlobalMuon && ' \
            'innerTrack.pt > 20. && ' \
            'abs(innerTrack.eta) < 2.1 && ' \
            'abs(dB) < 0.2 && ' \
            '(isolationR03.sumPt + isolationR03.emEt + isolationR03.hadEt) / innerTrack.pt < 0.15 && ' \
            'globalTrack.hitPattern.numberOfValidTrackerHits > 10 && ' \
            'globalTrack.hitPattern.numberOfValidPixelHits >= 1 && ' \
            'globalTrack.hitPattern.numberOfValidMuonHits > 0 && ' \
            'numberOfMatches >= 2'

## trigger_match = '(' \
##                 ' triggerMatchesByPath("HLT_Mu9").at(0).pt() > 15) || ' \
##                 '!triggerObjectMatchesByPath("HLT_Mu15_v1").empty()' \
##                 ')'

trigger_match = ' '

tight_cut = loose_cut + trigger_match


from MuonAnalysis.Examples.inclusiveMuonPlots_cfi import makeInclusiveMuonPlots;
process.globalMuons = cms.EDAnalyzer("InclusiveMuonPlots",
    makeInclusiveMuonPlots(0.4,20),
    #muons     = cms.InputTag('patMuonsWithTrigger',''),
    muons     = cms.InputTag('cleanPatMuonsTriggerMatch'),
    selection = cms.string("isGlobalMuon"),
    onlyLeadingMuon = cms.bool(False),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
)
## Now we make the sub-plots with classification by hits
process.globalMuonsGhost = process.globalMuons.clone(selection = "isGlobalMuon && (-5 <= userInt('classByHitsGlb') <= -1)")
process.globalMuonsPunch = process.globalMuons.clone(selection = "isGlobalMuon && ( 0 <= userInt('classByHitsGlb') <=  1)")
process.globalMuonsLight = process.globalMuons.clone(selection = "isGlobalMuon && (userInt('classByHitsGlb') == 2)")
process.globalMuonsHeavy = process.globalMuons.clone(selection = "isGlobalMuon && (userInt('classByHitsGlb') >= 3)")

process.globalMuonsVBTF = process.globalMuons.clone(selection = tight_cut)
process.globalMuonsVBTFGhost = process.globalMuons.clone(selection = tight_cut + " && isGlobalMuon && (-5 <= userInt('classByHitsGlb') <= -1)")
process.globalMuonsVBTFPunch = process.globalMuons.clone(selection = tight_cut + " && isGlobalMuon && ( 0 <= userInt('classByHitsGlb') <=  1)")
process.globalMuonsVBTFLight = process.globalMuons.clone(selection = tight_cut + " && isGlobalMuon && (userInt('classByHitsGlb') == 2)")
process.globalMuonsVBTFHeavy = process.globalMuons.clone(selection = tight_cut + " && isGlobalMuon && (userInt('classByHitsGlb') >= 3)")

from UserCode.Examples.inclusiveMuonPlotsMRTU_cfi import makeInclusiveMuonPlotsMRTU;
process.globalMuonsMRTU = cms.EDAnalyzer("InclusiveMuonPlotsMRTU",
    makeInclusiveMuonPlotsMRTU(),
    #muons     = cms.InputTag('patMuonsWithTrigger',''),
    muons     = cms.InputTag('cleanPatMuonsTriggerMatch'),
    selection = cms.string("isGlobalMuon"),
    #onlyLeadingMuon = cms.bool(False),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
)
process.globalMuonsMRTUVBTF = process.globalMuonsMRTU.clone(selection = tight_cut)

process.p = cms.Path(
#    process.recoMu +
#    process.bscMinBias +
    process.globalMuons +
    process.globalMuonsGhost +    
    process.globalMuonsPunch +    
    process.globalMuonsLight +    
    process.globalMuonsHeavy +
    process.globalMuonsVBTF +
    process.globalMuonsVBTFGhost +    
    process.globalMuonsVBTFPunch +    
    process.globalMuonsVBTFLight +    
    process.globalMuonsVBTFHeavy +
    process.globalMuonsMRTU +
    process.globalMuonsMRTUVBTF    
)

process.demo = cms.EDAnalyzer('CutAnalyzer',
  isMC = cms.untracked.bool(True),
  CrossSection = cms.untracked.double(1.),
  FilterEfficiency = cms.untracked.double(1.),
  TotalNevents = cms.untracked.double(1.),
  IntLuminosity = cms.untracked.double(1.),
  Muon = cms.untracked.InputTag("cleanPatMuonsTriggerMatch"),
  inspectVar = cms.untracked.string(''),
  maxAbsDxyInit = cms.untracked.double(0.2),
  maxAbsEtaInit = cms.untracked.double(2.1),
  minpTInit = cms.untracked.double(20.),
  maxRelCombIsoInit = cms.untracked.double(0.15),#0.15
)

process.pT = process.demo.clone()
process.pT.inspectVar = 'pT'
process.pT.minpTInit = 15. #with step 1 GeV
process.eta = process.demo.clone()
process.eta.inspectVar = 'eta'
process.eta.maxAbsEtaInit = 2.4  # with step 0.1 
process.dxy = process.demo.clone()
process.dxy.inspectVar = 'dxy'
process.dxy.maxAbsDxyInit = 0.25 # with step 0.02 cm
process.relCombIso = process.demo.clone()
process.relCombIso.inspectVar = 'relCombIso'
process.relCombIso.maxRelCombIsoInit = 0.2 # with step 0.01 

process.pCutAnalyzer = cms.Path(process.pT+process.eta+process.dxy+process.relCombIso)
