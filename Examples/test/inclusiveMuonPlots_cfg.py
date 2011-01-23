import FWCore.ParameterSet.Config as cms

process = cms.Process("MuonPlots")

# Messages
process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.MessageLogger.destinations += ['AnalyzerMessages']
process.MessageLogger.categories   += ['Dumper']
process.MessageLogger.debugModules += ['patTupleDumper','globalMuons']
process.MessageLogger.AnalyzerMessages = cms.untracked.PSet(
    threshold  = cms.untracked.string('DEBUG'),
    default    = cms.untracked.PSet(limit = cms.untracked.int32(0)),
    Dumper = cms.untracked.PSet(limit = cms.untracked.int32(-1))
    )


process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.GlobalTag.globaltag = 'START3X_V25B::All'

process.source = cms.Source(
    "PoolSource",
    noEventSort = cms.untracked.bool(True),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    fileNames = cms.untracked.vstring(
    '/store/user/tucker/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia/datamc_zmumu/b4341788d83565203f0d6250b5475e6e/pat_9_1_etd.root',
    #'/store/user/tucker/TTJets_TuneZ2_7TeV-madgraph-tauola/ttbar_merge/46b5d51805a79be16a6eb0a24eeb03e8/merged_2_1_Ggi.root',
    #'/store/user/aeverett//Mu_Pt0_500//pat_skim_101_1.root',
    #'/store/user/aeverett//Pi_Pt2_200//pat_skim_1_1.root','/store/user/aeverett//Pi_Pt2_200//pat_skim_54_1.root','/store/user/aeverett//Pi_Pt2_200//pat_skim_84_1.root','/store/user/aeverett//Pi_Pt2_200//pat_skim_6_1.root',
    )
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000) )
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True),
                                        fileMode    = cms.untracked.string("NOMERGE") )

## Define some utilities to declare bins easily
def nBins(n,min,max): return cms.vdouble(*[min + (max-min)/n*i for i in range(n+1)])
def evenBins(min,max,delta): 
    ret = cms.vdouble(min)
    x = min
    while x < max - 1e-4: # need a small hint otherwise for some numbers it will overstep due to numerical resolution
        x += delta
        ret.append(x)
    return ret 

process.TFileService = cms.Service('TFileService',
    fileName=cms.string('inclusiveMuonPlots_Data.test2.root')
)

###
# InclusiveMuonPlots can be run on PAT or RECO
# MuonAnalysis has general kinematic plots
# UserCode adds more figures about hits, segments, etc.
#from MuonAnalysis.Examples.inclusiveMuonPlots_cfi import makeInclusiveMuonPlots;
from UserCode.Examples.inclusiveMuonPlotsMRTU_cfi import makeInclusiveMuonPlotsMRTU;
commonInputs = cms.PSet(
    muons     = cms.InputTag('patMuons'),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
)
process.trackerMuons = cms.EDAnalyzer("InclusiveMuonPlotsMRTU",
    makeInclusiveMuonPlotsMRTU(),
    commonInputs,
    selection = cms.string("isTrackerMuon && muonID('TMLastStationAngTight')"),
)
process.globalMuons = process.trackerMuons.clone(
    selection = "isGlobalMuon"
    )
process.standAloneMuons = process.trackerMuons.clone(
    selection = "isStandAloneMuon"
    )
process.exclusiveGlobalMuons = process.trackerMuons.clone(
    selection = "isGlobalMuon && !isTrackerMuon"
    )
process.nonGlobalMuons = process.trackerMuons.clone(
    selection = "isTrackerMuon && isStandAloneMuon && !isGlobalMuon"
    )

import UserCode.Examples.patTupleDumper_cfi
process.patTupleDumper = UserCode.Examples.patTupleDumper_cfi.patTupleDumper.clone(muons = 'patMuons',selection = 'isGlobalMuon || isTrackerMuon')
process.patTupleDumperPunch = UserCode.Examples.patTupleDumper_cfi.patTupleDumper.clone(muons = 'patMuons', selection = '(isGlobalMuon && abs(userInt(\'classByHitsGlb\')) <= 1) || (isTrackerMuon && abs(userInt(\'classByHitsTM\')) <= 1)',outputFileName='ntuplePunch.root')

process.p = cms.Path(process.trackerMuons *
                     process.globalMuons *
                     process.standAloneMuons *
                     process.exclusiveGlobalMuons *
                     process.nonGlobalMuons
                     * process.patTupleDumper
                     * process.patTupleDumperPunch
                     )



if True: ## use when starting from pat::Muons
    process.trackerMuons.muons = 'patMuons'
    process.globalMuons.muons  = 'patMuons'


