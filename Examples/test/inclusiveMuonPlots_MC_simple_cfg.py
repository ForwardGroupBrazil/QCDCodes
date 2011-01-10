import FWCore.ParameterSet.Config as cms

process = cms.Process("MuonPlots")

# Messages
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

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
    'file:tupleMC.root'
    #'/store/user/tucker/InclusiveMu15/datamc_inclmu15/acbc16136637b3a76f8dd1607df2c067/pat_99_1_CJS.root'
    #'file:copy.root'
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.TFileService = cms.Service('TFileService', fileName=cms.string('inclusiveMuonPlots_MC_simple.class.root') )

from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
reskim  = hltHighLevel.clone(TriggerResultsTag = cms.InputTag('TriggerResults','',''))
process.recoMu     = reskim.clone(HLTPaths = ['skim_RecoMu'],)
process.bscMinBias = hltHighLevel.clone(HLTPaths = ['HLT_L1_BscMinBiasOR_BptxPlusORMinus'], )

## Common plots
from MuonAnalysis.Examples.inclusiveMuonPlots_cfi import makeInclusiveMuonPlots;
process.globalMuons = cms.EDAnalyzer("InclusiveMuonPlots",
    makeInclusiveMuonPlots(),
    muons     = cms.InputTag('patMuonsWithTrigger',''),
    #muons     = cms.InputTag('cleanPatMuonsTriggerMatch'),
    selection = cms.string("isGlobalMuon"),
    onlyLeadingMuon = cms.bool(False),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
)
## Now we make the sub-plots with classification by hits
process.globalMuonsGhost = process.globalMuons.clone(selection = "isGlobalMuon && (-5 <= userInt('classByHitsGlb') <= -1)")
process.globalMuonsPunch = process.globalMuons.clone(selection = "isGlobalMuon && ( 0 <= userInt('classByHitsGlb') <=  1)")
process.globalMuonsLight = process.globalMuons.clone(selection = "isGlobalMuon && (userInt('classByHitsGlb') == 2)")
process.globalMuonsHeavy = process.globalMuons.clone(selection = "isGlobalMuon && (userInt('classByHitsGlb') >= 3)")

from UserCode.Examples.inclusiveMuonPlotsMRTU_cfi import makeInclusiveMuonPlotsMRTU;
process.globalMuonsMRTU = cms.EDAnalyzer("InclusiveMuonPlotsMRTU",
    makeInclusiveMuonPlotsMRTU(),
    muons     = cms.InputTag('patMuonsWithTrigger',''),
    #muons     = cms.InputTag('cleanPatMuonsTriggerMatch'),
    selection = cms.string("isGlobalMuon"),
    #onlyLeadingMuon = cms.bool(False),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
)

process.p = cms.Path(
#    process.recoMu +
#    process.bscMinBias +
    process.globalMuons +
    process.globalMuonsMRTU +    
    process.globalMuonsGhost +    
    process.globalMuonsPunch +    
    process.globalMuonsLight +    
    process.globalMuonsHeavy     
)
