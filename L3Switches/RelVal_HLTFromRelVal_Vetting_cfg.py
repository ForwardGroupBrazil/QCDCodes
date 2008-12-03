import FWCore.ParameterSet.Config as cms

process = cms.Process("HLT2")

process.load("Configuration.StandardSequences.Services_cff")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
#process.load("RecoMuon.Configuration.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#    '/store/user/hyxu/TTbar-newStep3/newStep3-TTbar-0099.root',
'/store/user/aeverett/note2112/SingleMuPt500/SingleMuPt500_cfi_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_7.root',
    )
                            )

process.load("Configuration.StandardSequences.GeometryPilot2_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

# Conditions: fake or frontier
# process.load("Configuration.StandardSequences.FakeConditions_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'IDEAL_V9::All'

process.load("Configuration.StandardSequences.L1Emulator_cff")
# Choose a menu/prescale/mask from one of the choices
# in L1TriggerConfig.L1GtConfigProducers.Luminosity
process.load("Configuration.StandardSequences.L1TriggerDefaultMenu_cff")


# run HLT
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("HLTrigger.Configuration.HLT_2E30_cff")
#process.schedule = process.HLTSchedule

process.load("UserCode.L3Switches.SwitchToIOHit_cff")
#process.load("UserCode.L3Switches.SwitchToOIHit_cff")
#process.load("UserCode.L3Switches.SwitchToOIState_cff")


# If you want to remake the TrackingParticle collection
# redo the DigiLinks too
process.load("UserCode.L3Switches.TPLink_cfi")
process.simParticle_step = cms.Path(process.TPLink)


MuonHLTSchedule = cms.Schedule(
    process.simParticle_step,
    process.HLT_DoubleMu3,
    process.HLT_DoubleIsoMu3,
    process.HLT_Mu3, process.HLT_Mu5,
    process.HLT_Mu7, process.HLT_Mu9,
    process.HLT_Mu11, process.HLT_Mu13,
    process.HLT_Mu15,
    process.HLT_IsoMu9, process.HLT_IsoMu11,
    process.HLT_IsoMu15,
    process.HLT_L2Mu9,
    process.HLT_L1MuOpen, process.HLT_L1Mu,
    )
process.schedule = cms.Schedule()
process.schedule.extend( MuonHLTSchedule )

process.load("Validation.RecoMuon.muonValidationHLT_cff")
process.muonValidationHLT_step = cms.Path(process.recoMuonValidationHLT_seq)
process.schedule.append(process.muonValidationHLT_step)

process.hltL1gtTrigReport = cms.EDAnalyzer( "L1GtTrigReport",
    UseL1GlobalTriggerRecord = cms.bool( False ),
    L1GtRecordInputTag = cms.InputTag( "hltGtDigis" )
)
process.hltTrigReport = cms.EDAnalyzer( "HLTrigReport",
    HLTriggerResults = cms.InputTag( 'TriggerResults','','HLT2' )
)
process.HLTAnalyzerEndpath = cms.EndPath( process.hltL1gtTrigReport + process.hltTrigReport )
process.schedule.append(process.HLTAnalyzerEndpath)

process.load("Configuration.EventContent.EventContent_cff")
process.hltPoolOutput = cms.OutputModule("PoolOutputModule",
    process.FEVTDEBUGHLTEventContent,
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('RECO')
    ),
    basketSize = cms.untracked.int32(4096),
    fileName = cms.untracked.string('file:HLTFromDigiRaw.root')
)

process.load("DQMServices.Components.MEtoEDMConverter_cfi")
process.MEtoEDMConverter_step = cms.Path(process.MEtoEDMConverter)

process.schedule.append(process.MEtoEDMConverter_step)

# Bit Plotting
process.bitSummary = cms.EDAnalyzer(
    "BitPlotting",
    out = cms.untracked.string('file:bit.root'),
    HLTPaths = cms.vstring('HLT_L1MuOpen','HLT_L1Mu','HLT_L2Mu9',
                           'HLT_IsoMu15','HLT_IsoMu11','HLT_IsoMu9',
                           'HLT_Mu15','HLT_Mu13','HLT_Mu11','HLT_Mu9',
                           'HLT_Mu7','HLT_Mu5','HLT_Mu3','HLT_DoubleMu3',
                           'HLT_DoubleIsoMu3'),
    TriggerResultsTag = cms.InputTag('TriggerResults','','HLT2')
    )
process.BitSummaryEndPath = cms.EndPath( process.bitSummary )

process.schedule.append( process.BitSummaryEndPath )

# To include timing (via hltTimingSummary)
process.PathTimerService = cms.Service( "PathTimerService" )
process.timer = cms.EDProducer( "PathTimerInserter" )
#process.hltPoolOutput.outputCommands.append('drop *')
process.hltPoolOutput.outputCommands.append('keep HLTPerformanceInfo_*_*_*')
process.endp1 = cms.EndPath( process.timer + process.hltPoolOutput)

process.schedule.append( process.endp1  )

#process.muonCkfTrajectoryFilter.filterPset.maxNumberOfHits = 6
process.mergedtruth.vertexDistanceCut = 1000
process.muonTPSet.tip = 10000
process.muonTPSet.lip = 10000
process.l2MuonTrackV.tipTP = 10000
process.l2MuonTrackV.lipTP = 10000
process.l3MuonTrackV.tipTP = 10000
process.l3MuonTrackV.lipTP = 10000
process.l3TkMuonTrackV.tipTP = 10000
process.l3TkMuonTrackV.lipTP = 10000
process.l3TkMuonMuTrackV.tipTP = 10000
process.l3TkMuonMuTrackV.lipTP = 10000
process.l2MuonMuTrackV.tipTP = 10000
process.l2MuonMuTrackV.lipTP = 10000
process.l3MuonMuTrackV.tipTP = 10000
process.l3MuonMuTrackV.lipTP = 10000

process.TrackAssociatorByPosDeltaR.method = 'posdr'
process.TrackAssociatorByPosDeltaR.QCut = 0.1

import SimTracker.TrackAssociation.TrackAssociatorByPosition_cfi
process.TrackAssociatorByPosDeltaR2 = SimTracker.TrackAssociation.TrackAssociatorByPosition_cfi.TrackAssociatorByPosition.clone()
process.TrackAssociatorByPosDeltaR2.method = 'posdr'
process.TrackAssociatorByPosDeltaR2.QCut = 0.2
process.TrackAssociatorByPosDeltaR2.ComponentName = 'TrackAssociatorByDeltaR2'

process.tpToL2TrackAssociation.associator = 'TrackAssociatorByDeltaR2'
process.l2MuonTrackV.associators = 'TrackAssociatorByDeltaR2'
