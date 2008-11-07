import FWCore.ParameterSet.Config as cms

# Bit Plotting
bitSummary = cms.EDAnalyzer(
    "BitPlotting",
     out = cms.untracked.string('dqm.root'),
     HLTPaths = cms.vstring('HLT_L1MuOpen','HLT_L1Mu','HLT_L2Mu9',
                           'HLT_IsoMu15','HLT_IsoMu11','HLT_IsoMu9',
                           'HLT_Mu15','HLT_Mu13','HLT_Mu11','HLT_Mu9',
                           'HLT_Mu7','HLT_Mu5','HLT_Mu3','HLT_DoubleMu3',
                           'HLT_DoubleIsoMu3'),
     TriggerResultsTag = cms.InputTag('TriggerResults','','HLT2')
    )

