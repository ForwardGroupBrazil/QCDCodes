from RecoZ.MuonAnalyser.test.replaces.$cacheName_cff import *

process.load("Configuration.StandardSequences.Geometry_cff")

process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.load("RecoZ.MuonAnalyser.MuonAnalyser_cff")

process.load("RecoZ.MuonAnalyser.MuonQualityAnalyser_cff")

process.pMuonZ = cms.Path(process.MuonAnalyser)

process.GlobalTag.globaltag = 'IDEAL_V2::All'
process.MuonAnalyser.SimAssociator.SimHitsOccurrences = 0
process.MuonAnalyser.OutputFileName = '$outFileName_ric.root'
process.MuonQualityAnalyser.OutputFileName = '$outFileName_mqa.root'

#replaces
