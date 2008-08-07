from RecoZ.MuonAnalyser.test.replaces.mu1000_cff import *
from Configuration.StandardSequences.Geometry_cff import *
from Configuration.StandardSequences.MagneticField_cff import *
from Configuration.StandardSequences.FrontierConditions_GlobalTag_cff import *
from RecoZ.MuonAnalyser.MuonAnalyser_cff import *
from RecoZ.MuonAnalyser.MuonQualityAnalyser_cff import *
pMuonZ = cms.Path(MuonAnalyser)
GlobalTag.globaltag = 'IDEAL_V2::All'
MuonAnalyser.SimAssociator.SimHitsOccurrences = 0
MuonAnalyser.OutputFileName = '$outFileName_ric.root'
MuonQualityAnalyser.OutputFileName = '$outFileName_mqa.root'
