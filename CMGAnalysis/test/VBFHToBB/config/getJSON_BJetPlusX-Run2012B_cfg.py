import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.load('FWCore.MessageService.MessageLogger_cfi')

from CMGTools.Production.datasetToSource import *
from CMGTools.Common.Tools.applyJSON_cff import *

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(-1)
        )

process.source = datasetToSource('cmgtools','/BJetPlusX/Run2012B-13Jul2012-v1/AOD/PAT_CMG_V5_12_0','cmgTuple_.*.root')

applyJSON(process,'Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON_v2.txt')

#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
##-------------------- User analyzer  --------------------------------
process.json = cms.EDAnalyzer('JSONProducer',
  filename = cms.string("json_BJetPlusX-Run2012B.txt")
)

process.p = cms.Path(process.json)
