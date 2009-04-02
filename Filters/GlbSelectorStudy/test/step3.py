# Auto generated configuration file
# using: 
# Revision: 1.99.2.3 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: step2 -s RAW2DIGI,RECO,POSTRECO,ALCA:MuAlCalIsolatedMu+RpcCalHLT,VALIDATION --relval 25000,100 --datatier GEN-SIM-RECO --eventcontent RECOSIM --conditions FrontierConditions_GlobalTag,IDEAL_V9::All --filein file:XXX --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('VAL3')

## process.MessageLogger = cms.Service("MessageLogger",
##     detailedInfo = cms.untracked.PSet(
##         threshold = cms.untracked.string('INFO'),
##         INFO = cms.untracked.PSet(
##             limit = cms.untracked.int32(0)
##         ),
##         extension = cms.untracked.string('.txt')
##     ),
##     debug = cms.untracked.PSet(
##         threshold = cms.untracked.string('DEBUG'),
##         INFO = cms.untracked.PSet(
##             limit = cms.untracked.int32(0)
##         ),
##         SimLevelParent = cms.untracked.PSet(
##             limit = cms.untracked.int32(1000000)
##         ),
##         DEBUG = cms.untracked.PSet(
##             limit = cms.untracked.int32(0)
##         ),
##         MotherSearch = cms.untracked.PSet(
##             limit = cms.untracked.int32(1000000)
##         ),
##         GlbSelectorStudy = cms.untracked.PSet(
##             limit = cms.untracked.int32(1000000)
##         )
##     ),
##     debugModules  = cms.untracked.vstring('glbSelStudy',),
##     categories = cms.untracked.vstring('SimHitsAnlzrImproved', 
##         'GlbSelectorStudy', 
##         'MotherSearch', 
##         'SimLevelParent'),
##     destinations = cms.untracked.vstring('detailedInfo3')
## )

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.categories = ['MotherSearch', 'GlbSelectorStudy']
process.MessageLogger.debugModules = ['glbSelStudy']
process.MessageLogger.cout = cms.untracked.PSet(
    threshold = cms.untracked.string('DEBUG'),
    default = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
    ),
    MotherSearch = cms.untracked.PSet(
        limit = cms.untracked.int32(-1)
    ),
    GlbSelectorStudy = cms.untracked.PSet(
        limit = cms.untracked.int32(-1)
    ),
)
process.MessageLogger.cerr = cms.untracked.PSet(
    placeholder = cms.untracked.bool(True)
)



# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
#process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/MixingNoPileUp_cff')
process.load('Configuration/StandardSequences/GeometryPilot2_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
#process.load('Configuration/StandardSequences/RawToDigi_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/StandardSequences/PostRecoGenerator_cff')
#process.load('Configuration/StandardSequences/AlCaRecoStreams_cff')
process.load('Configuration/StandardSequences/Validation_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    annotation = cms.untracked.string('step2 nevts:1'),
    name = cms.untracked.string('PyReleaseValidation')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)
process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)
# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:step2_RAW2DIGI_RECO_POSTRECO_ALCA_VALIDATION.root')
)

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    outputCommands = process.RECOSIMEventContent.outputCommands,
    fileName = cms.untracked.string('step3.root'),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    )
)


# Additional output definition
process.output.outputCommands.append("keep *")
# Other statements
process.TFileService = cms.Service("TFileService", fileName = cms.string("histo3.root"))

process.GlobalTag.globaltag = 'IDEAL_V11::All'
process.validation = cms.Sequence(process.mix+process.globaldigisanalyze*process.globalhitsanalyze*process.globalrechitsanalyze*process.globalValidation)

process.load("UserCode.GlbSelectorStudy.glbselectorstudy_cfi")
process.glbSelStudy.doAssoc = False
process.glbSelStudy.trkMuAssocLabel = "tpToTkmuTrackAssociation"
process.glbSelStudy.staMuAssocLabel = "tpToStaTrackAssociation"
process.glbSelStudy.glbMuAssocLabel = "tpToGlbTrackAssociation"
process.p = cms.Path(process.glbSelStudy)


# Path and EndPath definitions
#process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path(process.reconstruction)
process.postreco_step = cms.Path(process.postreco_generator)

process.validation_step = cms.Path(process.validation)
process.endjob_step = cms.Path(process.endOfProcess)
process.out_step = cms.EndPath(process.output)

# Schedule definition
#process.schedule = cms.Schedule(process.raw2digi_step,process.reconstruction_step,process.postreco_step,process.pathALCARECORpcCalHLT,process.pathALCARECOMuAlCalIsolatedMu,process.validation_step,process.p,process.endjob_step,process.out_step)
process.schedule = cms.Schedule(process.p)
