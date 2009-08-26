# Auto generated configuration file
# using: 
# Revision: 1.123 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: MRTU -s L1,HLT:1E31 -n 100 --processName MHTU --conditions FrontierConditions_GlobalTag,MC_31X_V1::All --filein HLTDEBUG-file.root --no_output --customise Workspace/MuonRecoTreeUtility/muonRecoTreeUtility_customise.py --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('MRTU')

process.MessageLogger = cms.Service("MessageLogger",
    detailedInfo = cms.untracked.PSet(
        threshold = cms.untracked.string('INFO'),
        INFO = cms.untracked.PSet(
            limit = cms.untracked.int32(-1)
        ),
        extension = cms.untracked.string('.txt')
    ),
    debug = cms.untracked.PSet(
        threshold = cms.untracked.string('DEBUG'),
        extension = cms.untracked.string('.txt'),
        INFO = cms.untracked.PSet(
            limit = cms.untracked.int32(-1)
        ),
        DEBUG = cms.untracked.PSet(
            limit = cms.untracked.int32(100000)
        ),
        MuonRecoTreeUtility = cms.untracked.PSet(
            limit = cms.untracked.int32(100000)
        )        
    ),
    categories = cms.untracked.vstring(
     'MuonRecoTreeUtility',
     ),
    destinations = cms.untracked.vstring('detailedInfo','debug'),
    debugModules  = cms.untracked.vstring('hltMuonTreeMaker'),

)


# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
#process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/MixingNoPileUp_cff')
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/SimL1Emulator_cff')
process.load('Configuration/StandardSequences/RawToDigi_cff')
process.load('Configuration/StandardSequences/PostRecoGenerator_cff')
process.load('Configuration/StandardSequences/Validation_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.123 $'),
    annotation = cms.untracked.string('MHTU nevts:100'),
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
    fileNames = cms.untracked.vstring(
     '/store/user/aeverett/CMSSW_3_2_1/SingleMuPt200//aeverett//SingleMuPt200_CMSSW_3_2_1_step1//SingleMuPt200_CMSSW_3_2_1_step1//3f9251f9631b9dac27d828bfadc5f286//step1_29.root',
    )
)

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    fileName = cms.untracked.string('MRTU_Reco.root'),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        filterName = cms.untracked.string('')
    )
)

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'MC_31X_V1::All'
process.ttrhbwor.ComputeCoarseLocalPositionFromDisk = True
process.ttrhbwr.ComputeCoarseLocalPositionFromDisk = True

# Path and EndPath definitions
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.endjob_step = cms.Path(process.endOfProcess)
process.out_step = cms.EndPath(process.output)




# Schedule definition
#process.schedule = cms.Schedule()
process.schedule = cms.Schedule(process.L1simulation_step)
#bbb process.schedule.extend(process.HLTriggerFirstPath)
#bbb process.schedule.extend([process.HLTriggerFinalPath,process.HLTAnalyzerEndpath])
#process.schedule.extend(process.HLTSchedule)
##
process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path(process.reconstruction)
##
process.schedule.extend([process.raw2digi_step,process.reconstruction_step])
process.schedule.extend([process.endjob_step,process.out_step])

# Automatic addition of the customisation function
def customise(process):
     from Workspace.MuonRecoTreeUtility.muonRecoTreeUtility_cff import insertMHTU
     insertMHTU(process)
     return (process)

## def customise(process):
##     process.load('Workspace.MuonHLTTreeUtility.muonHLTTreeUtility_cff')
##     from Workspace.MuonHLTTreeUtility.muonHLTTreeUtility_cff import muonHLTrecoSequence
##     process.muonHLTreco = muonHLTrecoSequence(process)
##     process.MHTU_Path+=process.muonHLTreco
##     import FWCore.ParameterSet.SequenceTypes
##     for p in process.schedule:
##         if (p.__class__==FWCore.ParameterSet.SequenceTypes.EndPath):
##             process.schedule.insert(process.schedule.index(p), process.MHTU_Path )
##             break
##     process.schedule.append( process.MHTU_EndPath )
##     ##actually do the --no_output option
##     if (hasattr(process,"out_step")):
##         process.schedule.remove(process.out_step)
##     return (process) 
# End of customisation function definition

process = customise(process)
