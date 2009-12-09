# Auto generated configuration file
# using: 
# Revision: 1.138.2.1 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: run -s L1,HLT:1E31 -n 100 --conditions FrontierConditions_GlobalTag,IDEAL_31X::All --mc --filein /store/relval/CMSSW_3_1_0_pre6/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/IDEAL_31X_v1/0002/926A6C4C-FB32-DE11-A60D-001617C3B6CC.root --fileout refHLT.root --eventcontent FEVTDEBUGHLT --processName OfflineHLT --customise UserCode/L3Switches/customise.py --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('OfflineHLT')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/MixingNoPileUp_cff')
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/SimL1Emulator_cff')
process.load('HLTrigger/Configuration/HLT_1E31_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.138.2.1 $'),
    annotation = cms.untracked.string('run nevts:100'),
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
    '/store/relval/CMSSW_3_2_6/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/MC_31X_V8-v1/0013/5243BC31-7A9A-DE11-80E0-0019B9F581C9.root'
    )
                            )

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    fileName = cms.untracked.string('refHLT.root'),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        filterName = cms.untracked.string('')
    )
)

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'MC_31X_V8::All'

# Path and EndPath definitions
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.endjob_step = cms.Path(process.endOfProcess)
process.out_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(process.L1simulation_step)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.endjob_step,process.out_step])

# Automatic addition of the customisation function

def customise(process):
    #process.load(\'sample\')
    from Workspace.MuonHLTTreeUtility.muonHLTTreeUtility_cff import insertMHTU
    insertMHTU(process)
    process.hltMuonTreeMaker.outputFileName = "Tree_Iter3_OISOIH_muon.root"
    process.TFileService.fileName = "Bit_Iter3_OISOIH_muon.root"
    
    import UserCode.L3Switches.Switches as switch
    #switch.SwitchToBaseline(process)
    #switch.SwitchToBaselinePP(process)
    #switch.SwitchToOIState(process)
    #switch.SwitchToOIHit(process)
    #switch.SwitchToAllCombined(process)
    #switch.SwitchToOICombined(process)
    #switch.SwitchToIterative(process)
    switch.SwitchToIterative3(process)
     
    return(process)


# End of customisation function definition

process = customise(process)
