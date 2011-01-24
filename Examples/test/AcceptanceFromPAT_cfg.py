#!/usr/bin/env python

import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process

process.GlobalTag.globaltag = 'START38_V13::All'

#process.source.fileNames = ['/store/user/aeverett/ExoMu101102bis//ZP750//aeverett//ZprimeSSMToMuMu_M-750_7TeV-pythia6//ZP750_PATbis//3d1d7ee1efaf83bebd1a3c5daf37c48f//pat_3_1_WtD.root',]
process.source.fileNames = [$inputFileNames]
#process.source.fileNames = ['file:/home/ba01/u112/aeverett/scratch_rcac/110111Figs/pat_skim_Zmass__lowerMasses_Zmumu.root']
#process.source.fileNames = ['file:/home/ba01/u112/aeverett/scratch_rcac/110111Figs/pat_skim_Zmass_Zpeak_Zmumu.root']
#process.source.fileNames = ['file:/home/ba01/u112/aeverett/scratch_rcac/110111Figs/pat_skim_diMass_200up_Zmumu.root']
#process.source.eventsToProcess = [ '1:760580','1:760611','1:526898','1:158876', ]
input_is_MC = True

## process.MessageLogger.cerr.FwkReport.reportEvery = 1
## process.MessageLogger.destinations += ['AnalyzerMessages',]

## process.MessageLogger.categories   += ['ZP2M',"ZPMRTU"]
## process.MessageLogger.debugModules += ['genMuonsStoyan',]

## process.MessageLogger.AnalyzerMessages = cms.untracked.PSet(
##    threshold  = cms.untracked.string('DEBUG'),
##    default    = cms.untracked.PSet(limit = cms.untracked.int32(0)),
##    ZP2M = cms.untracked.PSet(limit = cms.untracked.int32(-1)),
##    ZPMRTU = cms.untracked.PSet(limit = cms.untracked.int32(-1))
##    )


process.load('SUSYBSMAnalysis.Zprime2muAnalysis.HistosFromPAT_cfi')
process.HistosFromPATExtract = process.HistosFromPAT.clone(leptonsFromDileptons = True)
process.p = cms.Path(process.goodDataFilter * process.Zprime2muAnalysisSequence * process.HistosFromPAT * process.HistosFromPATExtract)


##########################

def nBins(n,min,max): return cms.vdouble(*[min + (max-min)/n*i for i in range(n+1)])
def evenBins(min,max,delta): 
    ret = cms.vdouble(min)
    x = min
    while x < max - 1e-4: # need a small hint otherwise for some numbers it will overstep due to numerical resolution
        x += delta
        ret.append(x)
    return ret 

###########################

from SUSYBSMAnalysis.Zprime2muAnalysis.inclusiveDiMuonPlots_cfi import makeInclusiveDiMuonPlots
commonInputs = cms.PSet(
    dilepton_src = cms.InputTag('dimuons'),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
    selection = cms.string(""),
)
process.dileptonPlots = cms.EDAnalyzer("Zprime2muAnalysisPlots",
                                 makeInclusiveDiMuonPlots(),
                                 commonInputs)

process.MuonsFromDimuons = cms.EDProducer('MuonsFromDimuons', dimuon_src = cms.InputTag('dimuons'))

####
##
## Common plots
process.load("Configuration.StandardSequences.Reconstruction_cff")
from MuonAnalysis.Examples.inclusiveMuonPlots_cfi import makeInclusiveMuonPlots as makeInclusiveMuonPlotsGeneral;
process.globalMuons = cms.EDAnalyzer("InclusiveMuonPlots",
    makeInclusiveMuonPlotsGeneral(),
    muons     = cms.InputTag('cleanPatMuonsTriggerMatch'),
    selection = cms.string("isGlobalMuon"),
    onlyLeadingMuon = cms.bool(False),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
)
##
####

from UserCode.Examples.inclusiveMuonPlotsMRTU_cfi import makeInclusiveMuonPlotsMRTU;
commonInputsLepton = cms.PSet(
    muons     = cms.InputTag('cleanPatMuonsTriggerMatch'),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
)
process.usedMuonsMRTU = cms.EDAnalyzer("InclusiveMuonPlotsMRTU",
    makeInclusiveMuonPlotsMRTU(),
    commonInputsLepton,
    selection = cms.string(""),
)
process.usedMuonsMRTU.muons = 'MuonsFromDimuons'

process.gMuonsMRTU = cms.EDAnalyzer("InclusiveMuonPlotsMRTU",
    makeInclusiveMuonPlotsMRTU(),
    commonInputsLepton,
    selection = cms.string("isGlobalMuon")
    )
process.tMuonsMRTU = process.gMuonsMRTU.clone(
    selection = "isTrackerMuon"
    )

process.p3 = cms.Path(
    process.goodDataFilter *
    process.Zprime2muAnalysisSequence *
    process.MuonsFromDimuons *
    #adam process.dileptonPlots *
    process.usedMuonsMRTU *
    process.gMuonsMRTU *
    process.tMuonsMRTU *
    process.globalMuons
    )

######################

if input_is_MC:
    
    vbtfselection = cms.string(
        'isGlobalMuon && isTrackerMuon && '+
        'track.hitPattern.numberOfValidPixelHits > 0 && ' +
        'track.hitPattern.numberOfValidTrackerHits > 10 && '+
        'numberOfMatches > 1 && '+
        'abs(dB) < 0.2 && '+
        'globalTrack.normalizedChi2 < 10 && '+
        'globalTrack.hitPattern.numberOfValidMuonHits > 0 && '+
        'pt > 7. && abs(eta) < 2.1 && ' +
        '(isolationR03.sumPt + isolationR03.emEt + isolationR03.hadEt) / innerTrack.pt < 0.15 '
        )
    process.dyDimuons = cms.EDFilter("PATMuonRefSelector",
                                  src = cms.InputTag('cleanPatMuonsTriggerMatch'),
                                  cut = vbtfselection, 
                                  )
    process.dyz = cms.EDProducer("CandViewShallowCloneCombiner",
                                       decay = cms.string('dyDimuons@+ dyDimuons@-'),
                                       cut = cms.string('0.0 < mass < 20000.0'),
                                       name = cms.string('dyz'),
                                       roles = cms.vstring('muon1', 'muon2')
                                       )

    
    
    from UserCode.Examples.inclusiveMuonPlotsGENSIM_cfi import makeInclusiveMuonPlots;
    commonInputs = cms.PSet(
        muons     = cms.InputTag('cleanPatMuonsTriggerMatch'),
        particleSrc = cms.InputTag('prunedGenSimLeptons'), #genParticles
        dilepton_src = cms.InputTag('dimuons'), #('dimuons'), #('dyz'),
        primaryVertices = cms.InputTag("offlinePrimaryVertices"),
        weight = cms.untracked.double(1.0),
        eta_acc1 = cms.untracked.double(2.1),
        eta_acc2 = cms.untracked.double(2.4),
        pt_acc1 = cms.untracked.double(16.0),
        pt_acc2 = cms.untracked.double(7.0),
        mother = cms.untracked.int32(23),#Z: 23, ZPrime: 32
        daughter = cms.untracked.int32(13),
        daughterStatus = cms.untracked.int32(1),
        )

    process.allGenSim = cms.EDAnalyzer(
        "InclusiveMuonPlotsGENSIM",
        makeInclusiveMuonPlots(),
        commonInputs,
        selection = cms.string(""),
        selectionReco = cms.string(""),
        ptetaYBins = cms.uint32(200),
        ptetaYRange = cms.vdouble(0.,200.),
        ptetaXBins = cms.uint32(260),
        ptetaXRange = cms.vdouble(-2.6,2.6),
        petaYBins = cms.uint32(200),
        petaYRange = cms.vdouble(0.,200.),
        petaXBins = cms.uint32(260),
        petaXRange = cms.vdouble(-2.6,2.6),
        #massBins = nBins(100,0,200.),
        massBins = cms.vdouble(15, 20, 30, 40, 50, 60, 76, 86, 96, 106, 120, 150, 200, 600, 2000),
        )
    process.allGenSim.ptBins = evenBins(0,600.,5.)
    process.allGenSim.pBins = evenBins(0,600.,5.)

#    process.allGenSim3 = process.allGenSim.clone(
#        selection = "status == 3"
#        )
    process.genMuons = process.allGenSim.clone(
        selection = " mass > 12 && mass < 600",
        selectionReco = "isGlobalMuon && isTrackerMuon"
        )
    process.genMuonsStoyan = process.allGenSim.clone(
        #daughterStatus = 3,
        selection = " mass > 12 && mass < 600",
        selectionReco = "isGlobalMuon && isTrackerMuon",
        massBins = cms.vdouble(15,600),
        )
    process.genMuonsBin1 = process.allGenSim.clone(
        selection = "mass > 15 && mass < 20",
        selectionReco = "isGlobalMuon && isTrackerMuon"
        )
    process.genMuonsBin2 = process.allGenSim.clone(
        selection = "mass > 20 && mass < 30",
        selectionReco = "isGlobalMuon && isTrackerMuon"
        )
    process.genMuonsBin3 = process.allGenSim.clone(
        selection = "mass > 30 && mass < 40",
        selectionReco = "isGlobalMuon && isTrackerMuon"
        )
    process.genMuonsBin4 = process.allGenSim.clone(
        selection = "mass > 40  && mass < 50",
        selectionReco = "isGlobalMuon && isTrackerMuon"
        )
    process.genMuonsBin5 = process.allGenSim.clone(
        selection = "mass > 50 && mass < 60",
        selectionReco = "isGlobalMuon && isTrackerMuon"
        )
    process.genMuonsBin6 = process.allGenSim.clone(
        selection = "mass > 60 && mass < 76",
        selectionReco = "isGlobalMuon && isTrackerMuon"
        )
    process.genMuonsBin7 = process.allGenSim.clone(
        selection = "mass > 76 && mass < 86",
        selectionReco = "isGlobalMuon && isTrackerMuon"
        )
    process.genMuonsBin8 = process.allGenSim.clone(
        selection = "mass > 86 && mass < 96",
        selectionReco = "isGlobalMuon && isTrackerMuon"
        )
    process.genMuonsBin9 = process.allGenSim.clone(
        selection = "mass > 96 && mass < 106",
        selectionReco = "isGlobalMuon && isTrackerMuon"
        )
    process.genMuonsBin10 = process.allGenSim.clone(
        selection = "mass > 106 && mass < 120",
        selectionReco = "isGlobalMuon && isTrackerMuon"
        )
    process.genMuonsBin11 = process.allGenSim.clone(
        selection = "mass > 120 && mass < 150",
        selectionReco = "isGlobalMuon && isTrackerMuon"
        )
    process.genMuonsBin12 = process.allGenSim.clone(
        selection = "mass > 150 && mass < 200",
        selectionReco = "isGlobalMuon && isTrackerMuon"
        )
    process.genMuonsBin13 = process.allGenSim.clone(
        selection = "mass > 200 && mass < 600",
        selectionReco = "isGlobalMuon && isTrackerMuon"
        )
    process.genMuonsBin14 = process.allGenSim.clone(
        selection = "mass > 600 && mass < 2000",
        selectionReco = "isGlobalMuon && isTrackerMuon"
        )

    process.p2 = cms.Path(
        process.Zprime2muAnalysisSequence *
##        process.dyDimuons *
##        process.dyz * 
#adam        process.allGenSim *
        process.genMuons *
        process.genMuonsStoyan *
##         process.genMuonsBin1 *
##         process.genMuonsBin2 *
##         process.genMuonsBin3 *
##         process.genMuonsBin4 *
##         process.genMuonsBin5 *
##         process.genMuonsBin6 *
##         process.genMuonsBin7 *
         process.genMuonsBin8 #*
##         process.genMuonsBin9 *
##         process.genMuonsBin10 *
##         process.genMuonsBin11 *
##         process.genMuonsBin12 *
##         process.genMuonsBin13 *
##         process.genMuonsBin14
        )
    
#######################
##$outputFileName

trigger_match = ' && (' \
                '!triggerObjectMatchesByPath("HLT_Mu9").empty() || ' \
                '!triggerObjectMatchesByPath("HLT_Mu11").empty() || ' \
                '!triggerObjectMatchesByPath("HLT_Mu15_v1").empty()' \
                ')'

## trigger_match = ' && (' \
##                 '!triggerObjectMatchesByPath("HLT_DoubleMu3").empty() || ' \
##                 '!triggerObjectMatchesByPath("HLT_DoubleMu3_v2").empty()  ' \
##                 ')'



loose_cut = 'isGlobalMuon && isTrackerMuon && ' \
            'innerTrack.pt > 7. && ' \
            'abs(innerTrack.eta) < 2.4 && ' \
            'abs(dB) < 0.2 && ' \
            '(isolationR03.sumPt + isolationR03.hadEt) / innerTrack.pt < 0.15 && ' \
            'globalTrack.hitPattern.numberOfValidTrackerHits > 10 && ' \
            'globalTrack.hitPattern.numberOfValidPixelHits > 0 && ' \
            'globalTrack.hitPattern.numberOfValidMuonHits > 0 && ' \
            'numberOfMatches > 1 && ' \
            'globalTrack.normalizedChi2 < 10 ' # + trigger_match
       

tight_cut = 'abs(innerTrack.eta) < 2.1 && ' \
            'innerTrack.pt > 16. ' + trigger_match


process.allDimuons.loose_cut = loose_cut
process.allDimuons.tight_cut = tight_cut

process.leptons.muon_track_for_momentum = "global"
process.leptons.muon_cuts = loose_cut
#process.TFileService.fileName = "acc_skim_Zmass__lowerMasses_Zmumu.root"
#process.TFileService.fileName = "acc_skim_Zmass_Zpeak_Zmumu.root"
#process.TFileService.fileName = "acc_skim_diMass_200up_Zmumu.root"
