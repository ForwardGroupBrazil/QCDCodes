process.load("DQMServices.Components.EDMtoMEConverter_cff")

process.load("DQMServices.Components.DQMEnvironment_cfi")

process.postProcessor = cms.EDFilter("PostProcessor",
    commands = cms.vstring("RecoMuonV/MultiTrack R v_cotThetares_eta \'Res(cot(#theta)) vs #eta\' \'cotThetares_vs_eta\'", 
        "RecoMuonV/MultiTrack R v_cotThetares_pt \'Res(cot(#theta)) vs pt\' \'cotThetares_vs_pt\'", 
        "RecoMuonV/MultiTrack R v_dxy_Pull_eta \'Pull(dxy) vs #eta\' \'dxypull_vs_eta\'", 
        "RecoMuonV/MultiTrack R v_dxy_Resolution_eta \'Res(dxy) vs #eta\' \'dxyres_vs_eta\'", 
        "RecoMuonV/MultiTrack R v_dxy_Pull_pt \'Pull(dxy) vs pt\' \'dxypull_vs_pt\'", 
        "RecoMuonV/MultiTrack R v_dxy_Resolution_pt \'Res(dxy) vs pt\' \'dxyres_vs_pt\'", 
        "RecoMuonV/MultiTrack R v_dz_Pull_eta \'Pull(dz) vs #eta\' \'dzpull_vs_eta\'", 
        "RecoMuonV/MultiTrack R v_dz_Resolution_eta \'Res(dz) vs #eta\' \'dzres_vs_eta\'", 
        "RecoMuonV/MultiTrack R v_dz_Resolution_pt \'Res(dz) vs pt\' \'dzres_vs_pt\'", 
        "RecoMuonV/MultiTrack R v_eta_Resolution_eta \'Res(#eta) vs #eta\' \'etares_vs_eta\'", 
        "RecoMuonV/MultiTrack R v_phi_Pull_eta \'Pull(#phi) vs #eta\' \'phipull_vs_eta\'", 
        "RecoMuonV/MultiTrack R v_phi_Resolution_eta \'Res(#phi) vs #eta\' \'phires_vs_eta\'", 
        "RecoMuonV/MultiTrack R v_phi_Resolution_pt \'Res(#phi) vs pt\' \'phires_vs_pt\'", 
        "RecoMuonV/MultiTrack R v_Pt_Pull_eta \'Pull(pt) vs #eta\' \'ptpull_vs_eta\'", 
        "RecoMuonV/MultiTrack R v_Pt_Resolution_eta \'Res(pt) vs #eta\' \'ptres_vs_eta\'", 
        "RecoMuonV/MultiTrack R v_Pt_Resolution_pt \'Res(pt) vs pt\' \'ptres_vs_pt\'", 
        "RecoMuonV/MultiTrack R v_theta_Pull_eta \'Pull(#theta) vs #eta\' \'thetapull_vs_eta\'", 
        "RecoMuonV/MultiTrack E v_Efficiency_eta \'Efficiency Plot vs #eta\' \'num_assoc(simToReco)_eta\' \'num_simul_eta\'", 
        "RecoMuonV/MultiTrack E v_Efficiency_pt \'Efficiency Plot vs pt\' \'num_assoc(simToReco)_pT\' \'num_simul_pT\'", 
        "RecoMuonV/MultiTrack E v_Efficiency_hit \'Efficiency Plot vs hits\' \'num_assoc(simToReco)_hit\' \'num_simul_hit\'", 
        "RecoMuonV/MultiTrack E v_Fake_eta \'Fake Plot vs #eta\' \'num_assoc(recoToSim)_eta\' \'num_reco_eta\'", 
        "RecoMuonV/MultiTrack E v_Fake_pt \'Fake Plot vs pt\' \'num_assoc(recoToSim)_pT\' \'num_reco_pT\'", 
        "RecoMuonV/MultiTrack E v_Fake_hit \'Fake Plot vs hits\' \'num_assoc(recoToSim)_hit\' \'num_reco_hit\'"),
    subDir = cms.untracked.string('RecoMuonV/MultiTrack'),
    outputFileName = cms.untracked.string('$outFileName_post.root')
)

process.p1 = cms.Path(process.EDMtoMEConverter*process.postProcessor)

process.EDMtoMEConverter.convertOnEndRun = True
process.DQMStore.referenceFileName = ''
process.dqmSaver.convention = 'RelVal'
process.dqmSaver.workflow = '/Validation/Test/RECO'

