process.load("Configuration.StandardSequences.PostRecoGenerator_cff")

process.load("Configuration.StandardSequences.Validation_cff")

process.load("Configuration.StandardSequences.PostRecoMuon_cff")

process.load("Validation.Configuration.muonMultiValid_cff")

process.postreco_step = cms.Path(process.postreco_generator)
process.validation_step = cms.Path(process.validation)
process.postmuon_step = cms.Path(process.postreco_muon)
process.muonValidation_step = cms.Path(process.muonMultiValidSequence)

process.multiGlbTrackValidator.out = '$outFileName_dqm.root'
process.multiStaTrackValidator.out = '$outFileName_dqm.root'
process.multiGlbTrackValidatorPos.out = '$outFileName_dqm.root'
process.multiStaTrackValidatorPos.out = '$outFileName_dqm.root'
process.multiTkTrackValidator.out = '$outFileName_dqm.root'
