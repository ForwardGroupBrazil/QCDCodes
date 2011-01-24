def customise(process):
    from Workspace.MuonRecoTreeUtility.muonRecoTreeUtility_cff import insertMHTU
    insertMHTU(process)
    return (process)
