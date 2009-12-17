def customise(process):
    from Workspace.MuonRecoTreeUtility.muonRecoTreeUtility_cff import insertMRTU
    insertMRTU(process)
    return (process)
