def customise(process):
    from Workspace.MuonRecoTreeUtility.muonRecoTreeUtilityForData_cff import insertMRTU
    insertMRTU(process)
    return (process)
