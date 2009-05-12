import FWCore.ParameterSet.Config as cms

def customise(process):
    import UserCode.L3Switches.Switches as switch
    #switch.SwitchToOIState(process)
    #switch.SwitchToOIHit(process)
    #switch.SwitchToBaselinePP(process)
    #switch.SwitchToAllCombined(process)
    switch.SwitchToComboSeeds(process)
    
    return(process)
