import FWCore.ParameterSet.Config as cms

def customise(process):
    import UserCode.L3Switches.Switches as switch
    #switch.SwitchToBaseline(process)
    #switch.SwitchToBaselinePP(process)
    #switch.SwitchToOIState(process)
    #switch.SwitchToOIHit(process)
    #switch.SwitchToAllCombined(process)
    #switch.SwitchToOICombined(process)
    switch.SwitchToIterative3(process)
    
    return(process)
