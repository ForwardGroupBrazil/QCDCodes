import FWCore.ParameterSet.Config as cms

def customise(process):
    import UserCode.L3Switches.Switches as switch
    SwitchToBaseline(process)
    #SwitchToBaselinePP(process)
    #SwitchToOIState(process)
    #SwitchToOIHit(process)
    #SwitchToAllCombined(process)
    #SwitchToComboSeeds(process)
                            
    return(process)
