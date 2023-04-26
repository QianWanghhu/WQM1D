#! urs/bin/env python

import numpy as np

def caltran():
    """
    Function for calculating the physical transport of dissolved constituent.

    Returns:
    =========
    WQV: np.ndarray of (L, ND), water quality concentrations

    """
    
    pass



def calfqc():
    """
    CALCULATES MASS SOURCES AND SINKS ASSOCIATED WITH CONSTANT AND TIME SERIES 
    INFLOWS AND OUTFLOWS; CONTROL STRUCTURE INFLOWS AND OUTLOWS; WITHDRAWAL AND RETURN STRUCTURE  
    OUTFLOWS; AND  EMBEDED CHANNEL INFLOWS AND OUTFLOWS 

    Returns:
    =========
    FQC: np.ndarray of (L, ND), mass of external sources/sinks
    """
    pass

def calad(LM, U, Con1):
    """
    Calculate advective-diffusive processes of dissolved or suspended constituents.
    Parameters:
    ============
    LM: 

    Returns:
    =========

    """
    # TODO: Advective fluxes
    FU = np.array(L)
    for L in LM:
        FU[L] = U[L] * Con1[L]

    # TODO: Diffusive transport
    
    pass

