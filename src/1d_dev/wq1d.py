from __future__ import print_function, absolute_import, division
import _wq1d
import f90wrap.runtime
import logging
import numpy

def wq1d(istl_, is2tl_):
    """
    wq1d(istl_, is2tl_)
    
    
    Defined at wq1d.f lines 1-226
    
    Parameters
    ----------
    istl_ : int
    is2tl_ : int
    
    """
    _wq1d.f90wrap_wq1d(istl_=istl_, is2tl_=is2tl_)

