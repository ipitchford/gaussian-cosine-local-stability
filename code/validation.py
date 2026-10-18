"""Fail-closed finite-domain checks, active under optimized Python."""
import numpy as np

class NumericalDomainError(ValueError):
    """A valid mathematical frame exceeds this numerical implementation's coverage."""

def require(condition,message):
    if not condition:raise ValueError(message)

def finite(value,label):
    require(np.all(np.isfinite(value)),label+' must be finite')
    return value

def frame_input(U,exact_four=False):
    require(not np.iscomplexobj(U),'Only real frames are supported')
    U=np.asarray(U,dtype=float)
    require(U.ndim==2 and (U.shape[0]==4 if exact_four else 1<=U.shape[0]<=4) and 1<=U.shape[1]<=4,'Unsupported frame shape')
    finite(U,'frame')
    scales=np.max(np.abs(U),axis=1)
    require(np.all(scales>0),'Rows must be nonzero')
    U=U/scales[:,None]
    U/=np.linalg.norm(U,axis=1)[:,None]
    return finite(U,'normalized frame')
