"""Identical complete-fleet event pipeline for native and decomposed solvers."""
import time


def pipeline(module,operation,began):
    row=dict(detected=time.perf_counter()-began,reconstruction_started=time.perf_counter()-began)
    saved=module.validate_fleet
    def validate(*args,**kwargs):
        row['reconstructed']=time.perf_counter()-began
        result=saved(*args,**kwargs)
        row['validated']=time.perf_counter()-began
        if not result['valid']:raise ValueError(('physical validation',result['errors']))
        return result
    module.validate_fleet=validate
    try:result=operation()
    finally:module.validate_fleet=saved
    if 'validated'not in row:raise ValueError('missing physical validation phase')
    row['released']=time.perf_counter()-began;row['seconds']=row['released']
    return result,row
