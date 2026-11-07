"""Feasibility-gated near-integer normalization before exact carrier completion.

This adapter does not widen the original cost guard or change any MILP row.
Old source files and experiment snapshots remain untouched.
"""
from contextvars import ContextVar
import numpy as np

_rhs = ContextVar('carrier_actual_context_rhs', default=None)


def residual(m, x, rhs):
    eq = np.asarray(m.M.eq, dtype=bool)
    activity = m.A@x
    return float(max(np.max(abs(activity[eq]-rhs[eq]), initial=0),
                     np.max(activity[~eq]-rhs[~eq], initial=0),
                     np.max(-x, initial=0),
                     np.max(x-np.asarray(m.M.upper), initial=0)))


def normalize(m, vector, rhs=None, native_tolerance=1e-6, rounded_tolerance=1e-9):
    x = np.asarray(vector, dtype=float).copy()
    rhs = np.asarray(m.M.rhs if rhs is None else rhs, dtype=float)
    if not np.all(np.isfinite(x)):
        raise ValueError('nonfinite native primal point')
    raw_error = residual(m,x,rhs)
    if raw_error > native_tolerance:
        raise ValueError(('native primal residual before completion',raw_error))
    integer_distance = float(np.max(abs(x-np.rint(x)),initial=0))
    note = dict(changed=False,native_residual=raw_error,integer_distance=integer_distance)
    if integer_distance <= native_tolerance:
        rounded = np.rint(x)
        error = residual(m,rounded,rhs)
        note['rounded_residual'] = error
        if error <= rounded_tolerance:
            note.update(changed=bool(np.any(rounded!=x)),changed_entries=int(np.sum(rounded!=x)))
            return rounded,note
    return x,note


def install():
    """Apply the same reconstruction adapter to SK kernels and strong PQ MILPs."""
    from service_separator25_joint import model, kernel
    from revision24 import solver
    original_complete = model.complete
    original_solve = kernel.Kernel.solve
    notes = []
    def solve(self,start,end,consumed,deadline):
        token = _rhs.set((id(self.m),self.rhs(start,end,consumed)))
        try:
            return original_solve(self,start,end,consumed,deadline)
        finally:
            _rhs.reset(token)
    def complete(m,vector,start,end):
        current = _rhs.get()
        rhs = current[1] if current is not None and current[0]==id(m) else m.M.rhs
        x,note = normalize(m,vector,rhs)
        if note['changed']:
            notes.append(note)
        return original_complete(m,x,start,end)
    kernel.Kernel.solve = solve
    model.complete = complete
    solver.complete = complete
    return notes
