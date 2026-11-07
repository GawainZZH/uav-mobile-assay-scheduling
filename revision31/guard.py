"""Windows job owner and independent deadline process; imports no solver.

Job semantics: https://learn.microsoft.com/en-us/windows/win32/procthread/job-objects
The guard joins its own kill-on-close job before spawning the worker. Worker
descendants inherit membership, but no child inherits the owning job handle.
"""
import ctypes
from ctypes import wintypes as W
import math
import os
import subprocess
import sys
import time


def guarded_command(command,deadline,parent_pid):
    return [sys.executable,'-B','-X','utf8','-m','revision31.guard',repr(deadline),str(parent_pid),'--',*command]


class Basic(ctypes.Structure):
    _fields_=[('process_time',ctypes.c_int64),('job_time',ctypes.c_int64),('flags',W.DWORD),
              ('min_ws',ctypes.c_size_t),('max_ws',ctypes.c_size_t),('active',W.DWORD),
              ('affinity',ctypes.c_size_t),('priority',W.DWORD),('scheduling',W.DWORD)]

class Counters(ctypes.Structure):
    _fields_=[(name,ctypes.c_uint64)for name in ('read_ops','write_ops','other_ops','read_bytes','write_bytes','other_bytes')]

class Extended(ctypes.Structure):
    _fields_=[('basic',Basic),('io',Counters),('process_memory',ctypes.c_size_t),
              ('job_memory',ctypes.c_size_t),('peak_process',ctypes.c_size_t),('peak_job',ctypes.c_size_t)]


def supervise(command,deadline,parent_pid):
    if os.name!='nt':raise RuntimeError('Windows job guard required')
    if not math.isfinite(deadline)or time.perf_counter()>=deadline:os._exit(124)
    api=ctypes.WinDLL('kernel32',use_last_error=True)
    signatures={
        'CreateJobObjectW':([ctypes.c_void_p,W.LPCWSTR],W.HANDLE),
        'SetInformationJobObject':([W.HANDLE,ctypes.c_int,ctypes.c_void_p,W.DWORD],W.BOOL),
        'AssignProcessToJobObject':([W.HANDLE,W.HANDLE],W.BOOL),
        'GetCurrentProcess':([],W.HANDLE),
        'OpenProcess':([W.DWORD,W.BOOL,W.DWORD],W.HANDLE),
        'WaitForSingleObject':([W.HANDLE,W.DWORD],W.DWORD),
        'TerminateJobObject':([W.HANDLE,W.UINT],W.BOOL),
    }
    for name,(args,result)in signatures.items():
        func=getattr(api,name);func.argtypes=args;func.restype=result
    job=api.CreateJobObjectW(None,None)
    if not job:raise ctypes.WinError(ctypes.get_last_error())
    info=Extended();info.basic.flags=0x2000 # JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE
    if not api.SetInformationJobObject(job,9,ctypes.byref(info),ctypes.sizeof(info)):
        raise ctypes.WinError(ctypes.get_last_error())
    parent=api.OpenProcess(0x00100000,False,parent_pid) # SYNCHRONIZE
    if not parent:raise ctypes.WinError(ctypes.get_last_error())
    if not api.AssignProcessToJobObject(job,api.GetCurrentProcess()):
        raise ctypes.WinError(ctypes.get_last_error())
    # Owning handles deliberately remain open until this process exits. Closing
    # the guard (normally, by parent kill, or by parent death) reaps the job.
    child=subprocess.Popen(command,creationflags=subprocess.CREATE_NO_WINDOW,close_fds=True,
                           stdin=subprocess.DEVNULL,stdout=sys.stdout,stderr=sys.stderr,
                           env={**os.environ,'REVISION31_GUARD_DEADLINE':repr(deadline)})
    while True:
        code=child.poll()
        if code is not None:os._exit(code if 0<=code<256 else 125)
        parent_status=api.WaitForSingleObject(parent,0)
        if parent_status not in (0,258):raise ctypes.WinError(ctypes.get_last_error())
        if parent_status==0 or time.perf_counter()>=deadline:
            if not api.TerminateJobObject(job,124):raise ctypes.WinError(ctypes.get_last_error())
            os._exit(124)
        time.sleep(min(.01,max(0,deadline-time.perf_counter())))


if __name__=='__main__':
    try:
        if len(sys.argv)<6 or sys.argv[3]!='--':raise ValueError('deadline parent -- command required')
        supervise(sys.argv[4:],float(sys.argv[1]),int(sys.argv[2]))
    except BaseException:
        import traceback
        traceback.print_exc();sys.stderr.flush()
        os._exit(125)
