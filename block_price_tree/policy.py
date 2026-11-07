"""Point inspection and append-only learning, independent of native callbacks."""
import time
from fractions import Fraction as F
from block_flow_solver.root import check
from block_flow_solver.bands import exact,update
from block_flow_solver.matrix import band_rows


class Inspector:
    def __init__(self,m,query,validate,reconstruct,deadline):
        self.m=m;self.query=query;self.validate=validate;self.reconstruct=reconstruct;self.deadline=deadline
        self.cursor=len(m.rows);self.pending=[];self.queried=set();self.events=[];self.inserted=0
        self.error=None;self.reason=None;self.best=None;self.point=None;self.checks=0;self.validated=0

    def capture(self):
        if len(self.m.rows)<self.cursor:raise ValueError('learned matrix must be append-only')
        self.pending.extend(range(self.cursor,len(self.m.rows)));self.cursor=len(self.m.rows)

    def flush(self,insert):
        while self.pending:
            index=self.pending[0];insert(self.m.rows[index],index)
            self.pending.pop(0);self.inserted+=1

    def inspect(self,x):
        self.checks+=1
        try:
            if time.perf_counter()>=self.deadline:self.reason='global deadline';return 'unresolved'
            self.capture()
            if self.pending:return 'cuts'
            try:check(self.m,x)
            except ValueError:return 'linear_infeasible'
            if any(abs(x[j]-round(x[j]))>1e-6 for j,flag in enumerate(self.m.integer)if flag):return 'fractional'
            selected=[i for i,js in self.m.family_columns.items()if sum(x[j]for j in js)>.5]
            unresolved=False
            for i in selected:
                b=self.m.bands[i]
                if exact(b):continue
                if i in self.queried:unresolved=True;continue
                self.queried.add(i);old=dict(b);start=len(self.m.rows)
                try:
                    r=self.query(i);cross_rows=len(self.m.rows)-start
                    update(b,r,lambda plan:self.validate(i,plan))
                    if b['infeasible']or b['lower']!=old['lower']or b['upper']!=old['upper']:
                        self.m.rows.extend(band_rows(self.m,i))
                    self.events.append(dict(family=i,record=r,cross_family_rows=cross_rows,new_rows=len(self.m.rows)-start))
                finally:self.capture()
                # All new rows must reach the live tree before another query at
                # this point, even if only another family's interval improved.
                if self.pending:return 'cuts'
                if not exact(b):unresolved=True
            if unresolved:
                self.reason='selected family unresolved after once-only query';return 'unresolved'
            result=self.reconstruct(self.m,x);value=F(str(result['objective']))
            if abs(float(value)-sum(float(c)*v for c,v in zip(self.m.cost,x)))>1e-6:
                raise ValueError('physical cost differs from integer master point')
            self.validated+=1
            if self.best is None or value<F(str(self.best['objective'])):self.best=result;self.point=list(x)
            return 'feasible'
        except TimeoutError as exc:self.reason=str(exc);return 'unresolved'
        except Exception as exc:self.error=repr(exc);return 'error'
