from dataclasses import dataclass
from fractions import Fraction
import numpy as np


@dataclass
class Price:
    weights: np.ndarray
    constant: int
    scale: int

    def evaluate(self,rhs):
        rhs=np.asarray(rhs)
        if np.max(abs(rhs-np.rint(rhs)),initial=0)>0:raise ValueError('integer RHS required')
        return Fraction(sum(int(a)*int(b) for a,b in zip(self.weights,rhs))+self.constant,self.scale)


def price(A,c,upper,equalities,dual,scale=2**24):
    A=A.tocsc();c=np.asarray(c);upper=np.asarray(upper)
    if np.max(abs(A.data-np.rint(A.data)),initial=0)>0:raise ValueError('integer matrix required')
    if np.max(abs(c-np.rint(c)),initial=0)>0 or np.max(abs(upper-np.rint(upper)),initial=0)>0:
        raise ValueError('integer costs and finite caps required')
    y=np.rint(np.asarray(dual)*scale).astype(np.int64)
    y[~np.asarray(equalities,dtype=bool)]=np.minimum(0,y[~np.asarray(equalities,dtype=bool)])
    constant=0
    for j in range(A.shape[1]):
        lo,hi=A.indptr[j:j+2]
        rc=int(round(c[j]))*scale-sum(int(A.data[k])*int(y[A.indices[k]]) for k in range(lo,hi))
        constant+=int(upper[j])*min(0,rc)
    return Price(y,constant,scale)
