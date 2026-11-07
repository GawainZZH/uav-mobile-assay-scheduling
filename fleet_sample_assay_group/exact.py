"""Exact common dyadic scale for represented physical and dual coefficients."""
import math
from fractions import Fraction
from fleet_sample_assay_activated.kernel import down_float


def rational(v):
    if not math.isfinite(float(v)):raise ValueError('finite pricing coefficients required')
    return Fraction(v) if isinstance(v,int) else Fraction(float(v))


class Scale:
    def __init__(self,values):
        fractions=[rational(x) for x in values]
        self.denominator=max((x.denominator for x in fractions),default=1)
        if self.denominator&(self.denominator-1):raise ValueError('only dyadic input coefficients supported')

    def integer(self,v):
        f=rational(v)
        if self.denominator%f.denominator:raise ValueError('unregistered pricing coefficient')
        return f.numerator*(self.denominator//f.denominator)

    def result(self,v):
        f=Fraction(v,self.denominator)
        return dict(value=float(f),lower_bound=down_float(f),exact_value=(str(f.numerator),str(f.denominator)))
