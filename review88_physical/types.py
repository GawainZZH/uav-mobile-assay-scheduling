"""Immutable data for singleton sorties with complete service reservations."""
from dataclasses import dataclass
from fractions import Fraction
from math import ceil


def rational(value):
    return value if isinstance(value, Fraction) else Fraction(str(value))


def _integer(value, minimum, name):
    if type(value) is not int or value < minimum:
        raise ValueError(f'{name} must be an integer >= {minimum}')


@dataclass(frozen=True)
class Job:
    xy: tuple
    lo: int
    hi: int
    shelf: int
    deadline: int
    weight: object = 1
    name: str = ''

    def __post_init__(self):
        object.__setattr__(self, 'xy', tuple(self.xy))
        object.__setattr__(self, 'weight', rational(self.weight))
        for key in ('lo', 'hi', 'shelf'):
            _integer(getattr(self, key), 0, key)
        _integer(self.deadline, 1, 'deadline')
        if self.lo > self.hi or self.weight <= 0:
            raise ValueError('job window or weight')


@dataclass(frozen=True)
class Instance:
    name: str
    H: int
    pads: tuple
    jobs: tuple
    groups: int = 1
    drones: int = 1
    batteries: int = 6
    stations: int = 1
    sampling_duration: int = 1
    assay_duration: int = 6
    ground_duration: int = 4
    handoff_delay: int = 0
    endurance: int = 30
    flight_fixed: object = 2
    air_cost: object = 1
    road_cost: object = 1
    air_speed: object = 2
    flight_overhead: int = 0
    flight_times: object = None
    road_times: object = None

    def __post_init__(self):
        object.__setattr__(self, 'pads', tuple(tuple(p) for p in self.pads))
        object.__setattr__(self, 'jobs', tuple(self.jobs))
        for key in ('H', 'groups', 'drones', 'batteries', 'stations',
                    'sampling_duration', 'assay_duration', 'endurance'):
            _integer(getattr(self, key), 1, key)
        _integer(self.ground_duration, 0, 'ground_duration')
        _integer(self.handoff_delay, 0, 'handoff_delay')
        _integer(self.flight_overhead, 0, 'flight_overhead')
        if not self.pads or len(set(self.pads)) != len(self.pads):
            raise ValueError('distinct nonempty pads required')
        dimension = len(self.pads[0])
        if any(len(p) != dimension for p in self.pads) or any(len(j.xy) != dimension for j in self.jobs):
            raise ValueError('coordinate dimension')
        if any(j.hi >= self.H or j.deadline > self.H for j in self.jobs):
            raise ValueError('job outside horizon')
        for key in ('flight_fixed', 'air_cost', 'road_cost', 'air_speed'):
            value = rational(getattr(self, key))
            if value < 0 or (key == 'air_speed' and value == 0):
                raise ValueError('nonnegative costs and positive speed required')
            object.__setattr__(self, key, value)
        for key, size in (('flight_times', len(self.pads)+len(self.jobs)),
                          ('road_times', len(self.pads))):
            table = getattr(self, key)
            if table is not None:
                table = tuple(tuple(row) for row in table)
                if len(table) != size or any(len(row) != size for row in table):
                    raise ValueError(f'{key} dimensions')
                for i, row in enumerate(table):
                    for j, value in enumerate(row):
                        _integer(value, int(i != j), key)
                object.__setattr__(self, key, table)

    def xy(self, node):
        return self.pads[node] if node < len(self.pads) else self.jobs[node-len(self.pads)].xy

    def fly(self, a, b):
        if self.flight_times is not None:
            return self.flight_times[a][b]
        distance = sum(abs(rational(x)-rational(y)) for x, y in zip(self.xy(a), self.xy(b)))
        return max(1, ceil(distance/self.air_speed)) + self.flight_overhead

    def road(self, p, q):
        if p == q:
            return 0
        if self.road_times is not None:
            return self.road_times[p][q]
        distance = sum(abs(rational(x)-rational(y)) for x, y in zip(self.pads[p], self.pads[q]))
        if distance.denominator != 1:
            raise ValueError('default road reservations require integer coordinate distances')
        return max(1, int(distance))

    def latest_start(self, j, sample):
        job = self.jobs[j]
        return min(self.H-self.assay_duration, sample+job.shelf, job.deadline-self.assay_duration)


@dataclass(frozen=True, order=True)
class Route:
    p: int
    launch: int
    q: int
    recovery: int
    job: int
    sample: int

    @property
    def visits(self):
        return ((self.job, self.sample),)

    @property
    def context(self):
        return self.p, self.launch, self.q, self.recovery

    @property
    def launch_pad(self):
        return self.p

    @property
    def return_pad(self):
        return self.q

    @property
    def sample_start(self):
        return self.sample

    def window(self, inst):
        return self.recovery+inst.handoff_delay, inst.latest_start(self.job, self.sample)

    def cost(self, inst):
        return inst.flight_fixed + inst.air_cost*(self.recovery-self.launch) - inst.jobs[self.job].weight*self.sample
