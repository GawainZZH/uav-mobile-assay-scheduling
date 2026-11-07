from contextlib import contextmanager
from service_separator25_joint import monolithic
from service_separator25_joint.interface import certify as conservative
from referee_boundary21_v2.interface import certify as joint_rule


@contextmanager
def boundary_rule(joint=True):
    saved=monolithic.certify;monolithic.certify=joint_rule if joint else conservative
    try:yield
    finally:monolithic.certify=saved
