""" 
Spice21 Python 
"""
# Import and rename the protobuf-generated content
from . import protos
from .protos import Capacitor, Resistor, Diode, Mos
# Import the core interface
from .spice21py import *


def circuit(*args) -> protos.Circuit:
    """ Circuit generation from component list *args """
    from .spice21_pb2 import Circuit, Instance

    if len(args) == 1 and isinstance(args[0], (list, tuple)):
        comps = args[0]
    else:
        comps = args

    c = Circuit()
    for comp in comps:
        if isinstance(comp, Instance):
            c.comps.append(comp)
        elif isinstance(comp, Diode):
            c.comps.append(Instance(d=comp))
        elif isinstance(comp, Capacitor):
            c.comps.append(Instance(c=comp))
        elif isinstance(comp, Resistor):
            c.comps.append(Instance(r=comp))
        elif isinstance(comp, Mos):
            c.comps.append(Instance(mos=comp))
        else:
            raise TypeError(f"Invalid Circuit Component {comp}")
    return c
