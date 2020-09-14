""" 
Spice21 Python 
"""
from typing import List, Union, Any
# Import and rename the protobuf-generated content
from . import protos
from .protos import Capacitor, Resistor, Diode, Mos, Isrc, Vsrc
# Import from the core Rust interface
from .spice21py import health


def dcop(ckt: protos.Circuit) -> List[float]:
    """ DC Operating Point """
    from .spice21py import _dcop
    enc = ckt.SerializeToString()
    return _dcop(enc)


def tran(ckt: protos.Circuit) -> List[List[float]]:
    """ Transient Analysis """
    from .spice21py import _tran
    enc = ckt.SerializeToString()
    return _tran(enc)


def instance(comp) -> protos.Instance:
    """ Create an Instance from Component `comp`.
    Raises a TypeError if `comp` is not a Component or Instance. """
    from .protos import Instance
    if isinstance(comp, Instance):
        return comp
    elif isinstance(comp, Diode):
        return Instance(d=comp)
    elif isinstance(comp, Capacitor):
        return Instance(c=comp)
    elif isinstance(comp, Resistor):
        return Instance(r=comp)
    elif isinstance(comp, Mos):
        return Instance(mos=comp)
    elif isinstance(comp, Isrc):
        return Instance(i=comp)
    elif isinstance(comp, Vsrc):
        return Instance(v=comp)
    else:
        raise TypeError(f"Invalid Circuit Component {type(comp)}")


def circuit(*args) -> protos.Circuit:
    """ Circuit generation from component list *args """
    from .protos import Circuit, Instance

    if len(args) == 1 and isinstance(args[0], (list, tuple)):
        comps = args[0]
    else:
        comps = args

    c = Circuit()
    for comp in comps:
        c.comps.append(instance(comp))

    return c
