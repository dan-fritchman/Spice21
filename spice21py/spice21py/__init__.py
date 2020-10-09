""" 
Spice21 Python 
"""
from typing import List, Union, Any, Dict 

# Import protobuf-generated content
from . import protos
from .protos import Capacitor, Resistor, Diode, Mos, Isrc, Vsrc, Bsim4Model, Bsim4InstParams, MosType 

# Import from the core Rust interface
from .spice21py import health


def dcop(ckt: protos.Circuit) -> Dict[str, float]:
    """ DC Operating Point """
    from .spice21py import _dcop
    enc = ckt.SerializeToString()
    return _dcop(enc)


def tran(ckt: protos.Circuit) -> Dict[str, List[float]]:
    """ Transient Analysis """
    from .spice21py import _tran
    enc = ckt.SerializeToString()
    return _tran(enc)


def ac(ckt: protos.Circuit) -> Dict[str, List[complex]]:
    """ AC Analysis """
    from .spice21py import _ac
    enc = ckt.SerializeToString()
    res = _ac(enc)
    # Here we decode into Python `complex` values
    return {name: [complex(*v) for v in vals] for name, vals in res.items()}


def instance(comp: Any) -> protos.Instance:
    """ Create an Instance from Component `comp`.
    Raises a TypeError if `comp` is not a Component or Instance. """

    from .protos import Instance
    if isinstance(comp, Instance):
        return comp
    if isinstance(comp, Diode):
        return Instance(d=comp)
    if isinstance(comp, Capacitor):
        return Instance(c=comp)
    if isinstance(comp, Resistor):
        return Instance(r=comp)
    if isinstance(comp, Mos):
        return Instance(mos=comp)
    if isinstance(comp, Isrc):
        return Instance(i=comp)
    if isinstance(comp, Vsrc):
        return Instance(v=comp)

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
