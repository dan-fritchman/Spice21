""" 
Spice21 Python 
"""
from typing import List, Union, Any, Dict 

# Import protobuf-generated content
from . import protos
from .protos import Capacitor, Resistor, Diode, Mos, Isrc, Vsrc, MosType 
from .protos import Bsim4, Bsim4Model, Bsim4InstParams

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


def add(ckt: protos.Circuit, comp: Any) -> None:
    """ Organize and Add Circuit Items. 
    Visit by type and append to component and definition-lists, 
    inserting our enum-layers along the way. 
    Raises a TypeError for anything that can't get into a `protos.Circuit`. """

    # Instances 
    from .protos import Instance
    if isinstance(comp, Instance):
        return ckt.comps.append(comp)
    if isinstance(comp, Diode):
        return ckt.comps.append( Instance(d=comp))
    if isinstance(comp, Capacitor):
        return ckt.comps.append( Instance(c=comp))
    if isinstance(comp, Resistor):
        return ckt.comps.append( Instance(r=comp))
    if isinstance(comp, Mos):
        return ckt.comps.append( Instance(m=comp))
    if isinstance(comp, Isrc):
        return ckt.comps.append( Instance(i=comp))
    if isinstance(comp, Vsrc):
        return ckt.comps.append( Instance(v=comp))
    if isinstance(comp, Bsim4):
        return ckt.comps.append( Instance(bsim4=comp))
    
    # Definitions 
    if isinstance(comp, protos.Bsim4Model):
        return ckt.defs.append(protos.Def(bsim4model=comp))
    if isinstance(comp, protos.Bsim4InstParams):
        return ckt.defs.append(protos.Def(bsim4inst=comp))
    
    # Anything else got here by mistake
    raise TypeError(f"Invalid Circuit Component {type(comp)}")


def circuit(*args) -> protos.Circuit:
    """ Circuit generation from component list *args """
    from .protos import Circuit, Instance

    _args = args
    if len(args) == 1 and isinstance(args[0], (list, tuple)):
        _args = args[0]

    c = Circuit()
    for arg in _args:
        add(c, arg)

    return c
