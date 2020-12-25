""" 
Spice21 Python 
"""
from typing import List, Union, Any, Dict, Optional

# Import protobuf-generated content
from . import protos
from .protos import Capacitor, Resistor, Diode, Mos, Isrc, Vsrc, MosType
from .protos import TranOptions
from .protos import Bsim4Model, Bsim4InstParams

# Import from the core Rust interface
from .spice21py import health


def dcop(arg: protos.Op) -> Dict[str, float]:
    """ DC Operating Point """
    from .spice21py import _dcop

    enc = bytes(arg)
    rb = _dcop(enc)
    rv = protos.OpResult().parse(rb)
    return rv.vals


def tran(
    ckt: protos.Circuit, opts: Optional[protos.TranOptions] = None
) -> Dict[str, List[float]]:
    """ Transient Analysis """
    from .spice21py import _tran

    # Serialize our circuit 
    ckt_enc = ckt.SerializeToString()
    # Serialize options, if provided. Otherwise we provide empty `bytes`. 
    opts_enc = bytes() if opts is None else opts.SerializeToString()
    return _tran(ckt_enc, opts_enc)


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
        return ckt.comps.append(Instance(d=comp))
    if isinstance(comp, Capacitor):
        return ckt.comps.append(Instance(c=comp))
    if isinstance(comp, Resistor):
        return ckt.comps.append(Instance(r=comp))
    if isinstance(comp, Mos):
        return ckt.comps.append(Instance(m=comp))
    if isinstance(comp, Isrc):
        return ckt.comps.append(Instance(i=comp))
    if isinstance(comp, Vsrc):
        return ckt.comps.append(Instance(v=comp))

    # Definitions
    if isinstance(comp, protos.Bsim4Model):
        return ckt.defs.append(protos.Def(bsim4model=comp))
    if isinstance(comp, protos.Bsim4InstParams):
        return ckt.defs.append(protos.Def(bsim4inst=comp))
    if isinstance(comp, protos.Mos1Model):
        return ckt.defs.append(protos.Def(mos1model=comp))
    if isinstance(comp, protos.Mos1InstParams):
        return ckt.defs.append(protos.Def(mos1inst=comp))

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
