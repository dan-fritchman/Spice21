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


def dcop(
    arg: Optional[Union[protos.Circuit, protos.Op]] = None,
    *,
    ckt: Optional[protos.Circuit] = None,
    opts: Optional[protos.SimOptions] = None,
) -> Dict[str, float]:
    """ DC Operating Point """
    from .spice21py import _dcop

    # Discern whether our primary argument is an `Op` or a `Circuit`
    if arg is None:
        op = protos.Op(ckt=ckt, opts=opts)
    elif isinstance(arg, protos.Circuit):
        op = protos.Op(ckt=arg, opts=opts)
    elif isinstance(arg, protos.Op):
        op = arg
    else:
        raise TypeError

    # Encode
    enc = op.SerializeToString()
    # Do the real work
    rb = _dcop(enc)
    # Decode results
    rv = protos.OpResult()
    rv.ParseFromString(rb)
    # And return something in nicer dict-form
    return dict(rv.vals)


def tran(
    arg: Optional[Union[protos.Circuit, protos.Tran]] = None,
    *,
    ckt: Optional[protos.Circuit] = None,
    opts: Optional[protos.SimOptions] = None,
    args: Optional[protos.TranOptions] = None,
) -> Dict[str, List[float]]:
    """ Transient Analysis """
    from .spice21py import _tran

    # Discern whether our primary argument is a `Tran` or a `Circuit`
    if arg is None:
        tr = protos.Tran(ckt=ckt, opts=opts, args=args)
    elif isinstance(arg, protos.Circuit):
        tr = protos.Tran(ckt=arg, opts=opts, args=args)
    elif isinstance(arg, protos.Tran):
        tr = arg
    else:
        raise TypeError

    # Encode
    enc = tr.SerializeToString()
    # Do the real work
    rb = _tran(enc)
    # Decode results
    rv = protos.TranResult()
    rv.ParseFromString(rb)
    # And return something in nicer dict-form
    return {name: list(arr.vals) for name, arr in dict(rv.vals).items()}


def ac(
    arg: Optional[Union[protos.Circuit, protos.Ac]] = None,
    *,
    ckt: Optional[protos.Circuit] = None,
    opts: Optional[protos.SimOptions] = None,
    args: Optional[protos.AcOptions] = None,
) -> Dict[str, List[complex]]:
    """ AC Analysis """
    from .spice21py import _ac

    # Discern whether our primary argument is a `Tran` or a `Circuit`
    if arg is None:
        ac = protos.Ac(ckt=ckt, opts=opts, args=args)
    elif isinstance(arg, protos.Circuit):
        ac = protos.Ac(ckt=arg, opts=opts, args=args)
    elif isinstance(arg, protos.Ac):
        ac = arg
    else:
        raise TypeError

    # Encode
    enc = ac.SerializeToString()
    # Do the real work
    rb = _ac(enc)
    # Decode results
    rv = protos.AcResult()
    rv.ParseFromString(rb)
    # And return something in nicer dict-form
    return {
        name: [complex(c.re, c.im) for c in list(arr.vals)]
        for name, arr in dict(rv.vals).items()
    }


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
