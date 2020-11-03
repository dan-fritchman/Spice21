import spice21py
from spice21py import circuit, Resistor, Capacitor, Isrc, Vsrc, Diode


def test_health():
    """ Test the core library is "there". """
    from spice21py import health

    a = health()
    assert a == "alive"
    print(f"spice21py core is {a} from {spice21py}")


def test_dcop1():
    from spice21py import dcop

    c = circuit([Resistor(p="1", g=1e-3), Capacitor(p="1"), Isrc(p="1", dc=1e-3),])
    res = dcop(c)
    assert isinstance(res, dict)
    assert res["1"] == 1.0


def test_dcop2():
    from spice21py import dcop

    for v in (5, 6, 7, 8):
        c = circuit(Vsrc(name="vb", dc=0.1 * v, p="a"), Diode(p="a"))
        res = dcop(c)
        assert res["a"] == 0.1 * v


def test_tran1():
    from spice21py import tran

    c = circuit([Resistor(p="1", g=1e-3), Capacitor(p="1"), Isrc(p="1", dc=1e-3),])
    res = tran(c)
    assert isinstance(res, dict)
    assert isinstance(res["1"], list)
    for v in res["1"]:
        assert v == 1.0


def test_ac1():
    from spice21py import ac

    c = circuit(
        Resistor(p="inp", n="out", g=1e-3),
        Capacitor(p="out", c=1e-9),
        Vsrc(p="inp", dc=1e-3, acm=1),
    )
    res = ac(c)
    assert isinstance(res, dict)
    assert isinstance(res["inp"], list)
    assert isinstance(res["out"], list)
    for v in res["inp"]:
        assert isinstance(v, complex)
        assert v == complex(1.0, 0.0)
    # Checking the output is decreasing (low-pass)
    for k in range(1, len(res["out"])):
        assert abs(res["out"][k]) <= abs(res["out"][k - 1])


def test_json():
    """ Test a JSON round-trip via Protobuf """
    from google.protobuf import json_format
    from ..protos import Circuit

    c = circuit(
        Resistor(p="inp", n="out", g=1e-3),
        Capacitor(p="out", c=1e-9),
        Vsrc(p="inp", dc=1e-3, acm=1),
    )
    j = json_format.MessageToJson(c)
    assert isinstance(j, str)
    c2 = json_format.Parse(j, Circuit())
    assert isinstance(c2, Circuit)
    assert c == c2


def test_bsim4_model():
    from .. import Bsim4Model, MosType

    b = Bsim4Model()
    assert b.mos_type == MosType.NMOS
    b = Bsim4Model(mos_type=MosType.PMOS)
    assert b.mos_type == MosType.PMOS


def test_bsim4_inst_params():
    from .. import Bsim4InstParams

    p = Bsim4InstParams()
    assert isinstance(p, Bsim4InstParams)
    p.l.value = 1e-6
    p.w.value = 1e-6
    p.nf.value = 2
    assert p.l.value == 1e-6
    assert p.w.value == 1e-6
    assert p.nf.value == 2


def test_bsim4_ckt():
    from .. import Bsim4Model, Bsim4InstParams, Mos, circuit

    b = Mos(name="inst1")
    b.ports.g = "vdd"
    b.ports.d = "vdd"
    b.model = "default"
    b.params = "inst"
    c = circuit(
        Bsim4Model(name="default"),
        Bsim4InstParams(name="inst"),
        b,
        Vsrc(p="vdd", dc=1.0),
    )
    from .. import dcop

    res = dcop(c)
    # FIXME: checks


def test_bsim4_ro():
    from .. import circuit 
    from .. import Bsim4Model, Bsim4InstParams, Mos, MosType, TranOptions, Capacitor

    def inverter(name: str, inp: str, out: str) -> list:
        p = Mos(name=f"{name}p", model="pmos", params="default")
        p.ports.g = inp
        p.ports.d = out
        p.ports.s = "vdd"
        p.ports.b = "vdd"
        n = Mos(name=f"{name}n", model="nmos", params="default")
        n.ports.g = inp
        n.ports.d = out
        n.ports.s = "vss"
        n.ports.b = "vss"
        c = Capacitor(name=f"{name}c", p=out, n="vss", c=1e-12)
        return [p, n, c]

    nstg = 3
    insts = []
    for k in range(nstg):
        insts.extend(inverter(name=f"stg{k}", inp=f"ring{k}", out=f"ring{k+1}"))

    ckt = circuit(
        Bsim4Model(name="nmos", mos_type=MosType.NMOS),
        Bsim4Model(name="pmos", mos_type="PMOS"),
        Bsim4InstParams(name="default"),
        Vsrc(name="vvdd", p="vdd", n="vss", dc=1.0),
        Vsrc(name="vvss", p="vss", dc=0),
        Vsrc(name="vvfb", p="ring0", n=f"ring{nstg}", dc=0),
        *insts,
    )
    opts = TranOptions(tstop=1e-7, tstep=1e-12)
    opts.ic['ring0'] = 0.0

    from .. import tran

    res = tran(ckt, opts)
    assert isinstance(res, dict)
    # FIXME: checks


def test_mos1_ro():
    from .. import circuit 
    from .. import Mos, MosType, TranOptions, Capacitor
    from ..protos import Mos1Model, Mos1InstParams

    def inverter(name: str, inp: str, out: str) -> list:
        p = Mos(name=f"{name}p", model="pmos", params="default")
        p.ports.g = inp
        p.ports.d = out
        p.ports.s = "vdd"
        p.ports.b = "vdd"
        n = Mos(name=f"{name}n", model="nmos", params="default")
        n.ports.g = inp
        n.ports.d = out
        n.ports.s = "vss"
        n.ports.b = "vss"
        c = Capacitor(name=f"{name}c", p=out, n="vss", c=1e-12)
        return [p, n, c]

    nstg = 3
    insts = []
    for k in range(nstg):
        insts.extend(inverter(name=f"stg{k}", inp=f"ring{k}", out=f"ring{k+1}"))

    ckt = circuit(
        Mos1Model(name="nmos", mos_type=MosType.NMOS),
        Mos1Model(name="pmos", mos_type="PMOS"),
        Mos1InstParams(name="default"),
        Vsrc(name="vvdd", p="vdd", n="vss", dc=1.0),
        Vsrc(name="vvss", p="vss", dc=0),
        Vsrc(name="vvfb", p="ring0", n=f"ring{nstg}", dc=0),
        *insts,
    )
    opts = TranOptions(tstop=1e-7, tstep=1e-12)
    opts.ic['ring0'] = 0.0

    from .. import tran

    res = tran(ckt, opts)
    assert isinstance(res, dict)
    # FIXME: checks
