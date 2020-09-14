from spice21py import *


def test_health():
    """ Test the core library is "there". """
    a = health()
    assert a == 'alive'
    print(f"spice21 core is {a}")


def test_dcop1():
    c = circuit([
        Diode(p="1"),
        Resistor(p="1", g=1e3),
        Capacitor(p="1"),
        Isrc(p="1", dc=1.0),
    ])
    print(dcop(c))


def test_dcop2():
    for v in (5, 6, 7, 8):
        c = circuit(Vsrc(name='vb', dc=0.1 * v, p='a'), Diode(p='a'))
        print(dcop(c))
