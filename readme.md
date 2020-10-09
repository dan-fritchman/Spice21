
# Spice21

SPICE for the 21st century 

---

![test](https://github.com/HW21/Spice21/workflows/test/badge.svg)

--- 

## Circuit Simulation Circa 2020 

Spice21 is a circuit simulation library akin to the original Berkeley 
*Simulation Program with Integrated Circuit Emphasis*,
designed around a set of principles for 21st century users: 

### Spice is a *library*, not a *program*. 

Circuit simulation almost always is (or should be) embedded in a larger program - 
whether iterating over circuit conditions or performing design exploration. 
Spice21 is designed for such embedding. 
Circuits are data in the program - not schematics entombed in an 80s-grade GUI 
or intractable netlist-languages. 

You may ask: which programming language do I need for this embedding? 
How about *whatever language you like*. 

Spice21 is implemented in [Rust](https://www.rust-lang.org/) 
and uses Google's [Protocol Buffers](https://developers.google.com/protocol-buffers) 
for most input and output data. 
A set of thin wrappers provide bindings to just about any language that supports Protobuf. 

### Don't Do Stupid (Old) Shit 

Spice was born a generation ago, and had a lifetime to accumulate mistakes. 
Throw them out. Spice21 includes no: 

**Stupid file formats**. 
Spice21 has no custom netlist, output, or any other dreamt-up file-formats 
that will require anyone to parse anything, ever. 
Circuits are defined in the Protobuf schema language. 
All output lives in open, popular data formats. 

**Outdated device models and analyses.** 
Transistor models will come in two flavors: the simplest, and the newest and most relevant. 
To date no open-source SPICE supports [BSIM-CMG](http://bsim.berkeley.edu/models/bsimcmg/), 
the industry-standard model for now-industry-standard FinFETs. 
Remedying this is a near-term concern. 

**Intractable options nobody ever uses.** 
Careers-worth of circuit design and simulation experience have taught us how to set these things. 
We're throwing out all the lessons un-learned along the way. 

### No Shortcuts

Spice21 doesn't take any shortcuts to simulating transistor-level circuits.
There are no table-based "FastSpice" tactics, circuit simplifications, 
or attempts to categorize or tailor to particular circuit families. 
This is a fully general-purpose circuit solver; 
throw it whatever pile of transistors you like. 

### Open Souce 

Spice21 is distributed under a permissive open-source license for all. 

## Installation 

### [Rust](https://crates.io/crates/spice21) 

Add [spice21](https://crates.io/crates/spice21) to the Cargo.toml of any Rust project. 

```toml
[dependencies]
spice21 = "0.1.4"
```

### [Python](https://pypi.org/project/spice21py/)

`pip install spice21py` 

Pip will detect whether your combination of OS and Python interpreter have 
pre-compiler wheel-distributions of Spice21. If not, it'll run a source-build 
which will require installing the [Rust compiler](https://www.rust-lang.org/tools/install). 
(This will generally take < 1 minute.)

### [JavaScript](https://www.npmjs.com/package/spice21js)

```
npm install spice21js
```

Or 

```
yarn add spice21js
```

## Status

Spice21 is in early-days development. 
Versions 0.1.x support transient, operating-point, and AC analysis of circuits including 
level-one MOS, passives, diodes, and independent sources. 

## Thanks 

Spice21 is a descendant of a long line of circuit simulators, 
many borne at the University of California, Berkeley. 
Several component models are adapted from their original Berkeley SPICE implementations. 

