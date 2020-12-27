set -e 

# Spice21 Publishing Script 

# Publish Rust crates, in dependency order 
# Thus far this works better as a recipe of steps, rather than a script, 
# Since each requires waiting for the prior crates to be publicly available. 
cd spice21int   && cargo publish && cd .. 
cd spice21procs && cargo publish && cd .. 
cd spice21      && cargo publish && cd .. 

# Publish Python Bindings 
cd spice21py    && ./build.sh && maturin publish && cd .. 

# JavaScript will have to wait 
# cd spice21js  && neon build && ...

