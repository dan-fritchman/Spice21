
# Spice21 Build Script
# Including language bindings 
set -e 

# Rust core build, in root workspace-directory
cargo build 

# Python Bindings 
cd spice21py && ./build.sh && cd ..

# JavaScript Bindings 
cd spice21js && yarn install && yarn protoc && cd ..

