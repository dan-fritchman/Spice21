
# Rust core build
cd spice21 
cargo build 
cd .. 

cd spice21py
./build.sh
cd ..

# JavaScript Protobuf Compilation 
cd spice21js 
yarn 
yarn protoc 
cd ..

# Publishing
# cd spice21 && cargo publish
# cd spice21py && maturin publish 
# cd spice21js && ???

