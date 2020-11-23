
# Rust core build
cd spice21 
cargo build 
cd .. 

# Python Protobuf Compilation 
# Sadly protoc doesn't seem to know how Python3 imports work 
# https://github.com/protocolbuffers/protobuf/issues/1491
protoc -I=spice21/src --python_out=spice21py/spice21py spice21/src/*.proto 
2to3 -wn -f import spice21py/spice21py/*pb2* 

# JavaScript Protobuf Compilation 
cd spice21js 
yarn 
yarn protoc 
cd ..

# Publishing
# cd spice21 && cargo publish
# cd spice21py && maturin publish 
# cd spice21js && ???

