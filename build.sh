
# Rust core build
# cd spice21 
# cargo build 
# cd .. 

# Protobuf Compilation 
protoc -I=spice21/src --js_out=spice21js/spice21js spice21/src/*.proto 

# Sadly protoc doesn't seem to know how Python3 imports work 
# https://github.com/protocolbuffers/protobuf/issues/1491
protoc -I=spice21/src --python_out=spice21py/spice21py spice21/src/*.proto 
sed -i -r 's/^import (.+_pb2.*)/from . import \1/g' spice21py/spice21py/*_pb2*.py


# Publishing
# cd spice21 && cargo publish
# cd spice21py && maturin publish 
# cd spice21js && ???

