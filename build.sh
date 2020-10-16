
# Rust core build
# cd spice21 
# cargo build 
# cd .. 

# Protobuffers
protoc -I=spice21/src --python_out=spice21py/spice21py --js_out=spice21js/spice21js spice21/src/*.proto 

# Publishing
# cd spice21 && cargo publish
# cd spice21py && maturin publish 
# cd spice21js && ???

