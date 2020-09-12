
cd spice21 
cargo build 
cd .. 
protoc -I=spice21/src --python_out=spice21py/spice21 --js_out=spice21js/spice21 spice21/src/spice21.proto
