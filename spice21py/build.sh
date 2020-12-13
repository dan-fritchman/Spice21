
# Python (Better) Protobuf Compilation 
protoc -I=../spice21/protos --python_betterproto_out=spice21py/protos ../spice21/protos/*.proto
# Retain our __init__, which pulls these up to the 'protos' namespace 
# Better-Protobuf doesn't seem to want to
cp spice21py/protos/_init.py spice21py/protos/__init__.py 
