
# Python Protobuf Compilation 
protoc -I=../spice21/protos --python_out=spice21py/protos ../spice21/protos/*.proto
# Sadly protoc doesn't seem to know how Python3 imports work 
2to3 -wn -f import spice21py/protos/*.py
