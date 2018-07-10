#!/bin/bash
ls
# check to see if all required arguments were provided
if [ $# -eq 1 ]; then
    # assign the provided argument to variables
    input_data=$1
else
    # assign the default value to a variable
    input_data="../data/proto.txt"
fi

g++ -o main main.cpp 
echo "Running helloworld with argument $input_data:"
./main -file proto.txt -circulant 500 -upGirth 8 -emd 20 -seed 123 -numberOfMatrices 1 -girth 8 
