#!/bin/bash

cmake . 
make 
ctest -VV
