#!/bin/bash


for matrix in $(ls ../matrices/large/*) ; do

echo $matrix
echo ""

./lif1d --spfile $matrix --time=0.05

echo ""

done
