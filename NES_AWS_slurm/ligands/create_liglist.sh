#!/bin/bash

for file in *.pdb; do
filename=$(basename "$file" .pdb)
echo "$filename" >> liglist.txt
done
