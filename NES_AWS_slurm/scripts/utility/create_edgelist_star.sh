#!/bin/bash

input_file="liglist.txt"
output_file="edgelist.txt"

sed 's/^/lig_4-99,/' "$input_file" > "$output_file"
