#!/bin/bash

all_dirs=$(ls -d edge_lig_*_lig_*/state*/run*/)

results=$(find edge_lig_*_lig_*/state*/run*/ -name 'dhdl80*')

found_dirs=$(echo "$results" | awk -F'/' '{print $1"/"$2"/"$3"/"}' | sort | uniq)

not_found_dirs=$(comm -23 <(echo "$all_dirs" | sort) <(echo "$found_dirs" | sort))

echo "dir containing the file:"
echo "$found_dirs" | awk -F'/' '{print $1}' | sort | uniq -c

echo "dir does not contain the file:"
echo "$not_found_dirs" | awk -F'/' '{print $1"/"$2"/"$3}' | sort | uniq -c
