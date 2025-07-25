#!/bin/bash

files=$(find lig_*/run1*/cal03*/ -name 'md.log')

for file in $files; do
  last_line=$(tail -n 2 "$file")
  if [[ $last_line != Finished* ]]; then
    echo "$file"
  fi
done
