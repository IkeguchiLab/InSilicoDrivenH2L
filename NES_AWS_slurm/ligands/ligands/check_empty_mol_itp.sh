#!/bin/bash

for dir in */; do
    if [[ -f "${dir}/MOL.itp" ]]; then
        if [[ ! -s "${dir}/MOL.itp" ]] ; then
            echo "$dir"
	fi
    fi
done
