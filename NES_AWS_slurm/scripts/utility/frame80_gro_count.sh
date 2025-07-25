results=$(find edge_lig_*_lig_*/state*/run*/pre/ -name 'frame80.gro')

echo "$results" | awk -F'/' '{print $1}' | sort | uniq -c
