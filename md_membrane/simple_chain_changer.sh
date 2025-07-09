#!/bin/bash

# Simple PDB chain changer
# Usage: ./simple_chain_changer.sh input.pdb old_chain new_chain [output.pdb]

if [ $# -lt 3 ]; then
    echo "Usage: $0 <input.pdb> <old_chain> <new_chain> [output.pdb]"
    echo "Example: $0 protein.pdb A B protein_modified.pdb"
    echo "Use ' ' (space) for empty chain identifier"
    exit 1
fi

input_file="$1"
old_chain="$2"
new_chain="$3"
output_file="${4:-${input_file%.*}_modified.pdb}"

if [ ! -f "$input_file" ]; then
    echo "Error: Input file '$input_file' not found."
    exit 1
fi

if [ "$old_chain" = " " ] || [ "$old_chain" = "space" ] || [ "$old_chain" = "empty" ]; then
    old_chain=" "
fi

cp "$input_file" "$output_file"

if [ "$old_chain" = " " ]; then
    sed -i "s/^\(.\{21\}\) \(.*\)$/\1$new_chain\2/" "$output_file"
else
    sed -i "s/^\(.\{21\}\)$old_chain\(.*\)$/\1$new_chain\2/" "$output_file"
fi

old_count=$(grep "^ATOM\|^HETATM" "$input_file" | cut -c22 | grep -c "^$old_chain$")
echo "Modified $old_count atoms from chain '$old_chain' to chain '$new_chain'" 