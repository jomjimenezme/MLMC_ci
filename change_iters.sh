#!/bin/bash

# run with: ./update_iters.sh <new_number> <file>

NEW_VAL=$1
FILE=$2

# Update all lines with "trace max iters:" or "trace min iters:"
sed -i -E "s/(trace (max|min) iters:)[[:space:]]*[0-9]+/\1 $NEW_VAL/" "$FILE"

