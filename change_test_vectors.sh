#!/bin/bash

# Usage: ./update_test_vectors.sh <d0_value> <other_value> <file>

D0_VAL=$1
OTHER_VAL=$2
FILE=$3

# Set d0 value
sed -i -E "s/^(d0 test vectors:)[[:space:]]*[0-9]+/\1 $D0_VAL/" "$FILE"

# Set all other dx values (excluding d0)
sed -i -E "/^d[1-9][0-9]* test vectors:/s/(:)[[:space:]]*[0-9]+/\1 $OTHER_VAL/" "$FILE"

