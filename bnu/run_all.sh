#!/bin/bash
#
set -e

for stage in $( echo A??_*.?90 ); do
    echo "============================== Running $stage"
    entgvsd $stage
done

echo "=================== Finished run_all.sh"

