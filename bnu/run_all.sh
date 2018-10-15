#!/bin/bash
#
set -e

for stage in $( echo A??_*.?90 ); do
    echo "============================== Running $stage"
    entgvsd -d $stage
done

