#!/bin/bash
#
set -e

for stage in A05_trim_laidoy_1kmx1km.F90 A06_trim_lc_laimax_laimonth.F90; do
    echo "============================== Running $stage"
    xent -d $stage
done

echo "=================== Finished run_all.sh"

