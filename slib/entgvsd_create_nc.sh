#!/bin/sh -f
#
# sh ../slib/create_nc.sh ./test_out /home2/rpfische/entgvsd0/discover/Vegcover_1km/BNU/lc_lai_ent ./templates checksum EntMM29lc_lai_for_1kmx1km

set -e     # Quit on error
BASEDIR=$(dirname "$0")

OROOT=$1
OROOT_ORIG=$2
TROOT=$3        # Template directory
DIR=$4          # Directory of the file (after the root)
LEAF=$5         # Leaf name (without the .nc)

TMP=./tmp

mkdir -p $OROOT/$DIR

# Create the template if not already there
if [ ! -f $TROOT/$DIR/$LEAF.cdl ]; then
    mkdir -p $TROOT/$DIR

    # Unzip the original .nc file
    orig_nc=$OROOT_ORIG/$DIR/$LEAF.nc
    if [ ! -f $orig_nc ]; then
        mkdir -p $TMP
        if [ ! -f $orig_nc.gz ]; then
            echo "Unzipping $orig_nc.gz.gz"
            zcat $orig_nc.gz.gz >$TMP/x.nc.gz
            echo "Unzipping $orig_nc.gz"
            zcat $TMP/x.nc.gz | head -c 2000000 >$TMP/x.nc
        else
            echo "Unzipping $orig_nc.gz"
            zcat $orig_nc.gz >$TMP/x.nc
        fi
        orig_nc=$TMP/x.nc
    fi


    # echo ncdump -h $OROOT_ORIG/$DIR/$LEAF.nc >$TROOT/$DIR/$LEAF.cdl
    ncdump -h $orig_nc >$TROOT/$DIR/$LEAF.cdl
fi

# Create the output file from the template
# Use our ncgen replacement that compresses every variable
# ncgen -o $OROOT/$DIR/$LEAF.nc $TROOT/$DIR/$LEAF.cdl
python3 $BASEDIR/cdlparser.py $TROOT/$DIR/$LEAF.cdl $OROOT/$DIR/$LEAF.nc

