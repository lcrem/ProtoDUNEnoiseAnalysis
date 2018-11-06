#!/bin/bash

if [ -z "$1" ]
  then
    echo "./createFileList.sh [RUN] [SUBRUN]"
    echo "if no subrun is provided it will create a file list for the whole run"
    exit 1
fi

RUN=$1
SUBRUN=$(printf "%04d" $2)
outdir=/dune/data/users/lcremone/ProtoDuneNoise/filelists/
outname=$outdir/np04_raw_run00${RUN}_${SUBRUN}.list

if [ -z "$2" ]
    then
    SUBRUN='%'
    outname=$outdir/np04_raw_run00${RUN}_allsubruns.list
fi

filename=np04_raw_run00${RUN}_${SUBRUN}_dl%_decoder_%.root

echo $filename

sortedlist=`samweb list-files "file_name like $filename" | sort -V`

echo $sortedlist
prefix='enstore:'
for file in `echo $sortedlist`;
do 
    filepath=`samweb locate-file $file | grep pnfs`; 
    filepath=${filepath#"$prefix"}
    filepath=${filepath%(*}
    ls $filepath/*root
done > $outname

echo "Saved list in $outname"
