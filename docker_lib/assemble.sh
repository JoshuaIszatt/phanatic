#!/bin/bash
. ~/.bashrc
set -e
SECONDS=0

# Downloading checkV database if needed
if [ -d "/assemble/database/checkv-db-v1.5" ]; then
    echo "Found checkV database"
    export CHECKVDB=/assemble/database/$(basename /assemble/database/checkv*)
else
    echo "Could not find checkV database, downloading one now"
    checkv download_database /assemble/database
    echo "$(basename /assemble/database/che*) exists"
    export CHECKVDB=/assemble/database/$(basename /assemble/database/che*)
fi

# Running coordinator
echo "Running coordinator script"
python /assemble/bin/coordinator.py
chmod -R 777 /assemble/output/*

# Running mapping file
conda activate map
echo "Running mapping script"
python /assemble/bin/host_mapping.py
chmod -R 777 /assemble/output/*

# Running re-assembly
conda deactivate
echo "Running final assembly script"
python /assemble/bin/host_mapping.py
chmod -R 777 /assemble/output/*


if (($SECONDS > 3600)); then
    let "hours=SECONDS/3600"
    let "minutes=(SECONDS%3600)/60"
    let "seconds=(SECONDS%3600)%60"
    echo " " >>/assemble/output/phanatic_log.tsv
    echo "[hours:$hours,minutes:$minutes,seconds:$seconds]" >>/assemble/output/phanatic_log.tsv
elif (($SECONDS > 60)); then
    let "minutes=(SECONDS%3600)/60"
    let "seconds=(SECONDS%3600)%60"
    echo " " >>/assemble/output/phanatic_log.tsv
    echo "[minutes:$minutes,seconds:$seconds]" >>/assemble/output/phanatic_log.tsv
else
    echo " " >>/assemble/output/phanatic_log.tsv
    echo "[seconds:$SECONDS]" >>/assemble/output/phanatic_log.tsv
fi
