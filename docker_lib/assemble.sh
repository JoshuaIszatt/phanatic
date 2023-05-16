#!/bin/bash
. ~/.bashrc
set -e
SECONDS=0
conda activate assemble

python /assemble/bin/coordinator.py
chmod -R 777 /assemble/output/*

if (($SECONDS > 3600)); then
    let "hours=SECONDS/3600"
    let "minutes=(SECONDS%3600)/60"
    let "seconds=(SECONDS%3600)%60"
    echo " " >>/assemble/output/docker_log.tsv
    echo "[hours:$hours,minutes:$minutes,seconds:$seconds]" >>/assemble/output/docker_log.tsv
elif (($SECONDS > 60)); then
    let "minutes=(SECONDS%3600)/60"
    let "seconds=(SECONDS%3600)%60"
    echo " " >>/assemble/output/docker_log.tsv
    echo "[minutes:$minutes,seconds:$seconds]" >>/assemble/output/docker_log.tsv
else
    echo " " >>/assemble/output/docker_log.tsv
    echo "[seconds:$SECONDS]" >>/assemble/output/docker_log.tsv
fi
