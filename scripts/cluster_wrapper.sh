#!/usr/bin/env bash
# Should copy input to temporary location and handle removing it automatically.
# Useful because we want to temporarily put the file in /dev/shm but make sure 
# we're not leaving it there taking up memory from other cluster users.

function cleanup {
  rm -rf $2
}

trap cleanup EXIT
trap cleanup SIGINT

echo Creating directory `dirname $2`...
mkdir -p `dirname $2`
echo Copying $1 to $2...
cp $1 $2
echo Running create_uhgp_db.py...
scripts/create_uhgp_db.py -i $2 -o $3

