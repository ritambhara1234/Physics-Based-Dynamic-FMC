#!/bin/bash
dir=`pwd`
target=${dir##*/}

../Scripts/save_fdsinfo.sh
echo Building $target
make -j4 VPATH="/gpfs/home/rrd20z/WFDS_3/fds/Source" -f /gpfs/home/rrd20z/WFDS_3/fds/Build/makefile $target
