#!/bin/bash
#
#############################################################################
#
# name: parallel.sh
# date: 17 Sep 14
# auth: Zach Hartiwg
# mail: hartwig@psfc.mit.edu
#
# desc: Script to build the parallel version of ZKBrem. The script
#       (at present) because getting .PHONY targets to work within the
#       Geant4 GNU make system has proven difficult; hence, the script.
# 
#       To build ZKBrem with MPI parallelization: "$ ./parallel"
#       To clean ZKBrem parallel build directory: "$ ./parallel clean"
#
#############################################################################

PAR_TARGET=ZKBrem_MPI

BINARY=$G4WORKDIR/bin/$G4SYSTEM/$PAR_TARGET
BUILDDIR=$G4WORKDIR/tmp/$G4SYSTEM/$PAR_TARGET

# Remove the old binary
if [ -a $BINARY ]
    then
    rm $BINARY
fi

# Command-line capability to clean the build directory
if [ "$1" == 'clean' ];
then
    echo -e "\nCleaning up $PAR_TARGET build directory ...\n"
    rm -rf $BUILDDIR
    exit
else
    echo -e "\nBuilding the parallel version of ZKBrem ...\n"
fi

# Create the parallel implementation file
cp ZKBrem.cc $PAR_TARGET.cc

# If the parallel binary exists, remove it
if [ -a $BINARY ]; then rm $BINARY; fi

# Get the number of available effective cores 

if [[ "$OSTYPE" == "linux-gnu" ]]; then
        # linux
    CORES=$(grep -c ^processor /proc/cpuinfo)
elif [[ "$OSTYPE" == "darwin"* ]]; then
        # Mac OSX
    CORES=$(sysctl -a | grep core_count | sed 's/machdep.cpu.core_count: //')
elif [[ "$OSTYPE" == "freebsd"* ]]; then
        # ...
    CORES=$(sysctl -a | grep core_count | sed 's/machdep.cpu.core_count: //')
else
     echo 'Unknown OS';
     CORES=-1;
fi

# Build the parallel binary from source
make ARCH=parallel -j$CORES

# Test to ensure that the parallel binary was built successfully
if [ -a $BINARY ]
    then
    echo -e '\n************************************************************'
    echo -e '   ZKBrem was successfully built in parallel for '$CORES' cores!'
    echo -e '                MPI parallelization enabled!'
    echo -e '************************************************************\n'
fi

# Remove the parallel implementation file
rm -f $PAR_TARGET.cc
