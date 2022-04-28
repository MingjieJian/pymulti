#!/bin/bash
#
# script to run multi
# atom in parameter 1, atmosphere in parameter 2
#

atom=$1
atmos=$2
version=$3

MULTI_PATH=~/.pymulti/files/multi/mul23

# Usage text
# TO be added.

if [ "$atom" = "" ]; then
  echo -n "Give atom designation: "
  read atom
else
  echo "Atom designation is: $atom"
fi
if [ "$atmos" = "" ]; then
  echo -n "Give atmosphere designation: "
  read atmos
else
  echo "Atmos designation is: $atmos"
fi
if [ "$version" = "" ]; then
  echo -n "Give multi version (eg 20lus): "
  read version
else
  echo "Version is: $version"
fi

# Remove all the OUTPUT files 
rm -f DUM* INIT JNY JOBLOG NIIT OPC OUT PHI RSTRT2 \
               TIME IDL* DSCAL2
rm -f ABUND ABSDAT ABSLIN ABSMET ATOM ATOM2 ATMOS DSCALE INPUT

# Mandatory input
ln -s abund        ABUND
ln -s $MULTI_PATH/input/absdat       ABSDAT
ln -s ./atom.$atom   ATOM
ln -s ./atmos.$atmos ATMOS
ln -s ./input.$atom  INPUT
# if [ -e $MULTI_PATH/input/dscale.$atom"_"$atmos ]; then
#   ln -s $MULTI_PATH/input/dscale.$atom"_"$atmos DSCALE
# else
#   ln -s $MULTI_PATH/input/dscale.$atmos DSCALE
# fi
ln -s ./dscale.$atom"_"$atmos DSCALE

# Optional input
ln -s $MULTI_PATH/input/absmet       ABSMET
if [ -f $MULTI_PATH/input/atom2.$atom ]; then ln -s $MULTI_PATH/input/atom2.$atom ATOM2; fi
if [ -f $MULTI_PATH/input/abslin.$atom ]; then ln -s $MULTI_PATH/input/abslin.$atom ABSLIN; fi

# Run MULTI
ln -s $MULTI_PATH/source_dp/mul$version.x .
time ./mul$version.x
rm ./mul$version.x

# Manage the output files
if [ -f IDL1 ];   then mv IDL1   idl1.$atom"_"$atmos; fi
if [ -f IDLCNT ]; then mv IDLCNT idlcnt.$atom"_"$atmos; fi
if [ -f IDLOPC ]; then mv IDLOPC idlopc.$atom"_"$atmos; fi
if [ -f IDLNY ]; then mv JNY    jny.$atom"_"$atmos; fi
if [ -f IDLNY ]; then mv IDLNY  idlny.$atom"_"$atmos; fi
if [ -f RSTRT2 ]; then mv RSTRT2 rstrt.$atom"_"$atmos; fi
if [ -f DUMC ]; then mv DUMC   dumc.$atom"_"$atmos; fi
if [ -f NIIT ]; then mv NIIT   niit.$atom"_"$atmos; fi
if [ -f OUT ]; then mv OUT    out.$atom"_"$atmos; fi

