#!/bin/csh
#
# script to run multi_3d
# atom in parameter 1, version in parameter 2
#
set atom=$1
set version=$2

if ($atom == "") then
  echo -n "Give atom designation: "
  set atom=($<)
else
  echo "Atom designation is: $atom"
endif
if ($version == "") then
  echo -n "Give multi version (eg 23guc): "
  set version=($<)
else
  echo "Version is: $version"
endif

\rm -f DUM* INIT JNY JOBLOG NIIT OPC OUT PHI RSTRT2 \
               TIME IDL* DSCAL2
\rm -f ABUND ABSDAT ABSLIN ABSMET ATOM ATOM2 ATMOS DSCALE INPUT INPUT_3D
ln -s ../input/abund        ABUND
ln -s ../input/absdat       ABSDAT
ln -s ../input/absmet       ABSMET
ln -s ../input/atom.$atom   ATOM
if(-e ../input/atom2.$atom) ln -s ../input/atom2.$atom ATOM2
if(-e ../input/abslin.$atom) ln -s ../input/abslin.$atom ABSLIN
ln -s ../input/input.$atom  INPUT
ln -s ../input/input_3d.$atom  INPUT_3D
time ./mul$version"_3d.sx"

