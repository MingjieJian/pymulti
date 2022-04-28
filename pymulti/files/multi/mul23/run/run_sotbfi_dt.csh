#!/bin/csh
#
# script to run multi for response functions
#
set atom=$1
set atmos=$2
set version=$3
if ($atom == "") then
  echo -n "Give atom designation: "
  set atom=($<)
else
  echo "Atom designation is: $atom"
endif
if ($atmos == "") then
  echo -n "Give atmosphere designation: "
  set atmos=($<)
else
  echo "Atmos designation is: $atmos"
endif
if ($version == "") then
  echo -n "Give multi version (eg 20lus): "
  set version=($<)
else
  echo "Version is: $version"
endif
 
# run for perturbed atmospheres
\mv RSTRT2 RSTRT
ln -sf atmos.$atmos ATMOS
ln -sf ../input/input.sotbfi_nq384_init INPUT  # restart with iterations
ln -sf ../input/atom.sotbfi_init ATOM          # full NLTE
time ./mul$version.x
\mv RSTRT2 RSTRT
ln -sf ../input/input.sotbfi INPUT             # restart without iterations
ln -sf ../input/atom.sotbfi_nq384 ATOM         # smaller atom
time ./mul$version.x
if(-e IDL1)   \mv IDL1   idl1.$atom"_"$atmos




