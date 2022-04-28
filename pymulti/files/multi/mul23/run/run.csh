#!/bin/csh
#
# script to run multi
# atom in parameter 1, atmosphere in parameter 2
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
 
\rm -f DUM* INIT JNY JOBLOG NIIT OPC OUT PHI RSTRT2 \
               TIME IDL* DSCAL2
\rm -f ABUND ABSDAT ABSLIN ABSMET ATOM ATOM2 ATMOS DSCALE INPUT
ln -s ../input/abund        ABUND
ln -s ../input/absdat       ABSDAT
ln -s ../input/absmet       ABSMET
ln -s ../input/atom.$atom   ATOM
if(-e ../input/atom2.$atom) ln -s ../input/atom2.$atom ATOM2
if(-e ../input/abslin.$atom) ln -s ../input/abslin.$atom ABSLIN
ln -s ../input/atmos.$atmos ATMOS
if(-e ../input/dscale.$atom"_"$atmos) then
  ln -s ../input/dscale.$atom"_"$atmos DSCALE
else
  ln -s ../input/dscale.$atmos DSCALE
endif
ln -s ../input/input.$atom  INPUT
time ./mul$version.x
if(-e IDL1)   \mv IDL1   idl1.$atom"_"$atmos
if(-e IDLCNT) \mv IDLCNT idlcnt.$atom"_"$atmos
if(-e IDLOPC) \mv IDLOPC idlopc.$atom"_"$atmos
if(-e IDLNY)  \mv JNY    jny.$atom"_"$atmos
if(-e IDLNY)  \mv IDLNY  idlny.$atom"_"$atmos
if(-e RSTRT2) \mv RSTRT2 rstrt.$atom"_"$atmos
if(-e DUMC)   \mv DUMC   dumc.$atom"_"$atmos
if(-e NIIT)   \mv NIIT   niit.$atom"_"$atmos
if(-e OUT)    \mv OUT    out.$atom"_"$atmos

