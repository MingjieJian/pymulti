#!/bin/sh
#
# make double precision version of mul23 source code
# 22-04-25  Transfer to bash format (Mingjie Jian)
# 06-10-27  Mats Carlsson
#
cd ../source_dp
#
# delete all files in source_dp
#

force=0
while getopts 'f' FLAG; do
  case "$FLAG" in
    f)
      echo "Force mode"
      force=1
      ;;
    ?)
      echo "script usage: $(basename \$0) [-f]" >&2
      exit 1
      ;;
  esac
done

while [ $force -eq 0 ]; do
    read -p "Will delete all old mul23 files in directory source_dp, is that OK (y/n)? " yn
    case $yn in
        [Yy]* ) break;;
        [Nn]* ) echo 'Exit'; exit;;
        * ) echo "Please answer yes or no.";;
    esac
done

echo "Deleting all old mul23 files in directory source_dp."
\rm -f ../source_dp/mul23*.f
\rm -f ../source_dp/multi*dw2.f
\rm -f ../source_dp/transp_23.f
\rm -f ../source_dp/var_3d.f
\rm -f ../source_dp/nompi.f
\rm -f ../source_dp/C*
\rm -f ../source_dp/OSMPAR
\rm -f ../source_dp/PREC
\rm -f ../source_dp/PARAM*
\rm -f ../source_dp/makefile

# copy over include files and change REAL to REAL(8)

echo "Copying over include files and changing REAL to REAL(8)"
cp ../source/C* .
cp ../source/PARAM .
cp ../source/PARAMW .
cp ../source/PARAMO .
cp ../source/PREC.dp PREC
cp ../source/OSMPAR .

for f in C*
do
  cat $f | sed -e '/REAL /s/REAL /REAL(8) /' > a.jou
  mv a.jou $f
done

# copy over source files
echo "Linking over source files that are identical in sp and dp"
ln -sf ../source/mul23_bconst.f .
ln -sf ../source/mul23_bconst_3d.f .
ln -sf ../source/mul23_g.f .
ln -sf ../source/mul23_l.f .
ln -sf ../source/mul23_subg.f .
ln -sf ../source/mul23_subl.f .
ln -sf ../source/mul23_subb.f .
ln -sf ../source/mul23_sub.f .
ln -sf ../source/mul23_sub_noscratch.f .
ln -sf ../source/mul23_g_3d.f .
ln -sf ../source/mul23_l_3d.f .
ln -sf ../source/mul23_subg_3d.f .
ln -sf ../source/mul23_subl_3d.f .
ln -sf ../source/mul23_sub_noscratch_3d.f .
ln -sf ../source/multi_23_g_dw2.f .
ln -sf ../source/multi_23_l_dw2.f .
ln -sf ../source/mul23_sub_dw.f .
ln -sf ../source/nompi.f .
ln -sf ../source/mpif.h.nompi .
ln -sf ../source/makefile .
ln -sf ../source/transp_23.f .
ln -sf ../source/var_3d.f .
#
echo "Copying and changing source files that differ between sp and dp"
echo "(mul23_opacu.f mul23_writeidl.f)"
#
# use mul23_writeidldbl.f as mul23_writeidl.f
# 
ln -sf ../source/mul23_writeidldbl.f mul23_writeidl.f
#
# change REAL to REAL(8) in mul23_opacu.f (should be the only
# source file with non-include-file-variables that do not
# conform to the FORTRAN implicit typing)
#
cat ../source/mul23_opacu.f | sed -e '/REAL /s/REAL /REAL(8) /' > mul23_opacu.f
#
cd ../source