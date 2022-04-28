#!/bin/csh
#
# make double precision version of mul23 source code
# 06-10-27  Mats Carlsson
#
cd ../source_dp
#
# delete previous files
#
 echo -n "Will delete all old mul23 files in directory source_dp, is that OK (y/n)? "
if ($< != "y") then 
  exit
endif

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
#
# copy over include files and change REAL to REAL(8)
#
echo "Copying over include files and changing REAL to REAL(8)"
cp ../source/C* .
cp ../source/PARAM .
cp ../source/PARAMW .
cp ../source/PARAMO .
cp ../source/PREC.dp PREC
cp ../source/OSMPAR .

foreach f (C*)
cat $f | sed -e '/REAL /s/REAL /REAL(8) /' >! a.jou
\mv a.jou $f
end
#
# copy over source files
#
echo "Linking over source files that are identical in sp and dp"
ln -s ../source/mul23_bconst.f .
ln -s ../source/mul23_bconst_3d.f .
ln -s ../source/mul23_g.f .
ln -s ../source/mul23_l.f .
ln -s ../source/mul23_subg.f .
ln -s ../source/mul23_subl.f .
ln -s ../source/mul23_subb.f .
ln -s ../source/mul23_sub.f .
ln -s ../source/mul23_sub_noscratch.f .
ln -s ../source/mul23_g_3d.f .
ln -s ../source/mul23_l_3d.f .
ln -s ../source/mul23_subg_3d.f .
ln -s ../source/mul23_subl_3d.f .
ln -s ../source/mul23_sub_noscratch_3d.f .
ln -s ../source/multi_23_g_dw2.f .
ln -s ../source/multi_23_l_dw2.f .
ln -s ../source/mul23_sub_dw.f .
ln -s ../source/nompi.f .
ln -s ../source/mpif.h.nompi .
ln -s ../source/makefile .
ln -s ../source/transp_23.f .
ln -s ../source/var_3d.f .
#
echo "Copying and changing source files that differ between sp and dp"
echo "(mul23_opacu.f mul23_writeidl.f)"
#
# use mul23_writeidldbl.f as mul23_writeidl.f
# 
ln -s ../source/mul23_writeidldbl.f mul23_writeidl.f
#
# change REAL to REAL(8) in mul23_opacu.f (should be the only
# source file with non-include-file-variables that do not
# conform to the FORTRAN implicit typing)
#
cat ../source/mul23_opacu.f | sed -e '/REAL /s/REAL /REAL(8) /' >! mul23_opacu.f
#
cd ../source



