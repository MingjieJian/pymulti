FC=gfortran
LINK=gfortran
FOPT=-Ofast
OPT=$(FOPT)
#DEBUG=
#SHELL := /bin/zsh
LFLAGS=$(OPT) $(DEBUG)
FFLAGS=$(OPT) $(DEBUG)
.f.o:
	$(FC) $(FFLAGS) $*.f -c

# link lists for various versions
# MULvvops where vv is version (22=version 2.2, 23=version 2.3)
#                o is operator type (G=global, L=local, B=banded)
#                p is opacity package (U=Uppsala is only version supported)
#                s is scratch treatment (S=scratch files, C=scratch in core)
# special versions:
# _3D  reads atmosphere from multi3d format and solves for each column as 1D
# ops=DWG  multi-multi run with global operator
#     DWL  multi-multi run with local opeartor
#

MUL23GUS = mul23_g.o mul23_bconst.o mul23_subg.o mul23_sub.o mul23_opacu.o \
           mul23_writeidl.o transp_23.o

MUL23LUS = mul23_l.o mul23_bconst.o mul23_subl.o mul23_sub.o mul23_opacu.o \
           mul23_writeidl.o transp_23.o

MUL23GUC = mul23_g.o mul23_bconst.o mul23_subg.o mul23_sub_noscratch.o mul23_opacu.o \
           mul23_writeidl.o transp_23.o

MUL23LUC = mul23_l.o mul23_bconst.o mul23_subl.o mul23_sub_noscratch.o mul23_opacu.o \
           mul23_writeidl.o transp_23.o

MUL23BUS = mul23_g.o mul23_bconst.o mul23_subb.o mul23_sub.o mul23_opacu.o \
           mul23_writeidl.o transp_23.o

MUL23GUC_3D = mul23_g_3d.o mul23_bconst_3d.o mul23_subg_3d.o \
           mul23_sub_noscratch_3d.o mul23_opacu.o \
           mul23_writeidl.o transp_23.o var_3d.o

MUL23LUC_3D = mul23_l_3d.o mul23_bconst_3d.o mul23_subl_3d.o \
           mul23_sub_noscratch_3d.o mul23_opacu.o \
           mul23_writeidl.o transp_23.o var_3d.o

MUL23DWG = multi_23_g_dw2.o mul23_bconst.o mul23_subg.o mul23_sub_dw.o \
           mul23_opacu.o mul23_writeidl.o transp_23.o

MUL23DWL = multi_23_l_dw2.o mul23_bconst.o mul23_subl.o mul23_sub_dw.o \
           mul23_opacu.o mul23_writeidl.o transp_23.o


clean:
	\rm -f *.o *.mod

mul23gus.x: $(MUL23GUS)
	$(LINK) $(MUL23GUS) $(LFLAGS) -o mul23gus.x

mul23lus.x: $(MUL23LUS)
	$(LINK) $(MUL23LUS) $(LFLAGS) -o mul23lus.x

mul23guc.x: $(MUL23GUC)
	$(LINK) $(MUL23GUC) $(LFLAGS) -o mul23guc.x

mul23luc.x: $(MUL23LUC)
	$(LINK) $(MUL23LUC) $(LFLAGS) -o mul23luc.x

mul23bus.x: $(MUL23BUS)
	$(LINK) $(MUL23BUS) $(LFLAGS) -o mul23bus.x

mul23guc_3d.x: $(MUL23GUC_3D)
	$(LINK) $(MUL23GUC_3D) $(LFLAGS) $(XTRALIBS) -o mul23guc_3d.x

mul23luc_3d.x: $(MUL23LUC_3D)
	$(LINK) $(MUL23LUC_3D) $(LFLAGS) $(XTRALIBS) -o mul23luc_3d.x

mpif.h:
	cp mpif.h.nompi mpif.h

mul23guc_3d.sx: mpif.h $(MUL23GUC_3D)  nompi.o
	$(LINK) $(MUL23GUC_3D) $(LFLAGS) nompi.o -o mul23guc_3d.sx
	rm mpif.h

mul23luc_3d.sx: mpif.h $(MUL23LUC_3D)  nompi.o
	$(LINK) $(MUL23LUC_3D) $(LFLAGS) nompi.o -o mul23luc_3d.sx
	rm mpif.h

mul23dwg.x: $(MUL23DWG)
	$(LINK) $(MUL23DWG) $(LFLAGS) -o mul23dwg.x

mul23dwl.x: $(MUL23DWL)
	$(LINK) $(MUL23DWL) $(LFLAGS) -o mul23dwl.x

#dependencies:
 
mul23_bconst.o:PREC PARAM CTRAN CSLINE CCONST CBCNST C3D
mul23_bconst_3d.o:PREC PARAM CTRAN CSLINE CCONST CBCNST C3D
mul23_g.o:PREC PARAM PARAMW CATOM CATMOS CATMO2 CTRAN CSLINE CGAUSI CCONST \
         CINPUT CLGMX CLU COPCL
mul23_g_3d.o:PREC PARAM PARAMW CATOM CATMOS CATMO2 CTRAN CSLINE CGAUSI CCONST \
         CINPUT CLGMX CLU CMPI COPCL C3D var_3d.o
mul23_l.o:PREC PARAM PARAMW CATOM CATMOS CATMO2 CTRAN CSLINE CGAUSI CCONST \
         CINPUT CLGMX CLU COPCL
mul23_l_3d.o:PREC PARAM PARAMW CATOM CATMOS CATMO2 CTRAN CSLINE CGAUSI CCONST \
         CINPUT CLGMX CLU CMPI COPCL C3D var_3d.o
mul23_opacu.o:PREC PARAM CATOM CATMOS CATMO2 CSLINE CCONST CINPUT CLU C3D \
         OSMPAR COSM3 PARAMO CALIN CGENCL CTRAN CGAUSI CLGMX CTAUQQ COPCL \
         COSM1 COSM2
mul23_sub.o:PREC PARAM CATOM CATMOS CATMO2 CSLINE CLU CCONST CINPUT CFIX \
         CTERM CTRAN CGAUSI CLGMX CTAUQQ CCNTRB CATOM23 CIMIN COPCL CBCNST
mul23_sub_dw.o:PREC PARAM CATOM CATMOS CATMO2 CSLINE CLU CCONST CINPUT CFIX \
         CTERM CTRAN CGAUSI CLGMX CTAUQQ CCNTRB CATOM23 CIMIN COPCL CDW \
         CBCNST
mul23_sub_noscratch.o:PREC PARAM CATOM CATMOS CATMO2 CSLINE CLU CCONST CINPUT \
         CFIX CTERM CTRAN CGAUSI CLGMX CTAUQQ CCNTRB CATOM23 CIMIN COPCL \
         CBCNST CXJP
mul23_sub_noscratch_3d.o:PREC PARAM CATOM CATMOS CATMO2 CCONST CSLINE CLU \
         CMPI CINPUT CFIX CTERM CTRAN CGAUSI CLGMX CTAUQQ CCNTRB CATOM23 \
         CIMIN COPCL C3D CBCNST CXJP
mul23_subb.o:PREC PARAM CATOM CATMOS CTRAN CSLINE CGAUSI CINPUT CATOM23 \
         CCONST CLU CIMIN
mul23_subg.o:PREC PARAM CATOM CATMOS CTRAN CSLINE CGAUSI CINPUT CATOM23 \
         CCONST CLU CIMIN
mul23_subg_3d.o:PREC PARAM CATOM CATMOS CTRAN CSLINE CGAUSI CINPUT CATOM23 \
         CCONST CLU CIMIN
mul23_subl.o:PREC PARAM CATOM CATMOS CTRAN CSLINE CGAUSI CINPUT CATOM23 \
         CCONST CLU CIMIN
mul23_subl_3d.o:PREC PARAM CATOM CATMOS CTRAN CSLINE CGAUSI CINPUT CATOM23 \
         CCONST CLU CIMIN
mul23_writeidl.o:PREC PARAM CATOM CATMOS CTRAN CSLINE CGAUSI CINPUT CCONST \
         CLGMX CTAUQQ CCNTRB CLU CATMO2
multi_23_g_dw2.o:PREC PARAM PARAMW CATOM CATMOS CATMO2 CTRAN CSLINE CGAUSI \
         CCONST CINPUT CLU COPCL
multi_23_l_dw2.o:PREC PARAM PARAMW CATOM CATMOS CATMO2 CTRAN CSLINE CGAUSI \
         CCONST CINPUT CLU COPCL
transp_23.o:PREC PARAM CATMOS CTRAN CSLINE CINPUT CGAUSI CLU
