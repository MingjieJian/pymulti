      SUBROUTINE WIDLCN(KR)
C
C  WRITES INTENSITY CONTRIBUTION FUNCTIONS
C:
C: WIDLCN 95-08-16  NEW ROUTINE: (MATS CARLSSON)
C:        WRITES INTENSITY CONTRIBUTION FUNCTIONS
C:        SAME AS IDL VERSION OF WCNTRB IN VERSION 2.1 AND EARLIER
C:
      INCLUDE 'PREC'
      INCLUDE 'PARAM'
      INCLUDE 'CATOM'
      INCLUDE 'CATMOS'
      INCLUDE 'CTRAN'
      INCLUDE 'CSLINE'
      INCLUDE 'CGAUSI'
      INCLUDE 'CINPUT'
      INCLUDE 'CCONST'
      INCLUDE 'CLGMX'
      INCLUDE 'CTAUQQ'
      INCLUDE 'CCNTRB'
      INCLUDE 'CLU'
C
      SAVE LCNTR
C
      IF(IDLCNT.EQ.0) RETURN
C
      IF(KR.EQ.1) THEN
        CALL MOPEN(LCNTR,'IDLCNT',0,'NEW')
        WRITE(LCNTR) NDEP,NLINE,NRAD,MQ
        WRITE(LCNTR) (NQ(MR),MR=1,NRAD)
      ENDIF
      DO 200 NY=1,NQ(KR)
        WRITE(LCNTR) (REAL(CNTRBI(NY,K)),K=1,NDEP)
        WRITE(LCNTR) (REAL(CNTRBF(NY,K)),K=1,NDEP)
  200 CONTINUE
C
      IF(KR.LE.NLINE) THEN
        DO 300 NY=1,NQ(KR)
          WRITE(LCNTR) (REAL(CNTRBR(NY,K)),K=1,NDEP)
  300   CONTINUE
      ENDIF
C
      IF(KR.EQ.NRAD) THEN
        CALL CLOSE(LCNTR)
      ENDIF
C
      RETURN
      END
C
C***********************************************************************
C
  
      SUBROUTINE WIDL1
C
C  WRITES MOST COMMON-BLOCK VARIABLES TO UNFORMATTED FILE IDL1
C  IT IS ASSUMED THAT NLINE.GT.0
C
C:
C: WIDL1  95-08-16  NEW ROUTINE: (MATS CARLSSON)
C:        WRITES MOST COMMON-BLOCK VARIABLES TO UNFORMATTED FILE IDL1
C:        SAME AS IDL VERSION OF WRAD IN VERSIONS 2.1 AND EARLIER
C:
C:        07-02-15  MODIFICATIONS: (MATS CARLSSON)
C:        SETS NWIDE TO NRAD-NLINE+NUMBER OF WIDE LINES
C:
      INCLUDE 'PREC'
      INCLUDE 'PARAM'
      INCLUDE 'CATOM'
      INCLUDE 'CATMOS'
      INCLUDE 'CATMO2'
      INCLUDE 'CTRAN'
      INCLUDE 'CSLINE'
      INCLUDE 'CGAUSI'
      INCLUDE 'CCONST'
      INCLUDE 'CINPUT'
      INCLUDE 'CLU'
C
      DIMENSION KRSEL(MRAD)
      SAVE ICALL,NSEL,KRSEL
      DATA ICALL/0/
C
      IF(IDL1.EQ.0) RETURN
      ICALL=ICALL+1
      IF(IDL1.EQ.2 .AND. ICALL.EQ.1) THEN
        CALL MOPEN(LKRSEL,'KRSEL',1,'OLD')
        READ(LKRSEL,*) NSEL
        READ(LKRSEL,*) (KRSEL(KR),KR=1,NSEL)
        CALL CLOSE(LKRSEL)
      ENDIF
C
      CALL MOPEN(LIDL1,'IDL1',0,'NEW')
C
      NWIDE=NRAD-NLINE
      DO 100 KR=1,NLINE
        IF(IWIDE(KR)) NWIDE=NWIDE+1
  100 CONTINUE
      IF(IDL1.NE.2) THEN
        IFMT=0
        WRITE(LIDL1) NDEP,NK,NLINE,NWIDE,NRAD,NRFIX,NMU,MQ
      ELSE
        IFMT=2
        WRITE(LIDL1) 0,   0, 0,    0,    0,   0,    0,  0
        WRITE(LIDL1) IFMT
        WRITE(LIDL1) NSEL
        WRITE(LIDL1) (KRSEL(KR),KR=1,NSEL)
        WRITE(LIDL1) NDEP,NK,NLINE,NWIDE,NRAD,NRFIX,NMU,MQ
      ENDIF
      IF(NRAD.GT.0) WRITE(LIDL1) (NQ(MR),MR=1,NRAD)
C
C  COMMON BLOCK CATOM
C
      WRITE(LIDL1) REAL(QNORM)
      WRITE(LIDL1) REAL(ABND),REAL(AWGT)
      WRITE(LIDL1) (REAL(EV(I)),I=1,NK)
      WRITE(LIDL1) (REAL(G(I)),I=1,NK)
      WRITE(LIDL1) (ION(I),I=1,NK)
      WRITE(LIDL1) REAL(HN3C2)
      IF(NRAD.GT.0) THEN
        WRITE(LIDL1) (KTRANS(MR),MR=1,NRAD)
        WRITE(LIDL1) (JRAD(MR),MR=1,NRAD)
        WRITE(LIDL1) (IRAD(MR),MR=1,NRAD)
        WRITE(LIDL1) (REAL(F(MR)),MR=1,NRAD)
        WRITE(LIDL1) (IWIDE(MR),MR=1,NRAD)
        WRITE(LIDL1) (REAL(GA(MR)),MR=1,NRAD)
        WRITE(LIDL1) (REAL(GW(MR)),MR=1,NRAD)
        WRITE(LIDL1) (REAL(GQ(MR)),MR=1,NRAD)
      ENDIF
      IF(IFMT.EQ.0) WRITE(LIDL1) ((KRAD(I,J),I=1,NK),J=1,NK)
      IF(IFMT.EQ.0) WRITE(LIDL1) (REAL(Z(K)),K=1,NDEP)
      IF(NWIDE.GT.0) WRITE(LIDL1) ((REAL(ALFAC(NU,KT)),
     * NU=1,MQ),KT=1,NWIDE)
      WRITE(LIDL1) REAL(HNY4P)
      IF(NRAD.GT.0) WRITE(LIDL1) (REAL(ALAMB(MR)),MR=1,NRAD)
      IF(NLINE.GT.0) WRITE(LIDL1) (REAL(A(MR)),MR=1,NLINE)
      IF(IFMT.EQ.0) WRITE(LIDL1) ((REAL(B(I,J)),I=1,NK),J=1,NK)
      WRITE(LIDL1) (REAL(TOTN(K)),K=1,NDEP)
      IF(NRAD.GT.0 .AND. IFMT.EQ.0) 
     * WRITE(LIDL1) ((REAL(BP(K,MR)),K=1,NDEP),MR=1,NRAD)
      WRITE(LIDL1) ((REAL(NSTAR(I,K)),I=1,NK),K=1,NDEP)
      WRITE(LIDL1) ((REAL(N(I,K)),I=1,NK),K=1,NDEP)
      IF(IFMT.EQ.0) 
     * WRITE(LIDL1) (((REAL(C(I,J,K)),I=1,NK),J=1,NK),K=1,NDEP)
      IF(NRFIX.GT.0) THEN
        WRITE(LIDL1) (JFX(MR),MR=1,NRFIX)
        WRITE(LIDL1) (IFX(MR),MR=1,NRFIX)
        WRITE(LIDL1) (IPHO(MR),MR=1,NRFIX)
        WRITE(LIDL1) (REAL(A0(MR)),MR=1,NRFIX)
        WRITE(LIDL1) (REAL(TRAD(MR)),MR=1,NRFIX)
        WRITE(LIDL1) (ITRAD(MR),MR=1,NRFIX)
      ENDIF
      WRITE(LIDL1) (REAL(DNYD(K)),K=1,NDEP)
      IF(NLINE.GT.0 .AND. IFMT.EQ.0) WRITE(LIDL1) ((REAL(ADAMP(K,MR)),
     * K=1,NDEP),MR=1,NLINE)
      WRITE(LIDL1) (LABEL(I),I=1,NK)
      WRITE(LIDL1) ATOMID
      WRITE(LIDL1) CROUT
C
C  COMMON BLOCK CATMOS
C
      WRITE(LIDL1) REAL(GRAV)
      WRITE(LIDL1) (REAL(CMASS(K)),K=1,NDEP)
      WRITE(LIDL1) (REAL(TEMP(K)),K=1,NDEP)
      WRITE(LIDL1) (REAL(NE(K)),K=1,NDEP)
      WRITE(LIDL1) (REAL(VEL(K)),K=1,NDEP)
      WRITE(LIDL1) (REAL(TAU(K)),K=1,NDEP)
      WRITE(LIDL1) (REAL(XNORM(K)),K=1,NDEP)
      WRITE(LIDL1) (REAL(HEIGHT(K)),K=1,NDEP)
      WRITE(LIDL1) ATMOID,DPID,DPTYPE
C
C  COMMON BLOCK CATMO2
C
      WRITE(LIDL1) (REAL(VTURB(K)),K=1,NDEP)
      IF(IFMT.EQ.0) WRITE(LIDL1) ((REAL(BH(J,K)),J=1,5),K=1,NDEP)
      WRITE(LIDL1) ((REAL(NH(J,K)),J=1,6),K=1,NDEP)
      IF(IFMT.EQ.0) WRITE(LIDL1) (REAL(RHO(K)),K=1,NDEP)
C
C  COMMON BLOCK CSLINE
C
      IF(NRAD.GT.0) THEN
        WRITE(LIDL1) (REAL(QMAX(MR)),MR=1,NRAD)
        WRITE(LIDL1) (REAL(Q0(MR)),MR=1,NRAD)
        WRITE(LIDL1) (IND(MR),MR=1,NRAD)
        WRITE(LIDL1) REAL(DIFF)
        WRITE(LIDL1) ((REAL(Q(NU,MR)),NU=1,MQ),MR=1,NRAD)
        IF(IFMT.EQ.0) WRITE(LIDL1) ((REAL(WQ(NU,MR)),NU=1,MQ),MR=1,NRAD)
      ENDIF          
      WRITE(LIDL1) REAL(WQMU)
      IF(NWIDE.GT.0) WRITE(LIDL1) ((REAL(FRQ(NU,MR)),
     *NU=0,MQ),MR=1,NWIDE)
      IF(NRAD.GT.0) THEN
        IF(IFMT.EQ.0) 
     *   WRITE(LIDL1) ((REAL(WPHI(K,MR)),K=1,NDEP),MR=1,NRAD)
        WRITE(LIDL1) ((REAL(SL(K,MR)),K=1,NDEP),MR=1,NRAD)
      ENDIF
      IF(NLINE.GT.0) THEN
        WRITE(LIDL1) (REAL(WEQLTE(MR)),MR=1,NLINE)
        WRITE(LIDL1) (REAL(WEQ(MR)),MR=1,NLINE)
      ENDIF
      IF(NRAD.GT.0) THEN
        IF(IFMT.EQ.0) THEN
          WRITE(LIDL1) ((REAL(RIJ(K,MR)),K=1,NDEP),MR=1,NRAD)
          WRITE(LIDL1) ((REAL(RJI(K,MR)),K=1,NDEP),MR=1,NRAD)
        ENDIF
        IF(IFMT.NE.2) THEN
          WRITE(LIDL1) ((REAL(FLUX(NU,MR)),NU=0,MQ),MR=1,NRAD)
          WRITE(LIDL1) (((REAL(OUTINT(NU,MU,MR)),NU=0,MQ),
     *     MU=1,NMU),MR=1,NRAD)
        ELSE
          DO 200 KR=1,NSEL
            WRITE(LIDL1) (REAL(FLUX(NU,KRSEL(KR))),NU=0,MQ)
  200     CONTINUE
          DO 300 KR=1,NSEL
            WRITE(LIDL1) ((REAL(OUTINT(NU,MU,KRSEL(KR))),NU=0,MQ),
     *       MU=1,NMU)
  300     CONTINUE
        ENDIF
        IF(IFMT.EQ.0) 
     *   WRITE(LIDL1) ((REAL(COOL(K,MR)),K=1,NDEP),MR=1,NRAD)
      ENDIF
C
C  COMMON BLOCK CGAUSI
C
      WRITE(LIDL1) (REAL(XMU(MU)),MU=1,NMU)
      WRITE(LIDL1) (REAL(WMU(MU)),MU=1,NMU)
C
C  COMMON BLOCK CCONST
C
      WRITE(LIDL1) REAL(EE),REAL(HH),REAL(CC),REAL(BK),REAL(EM),
     * REAL(UU),REAL(HCE),REAL(HC2),REAL(HCK),REAL(EK),REAL(PI)

      CALL CLOSE(LIDL1)
C
      END
C
C***********************************************************************
C
      SUBROUTINE WIDLNY(KR,NY)
C
C  WRITES MOST COMMON-BLOCK VARIABLES TO UNFORMATTED FILE IDLNY
C  REGULATED BY PRINTOUT VARIABLE IDLNY
C:
C: WIDLNY 95-08-16  NEW ROUTINE: (MATS CARLSSON)
C:        WRITES A FILE IDLNY WITH DATA SUITABLE FOR INPUT TO IDL
C:        SAME AS IDL VERSION OF WTEST IN VERSION 2.1 AND EARLIER
C:
      INCLUDE 'PREC'
      INCLUDE 'PARAM'
      INCLUDE 'CATOM'
      INCLUDE 'CATMOS'
      INCLUDE 'CATMO2'
      INCLUDE 'CTRAN'
      INCLUDE 'CSLINE'
      INCLUDE 'CGAUSI'
      INCLUDE 'CCONST'
      INCLUDE 'CINPUT'
      INCLUDE 'CLU'
C
      SAVE LIDLNY,IREC
C
      IF(IDLNY.EQ.0) RETURN
C
      IF(KR.EQ.1 .AND. NY.EQ.1) THEN
        NREC=NDEP*11
        CALL MOPEN(LIDLNY,'IDLNY',(NREC+1)/2,'NEW')
        IREC=0
      ENDIF
C
C  NY DEPENDENT VARIABLES FROM COMMON BLOCK CTRAN AND CSLINE
C
      IREC=IREC+1
      WRITE(LIDLNY,REC=IREC) (REAL(PMS(K)),K=1,NDEP),
     * (REAL(IPLUS(K)),K=1,NDEP),(REAL(IMINUS(K)),K=1,NDEP),
     * (REAL(P(K)),K=1,NDEP),(REAL(S(K)),K=1,NDEP),
     * (REAL(TAUQ(K)),K=1,NDEP),(REAL(DTAUQ(K)),K=1,NDEP),
     * (REAL(XCONT(K)),K=1,NDEP),(REAL(SC(K)),K=1,NDEP),
     * (REAL(SCAT(K)),K=1,NDEP),(REAL(X(K)),K=1,NDEP)

      IF(KR.EQ.NRAD .AND. NY.EQ.NQ(KR)) THEN
        CALL CLOSE(LIDLNY)
      ENDIF
C
      END
C
C***********************************************************************
C
