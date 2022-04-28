      SUBROUTINE BMAT(W,NK1,NDEP1,J,I)
C
C  FILLS THE DELTA LAMBDA ELEMENTS IN THE PRECONDITIONED MATRIX W
C  IN THE EQUATION
C
C    W*DN/N = E
C
C  DN IS THE APPROXIMATE CORRECTION, REQUIRED TO MAKE THE POPULATIONS
C  SATISFY THE EQUATIONS OF STATISTICAL EQUILIBRIUM.
C
C:
C: BMAT   90-06-08  MODIFICATIONS: (MARTIN J STIFT)
C:                  PURELY DIAGONAL LAMBDA OPERATOR
C:
C:        90-08-14  MODIFICATIONS: (MATS CARLSSON)
C:        CHANGED TO CONTAIN DECLARATION OF W 
C:        W PASSED IN ARGUMENT LIST INSTEAD OF IN COMMON TO
C:        AVOID NK=MK1 AND NDEP=MDEP1 REQUIREMENT
C:
C:        90-09-11  MODIFICATIONS: (MATS CARLSSON)
C:        CHANGED TO DO THE INTEGRATION OF ALI INSIDE LOOP
C:
C:        91-02-01  MODIFICATIONS: (MATS CARLSSON)
C:        ERROR IN LOCAL OPERATOR CORRECTED
C:
C:        97-06-27  MODIFICATIONS: (MATS CARLSSON)
C:        CHANGED TO USE RYBICKI-HUMMER OPERATOR
C:
      INCLUDE 'PREC'
      INCLUDE 'PARAM'
      INCLUDE 'CATOM'
      INCLUDE 'CATMOS'
      INCLUDE 'CTRAN'
      INCLUDE 'CSLINE'
      INCLUDE 'CGAUSI'
      INCLUDE 'CINPUT'
C
      COMMON /CRYBH/ SP21,SP31,SP1N,SP2N
C
      DIMENSION XTEMP(MDEP),D1(MDEP),E1(MDEP),ALOP(MDEP)
      DIMENSION W(NK1,NK1,NDEP1)
C
      WC = -0.5 * WQMU / HNY4P
C
C  FIND COEFFICIENTS FOR RYBICKI HUMMER OPERATOR
C
      D1(1)=SP31/(SP31-SP21)
      E1(NDEP)=SP1N/(SP1N-SP2N)
      DO 410 K=2,NDEP
        XTEMP(K)=1./(1.+C1(K)+A1(K)*(1-D1(K-1)))
        D1(K)=C1(K)*XTEMP(K)
  410 CONTINUE
      DO 420 K=NDEP-1,1,-1
        E1(K)=A1(K)/(1.+A1(K)+C1(K)*(1.-E1(K+1)))
  420 CONTINUE
C
      ALOP(1)=1./(1.-D1(1)*E1(2))/(SP21-SP31)
      ALOP(NDEP)=1./(SP2N-SP1N*(1.-D1(NDEP-1)))
      DO 440 K=2,NDEP-1
        ALOP(K)=XTEMP(K)/(1.-D1(K)*E1(K+1))
  440 CONTINUE
C
      DO 500 K=1,NDEP-1
        DW = WC * ALFA(K) * Z(K) * ALFA(K) / (X(K) * XNORM(K))
        ALIDW = DW*(ALOP(K)-1.)
        W(I,J,K) = W(I,J,K) + GIJ(K)*ALIDW*
     *   (IPLUS(K) + IMINUS(K) + 2.*HN3C2)
        W(J,I,K) = W(J,I,K) + ALIDW*(IPLUS(K) + IMINUS(K))
  500 CONTINUE
C
C  BOTTOM BOUNDARY, SET J=S  ==  FILL NO ELEMENTS
C
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE ITER(W,E,NK1,NDEP1)
C
C  ADMINISTERS ITERATION.
C  FILLS THE MATRICES AND SOLVES THE SYSTEM OF EQUATIONS.
C:
C: ITER   87-12-10  MODIFICATIONS: (MATS CARLSSON)
C:        IF ISUM=0 ISUM IS SET TO THE LEVEL WITH THE LARGEST
C:        POPULATION AT EACH DEPTH-POINT. THIS MODIFICATION ALSO
C:        AFFECTS VALCHK AND STATEQ
C:
C:        88-04-21  MODIFICATIONS: (MATS CARLSSON)
C:        VERSION THAT IMPLEMENTS COLLISIONAL-RADIATIVE SWITCHING,
C:        HUMMER,D.G., VOELS,S.A.: 1988, ASTRON. ASTROPHYS. 192,279
C:        SWITCHING IS REGULATED BY ICONV.
C:        ICONV=2 GIVES INTERACTIVE MODE
C:              3 GIVES AUTOMATIC MODE
C:
C:        89-03-18  MODIFICATIONS: (MATS CARLSSON)
C:        EMAX PRINTOUT CHANGED TO INCLUDE SIGN AND
C:        LEVEL/TRANSITION AND DEPTH 
C:
C:        89-06-06  MODIFICATIONS: (MATS CARLSSON)
C:        HSEINT ONLY CALLED IF SWITCH=1.0
C:
C:        89-06-07  MODIFICATIONS: (MATS CARLSSON)
C:        WHEN SWITCHING WAS USED COLLISIONAL RATES WERE TAKEN
C:        FROM FIRST CALCULATED VALUES THROUGH VARIABLE COL READ
C:        FROM DUMC FILE. THIS GAVE ERRORS FOR HYDROGEN WHEN HSE
C:        INTEGRATIONS WERE PERFORMED SINCE THE UPDATED COLLISIONAL 
C:        RATES WERE OVERWRITTEN WITH THE OLD VALUES. FIXED SO
C:        THAT COL IS ONLY USED WHEN SWITCH.GT.1.0 AND THE FIRST
C:        TIME SWITCH=1.0. SECOND CASE CONTROLLED BY VARIABLE OLDSW
C:
C:        90-01-11  MODIFICATIONS: (MATS CARLSSON)
C:        X(K) MAY BE 0 IF STIMULATED EMISSION EXCEEDS ABSORPTION
C:        X(K) IS SET TO SIGN(1,X(K))*MAX(ABS(X(K)),1.E-4*XCONT(K))
C:
C:        90-06-05  MODIFICATIONS: (MARTIN J. STIFT)
C:        IMPLEMENTED OPTIMUM DIAGONAL LAMBDA OPERATOR ACCORDING 
C:        TO OLSON, AUER AND BUCHLER JQSRT 35,431 (1986) (SEE ALSO
C:        PULS AND HERRERO A&A 204, 219 (1988))
C:        WITH NG ACCELERATION. 
C:        CODING MADE COMPATIBLE WITH OLD MULTI FORMAT BY MC.
C:
C:        90-08-14  MODIFICATIONS: (MATS CARLSSON)
C:        CHANGED TO CONTAIN DECLARATION OF W AND E
C:        THESE PASSED IN ARGUMENT LIST INSTEAD OF IN COMMON TO
C:        AVOID NK=MK1 AND NDEP=MDEP1 REQUIREMENT
C:
C:        90-10-11  MODIFICATIONS: (MATS CARLSSON)
C:        ALI INTEGRATION MOVED TO BMAT
C:
C:        90-10-18  MODIFICATIONS: (MATS CARLSSON)
C:        NG ACCELERATION ALSO INCLUDED IN HSE INTEGRATION.
C:        NG ACCELERATION IS SWITCHED ON WHEN HSE STEPS ARE
C:        PERFORMED FOR EACH ITERATION (EACH HSE STEP GIVES RISE
C:        TO EMAX.LT.ELIM1). THIS IS ACHIEVED THROUGH ITERATION
C:        COUNT VARIABLE ITHSE
C:
C:        92-06-05  MODIFICATIONS: (MATS CARLSSON)
C:        TEST OF KREJ AND KEJ INCLUDED TO AVOID ACCESSING
C:        EJ(0,0) WHEN NRAD=0
C:
C:        92-10-08  MODIFICATIONS: (MATS CARLSSON)
C:        INCIDENT RADIATION FIELD IS SET IF ITRAN.GE.10
C:
C:        95-08-16  MODIFICATIONS: (MATS CARLSSON)
C:        BOTH OLD PRINTOUT ROUTINES (TO FILE OUT) AND IDL PRINTOUT
C:        ROUTINES ARE CALLED
C:
C:        95-08-17  MODIFICATIONS: (MATS CARLSSON)
C:        INCOMING RADIATION SWITCHED ON BY INCRAD.NE.0 INSTEAD OF
C:        ITRAN=10-14
C:
C:        95-08-21  MODIFICATIONS: (MATS CARLSSON)
C:        WARNING MESSAGES REGULATED BY IWARN
C:
C:        95-08-22  MODIFICATIONS: (MATS CARLSSON)
C:        COLLISIONAL RADIATIVE SWITCHING NOW SWITCHED ON WITH ICRSW
C:        ICRSW=-2 GIVES AUTOMATIC MODE
C:        ICRSW=-1 GIVES INTERACTIVE MODE
C:        ICRSW.GT.0 GIVES DECREASE OF
C:        THE SWITCHING PARAMETER WITH ICRSW STEPS PER DECADE
C:        IN THE INTERACTIVE MODE, A REDO OF THE PREVIOUS ITERATION IS ASKED
C:        FOR WITH A NEGATIVE VALUE OF SWITCH. ZERO CAUSES SWITCH TO 
C:        AUTOMATIC MODE
C:
C:        99-10-28  MODIFICATIONS: (MATS CARLSSON)
C:        C(I,J,K)=0 IS TESTED FOR TO AVOID DIVISION BY ZERO IN SETTING
C:        FIRST SWITCHING VALUE
C:
C:        03-08-12  MODIFICATIONS: (MATS CARLSSON)
C:        CHANGED FORMAT FOR WARNING PRINTOUT ABOUT NEGATIVE OPACITIES
C:
C:        05-01-28  MODIFICATIONS: (MATS CARLSSON)
C:        AVOID USING DUMC FILE IF NRFIX=0
C:
C:        06-10-27  MODIFICATIONS: (MATS CARLSSON)
C:        AVOID USING DUMC FILE, COL AND FIX STORED
C:
C:        08-03-25  MODIFICATIONS: (MATS CARLSSON)
C:        ONLY COUNT ITERATIONS WHEN CRSW=1
C:
C:        11-06-03  MODIFICATIONS: (MATS CARLSSON)
C:        CHECK FOR NAN, MODIFIED SWITCH BEHAVIOUR
C:
      INCLUDE 'PREC'
      INCLUDE 'PARAM'
      INCLUDE 'CATOM'
      INCLUDE 'CATOM23'
      INCLUDE 'CATMOS'
      INCLUDE 'CTRAN'
      INCLUDE 'CSLINE'
      INCLUDE 'CGAUSI'
      INCLUDE 'CCONST'
      INCLUDE 'CINPUT'
      INCLUDE 'CLU'
      INCLUDE 'CIMIN'
      COMMON /CNGACC/ NGACC
      LOGICAL NGACC
C
      DIMENSION W(NK1,NK1,NDEP1),E(NK1,NDEP1)
      CHARACTER*16 ETEXT
      DIMENSION DEIK(MDEP)
      DIMENSION AN(MK,MDEP)
      LOGICAL NEWMAT,CONT,HSE
      LOGICAL LISUM,NGON,NGHSE,LWARN
      INTEGER ISUMA(MDEP)
C
      CALL MCTIME('SUM IT',0,0,3)
      ONE=1.0
      HSE=ATOMID.EQ.'H   ' .AND. IHSE.NE.0
      EMAX=10.
      SWITCH=1.
      OLDSW=SWITCH
      EREDO=0.0
      FS0=9.
      DFSMAX=1.3
      NGON=.FALSE.
      NGHSE=.FALSE.
      ITHSE=0
C
C  SWITCHING INITIALIZATION
C
      IF(ICRSW.NE.0) THEN
C
C  CALL TRPT TO GET RIJ AND RJI
C
        CALL TRPT
C
C  FIND INITIAL SWITCHING PARAMETER VALUE
C
        DO 150 KR=1,NLINE
          I=IRAD(KR)
          J=JRAD(KR)
          DO 100 K=1,NDEP
            IF(COL(I,J,K).GT.0.0 .AND. COL(I,J,K).GT.0.0) THEN
              SWITCH=MAX(SWITCH,RIJ(K,KR)/COL(I,J,K))
              SWITCH=MAX(SWITCH,RJI(K,KR)/COL(J,I,K))
            ENDIF
  100     CONTINUE
  150   CONTINUE
        IF(ICRSW.EQ.-1) THEN
          WRITE(*,160) ' MAX(Rij/Cij) IS:',SWITCH
  160     FORMAT(A,1P,E10.2)
          WRITE(*,*) ' GIVE INITIAL SWITCH VALUE TO BE USED'
          READ(*,*) SWITCH
        ELSE
          SWFAC=10.**(-1./ICRSW)
          SWITCH=SWITCH*10.
        ENDIF
      ENDIF
C
      IF(ICRSW.EQ.0) THEN
        ITMAX1=ITMAX
      ELSE
        ITMAX1=1000
      ENDIF
      IT2=1
      DO 9000 IT=1,ITMAX1
        ITHSE=ITHSE+1
  170   CONTINUE
        IF(SWITCH.GT.1.0 .OR. OLDSW.NE.1.0) THEN
C
C  CALCULATE SWITCH SCALED C
C
          DO 200 K=1,NDEP
            DO 190 J=1,NK
              DO 180 I=1,NK
                C(I,J,K)=COL(I,J,K)*SWITCH + FIX(I,J,K)
  180         CONTINUE
  190       CONTINUE
  200     CONTINUE
        ENDIF
C
        CALL MCTIME('ITER  ',0,0,2)
        EMAXJ=0.0
        NEWMAT=EMAX.GT.ELIM1 .OR. SWITCH.NE.1.
        CALL MCTIME('ZERO W',0,0,1)
C
C  ZEROING OF MATRICES E AND W IN ROUTINE ZMAT
C
        IF(NEWMAT) THEN
          CALL ZMAT(W,NK*NK*NDEP)
        ENDIF
        CALL ZMAT(E,NK*NDEP)
        CALL MCTIME('ZERO W',0,1,1)
C
C  **** START OF DETAILED TRANSITIONS ****
C  CALCULATE THE ERROR TERM E BY SOLVING RADIATIVE TRANSFER IN DETAIL
C  CALCULATE MATRIX W
C
        CALL MCTIME('BMATTR',0,0,1)
        CALL REWIND(LOPC)
        CALL REWIND(LPHI)
        IREC=0
        DO 400 KR=1,NRAD
          LWARN=.FALSE.
          I=IRAD(KR)
          J=JRAD(KR)    
C
C  CALCULATE SOURCE FUNCTIONS, LINE OPACITIES AND CROSS-SECTIONS
C
          DO 270 K=1,NDEP
            DEIK(K)=0.
            EJ(KR,K)=0.0
  270     CONTINUE
          IF(KR.GT.NLINE) THEN
            CONT=.TRUE.
            DO 275 K=1,NDEP
              WPHI(K,KR)=1.0
  275       CONTINUE
          ELSE
            CONT=.FALSE.
            GIJK=G(I)/G(J)
            HN3C2=A(KR)/B(J,I)
            DO 280 K=1,NDEP
              Z(K)=N(I,K)-GIJK*N(J,K)
              IF(Z(K).LT.0.0) THEN
                IF(.NOT.LWARN) THEN
                  LWARN=.TRUE.
                  IF(IWARN.GE.2) WRITE(LJOBLO,278) KR,IT
  278             FORMAT(' NEGATIVE OPACITIES, KR=',I5,' IT=',I5)
                ENDIF
              ENDIF
              SL(K,KR)=HN3C2*N(J,K)*GIJK/Z(K)
              GIJ(K)=GIJK
  280       CONTINUE
          ENDIF
          IF(IWIDE(KR)) KT=KTRANS(KR)
          DO 380 NY=1,NQ(KR)
            IREC=IREC+1
            CALL READX
            CALL READJ(IREC)
            DO 285 K=1,NDEP
              SBCK(K)=SC(K)+SCAT(K)*JNY(K)
              JNYOLD(K)=JNY(K)
              JNY(K)=0.
  285       CONTINUE
            IF(CONT) THEN
              HN3C2=2.*HH*FRQ(NY,KT)/CC*FRQ(NY,KT)/CC*FRQ(NY,KT)
              DO 290 K=1,NDEP
                GIJ(K)=NSTAR(I,K)/NSTAR(J,K)*
     *                  EXP(-HH*FRQ(NY,KT)/BK/TEMP(K))
                ALFA(K)=ALFAC(NY,KR-NLINE)
                Z(K)=N(I,K)-GIJ(K)*N(J,K)
                SL(K,KR)=HN3C2*N(J,K)*GIJ(K)/Z(K)
  290         CONTINUE
            ENDIF
            DO 370 MU=1,NMU
              XMY=XMU(MU)
              WQMU=WQ(NY,KR)*WMU(MU)
              IF(.NOT.CONT .AND. (IND(KR).EQ.2 .OR. MU.EQ.1)) THEN
                CALL READP
                DO 300 K=1,NDEP
                  ALFA(K)=B(I,J)*PHI(K)*HNY4P
  300           CONTINUE
              ENDIF
              IF(IND(KR).EQ.2 .OR. MU.EQ.1) THEN
                DO 310 K=1,NDEP
                  X(K)=Z(K)*ALFA(K)/XNORM(K)+XCONT(K)
                  RNY(K)=XCONT(K)/X(K)
                  X(K)=SIGN(ONE,X(K))*MAX(ABS(X(K)),1.E-4*XCONT(K))
                  S(K)=(1.-RNY(K))*SL(K,KR)+RNY(K)*SBCK(K)
  310           CONTINUE
              ENDIF
C
C  CALCULATE INTENSITIES AND FILL NON-LOCAL ELEMENTS
C
              IF(INCRAD.NE.0) IMINUS(0)=XIMIN(NY,MU,KR)
              CALL TRANSP
              IF(NEWMAT) CALL BMAT(W,NK1,NDEP1,J,I)
C
C  RADIATION PART TO DIAGONAL IN W
C
              IF(CONT) THEN
                WQMUH=WQMU/HNY4P/FRQ(NY,KT)*FRQ(0,KT)
              ELSE
                WQMUH=WQMU/HNY4P
              ENDIF
              IF(NEWMAT) THEN
                DO 340 K=1,NDEP
                  W(I,J,K)=W(I,J,K)+WQMUH*ALFA(K)*GIJ(K)*
     *             RNY(K)*(P(K)+HN3C2)
                  W(J,I,K)=W(J,I,K)+WQMUH*ALFA(K)*
     *             RNY(K)*P(K)
  340           CONTINUE
              ENDIF
C
C  RADIATION PART TO ERROR E
C
              DO 360 K=1,NDEP
                DEIK(K)=DEIK(K)+WQMUH*ALFA(K)*(Z(K)*(PMS(K)+
     *           RNY(K)*SBCK(K))-N(J,K)*HN3C2*GIJ(K)*RNY(K))
                JNY(K)=JNY(K)+WMU(MU)*P(K)
  360         CONTINUE
  370       CONTINUE
            IF(IWTEST.LT.0) CALL WTEST(KR,NY)
            IF(IDLNY.LT.0) CALL WIDLNY(KR,NY)
C
C  CHECK CONVERGENCE OF JNY AND WRITE NEW JNY TO FILE
C
            DO 375 K=1,NDEP
CPGJ  IF-LOOP ADDED BY PGJ:
              IF(JNY(K) .GT. 0.)THEN
                DJ=(JNY(K)-JNYOLD(K))/JNY(K)
              ELSE
                DJ=0.
              ENDIF
              IF(ABS(DJ).GT.ABS(EJ(KR,K))) EJ(KR,K)=DJ
  375       CONTINUE
            CALL WRITEJ(IREC)
C
  380     CONTINUE
          DO 390 K=1,NDEP
            E(I,K)=E(I,K)+DEIK(K)*WPHI(K,KR)
            E(J,K)=E(J,K)-DEIK(K)*WPHI(K,KR)
  390     CONTINUE
          IF(IWTEST.LT.0 .AND. IWEVEC.GT.0) CALL WEVEC(E,NK1,NDEP1)
  400   CONTINUE
        CALL MCTIME('BMATTR',0,1,1)
C*
C* 89-03-18 START MODIFICATION
C  STORE KR INDEX FOR MAX CHANGE IN KREJ, K INDEX IN KEJ
C
        KREJ=0
        KEJ=0
        DO 404 KR=1,NRAD
          DO 402 K=1,NDEP
            IF(ABS(EJ(KR,K)).GT.EMAXJ) THEN
              KREJ=KR
              KEJ=K
              EMAXJ=ABS(EJ(KR,K))
            ENDIF
  402     CONTINUE
  404   CONTINUE
        IF(KREJ.GT.0 .AND. KEJ.GT.0) EMAXJ=EJ(KREJ,KEJ)
        WRITE(ETEXT,405) KREJ,KEJ
  405   FORMAT('EMAXJ(',I4,',',I4,')')
C
        IF(IT.GT.1) THEN
          CALL WCHANG(E,NK1,NDEP1,2)
          CALL WEMAX(ETEXT,EMAXJ)
        ENDIF
        EMAXJ=ABS(EMAXJ)
C* 89-03-18 END MODIFICATION
C*
C* CALCULATE ISUM
C*
        LISUM=ISUM.EQ.0
        IF(LISUM) THEN
          DO 408 K=1,NDEP
            POPMAX=N(1,K)
            ISUMA(K)=1
            DO 406 I=2,NK
              IF(N(I,K).GT.POPMAX) THEN
                POPMAX=N(I,K)
                ISUMA(K)=I
              ENDIF
  406       CONTINUE
  408     CONTINUE
        ENDIF
C
C  NORMALISE AND ADD INTO SMALL DIAGONALS
C
        IF(NEWMAT) THEN
          DO 440 I=1,NK
            DO 430 J=1,NK
              IF(KRAD(I,J).EQ.0) GOTO 430
              KR=KRAD(I,J)
              DO 410 K=1,NDEP
                W(I,J,K)=W(I,J,K)*WPHI(K,KR)*N(J,K)
                W(J,J,K)=W(J,J,K)-W(I,J,K)
  410         CONTINUE
  430       CONTINUE
  440     CONTINUE
          CALL MCTIME('NORM  ',0,1,1)
C
C  **** END OF DETAILED TRANSITIONS ****
C
C  FIXED TRANSITIONS
C  W MATRIX
C  C(I,I,K) HAS BEEN INITIALIZED TO ZERO IN COLRAT
C  AN IF-STATEMENT TO SKIP THE CASES WHEN I=J IS THUS
C  UNNECESSARY
C
          DO 575 K=1,NDEP
            DO 570 J=1,NK
              DO 560 I=1,NK
                W(I,I,K)=W(I,I,K)-C(I,J,K)*N(I,K)
                W(I,J,K)=W(I,J,K)+C(J,I,K)*N(J,K)
  560         CONTINUE
  570       CONTINUE
  575     CONTINUE
          CALL MCTIME('WCOLL ',0,1,1)
C
C  PARTICLE CONSERVATION EQUATION
C
          DO 600 K=1,NDEP
            IF(LISUM) ISUM=ISUMA(K)
            DO 590 J=1,NK
              W(ISUM,J,K)=N(J,K)
  590       CONTINUE
  600     CONTINUE
        ENDIF
C
C  ERROR VECTOR. FOR I=J SEE COMMENT ABOVE
C
        DO 700 K=1,NDEP
          DO 670 J=1,NK
            DO 660 I=1,NK
              E(I,K)=E(I,K)+N(I,K)*C(I,J,K)-N(J,K)*C(J,I,K)
  660       CONTINUE
  670     CONTINUE
          IF(LISUM) ISUM=ISUMA(K)
          E(ISUM,K)=TOTN(K)
          DO 690 J=1,NK
            E(ISUM,K)=E(ISUM,K)-N(J,K)
  690     CONTINUE
  700   CONTINUE
        IF(LISUM) ISUM=0
        CALL MCTIME('ERCOLL',0,1,1)
        IF(IWEVEC.GT.0) CALL WEVEC(E,NK1,NDEP1)
C
C  OUTPUT TO FILE WMAT
C
        IF(NEWMAT) CALL WWMAT(W,E,NK1,NDEP1,EMAX)
C
C  CALCULATE THE NEW POPULATION NUMBERS
C  PLOT LOGARITHMIC RELATIVE CHANGE
C
        CALL MCTIME('ENEQ  ',0,0,1)
        CALL ENEQ  (NK1,NDEP1,W,E,NEWMAT)
        CALL MCTIME('ENEQ  ',0,1,1)
C
        CALL WCHANG(E,NK1,NDEP1,1)
C
C  STORE I INDEX FOR MAX CHANGE IN IE, K INDEX IN KE
C
        EMAX=0.0
        IE=1
        KE=1
        DO 800 K=1,NDEP
          DO 780 I=1,NK
            IF(ABS(E(I,K)).GT.EMAX) THEN
              IE=I
              KE=K
              EMAX=ABS(E(I,K))
            ENDIF
C  ABORT IF NAN
            IF(N(I,K).NE.N(I,K)) THEN
              ICONV=0
              WRITE(ETEXT,810) I,K
              CALL WEMAX('      ',EMAX)
              CALL WEMAX(ETEXT,N(I,K))
              GOTO 9100
            ENDIF
  780     CONTINUE
  800   CONTINUE
        EMAX=E(IE,KE)
        WRITE(ETEXT,810) IE,KE
  810   FORMAT('EMAX (',I4,',',I4,')')
        CALL WNIIT(1)
        CALL WEMAX('      ',EMAX)
        IF(ICRSW.NE.0) CALL WEMAX('SWITCH        ',SWITCH)
        CALL WEMAX(ETEXT,EMAX)
        EMAX=ABS(EMAX)
C
C  IF SWITCHING IS ENABLED IN AUTOMATIC MODE:
C  CHECK TO SEE IF REDO IS TO BE DONE
C
        IF(ICRSW.EQ.-2) THEN
          IF(LOG10(EMAX).GT.EREDO .AND. OLDSW.NE.1.0) THEN
            DIVFAC=MAX(ONE+0.1,(1.1+(LOG10(EMAX)-EREDO)*2.))
            FS0=FS0/DIVFAC
            IF(FS0.LT.1.E-4) THEN
              CALL STOP('ITER: SWITCHING NOT CONVERGING')
            ENDIF
            SWITCH=MAX(ONE,OLDSW/(1.+FS0))
            IF(IWEMAX.NE.0) THEN
              WRITE(*,815) FS0
              WRITE(LOUT,815) FS0
  815         FORMAT(' REDOING ITERATION, FS0= ',1P,E10.2,0P)
            ENDIF
            GOTO 9000
          ENDIF
        ELSE IF(ICRSW.GT.0 .AND. IT.EQ.1 .AND. LOG10(EMAX).GT.EREDO)THEN
          SWITCH=SWITCH*100.
          GOTO 170
        ENDIF
C
        IF(ICRSW.NE.0) OLDSW=SWITCH
        IF(ICRSW.EQ.-2) THEN
          FS0=MIN(DFSMAX,(1.+(-LOG10(EMAX)+EREDO)*0.2))*FS0
          SWITCH=MAX(ONE,SWITCH/(1.+FS0))
        ELSE IF(ICRSW.EQ.-1) THEN
          WRITE(*,*) ' GIVE NEW SWITCH VALUE (<0 TO REDO ITERATION)'
          READ(*,*) SWITCH
          IF(SWITCH.LT.0.) THEN
            SWITCH=-SWITCH
            GOTO 9000
          ELSE IF(SWITCH.EQ.0.) THEN
            WRITE(*,*)    ' SWITCHING TO AUTOMATIC MODE'
            WRITE(LOUT,*) ' SWITCHING TO AUTOMATIC MODE'
            SWITCH=OLDSW
            ICRSW=-2
          ENDIF
        ELSE IF(ICRSW.GT.0) THEN
          IF(ABS(EMAX).LT.10.) THEN
            SWITCH=MAX(SWITCH*SWFAC,ONE)
          ELSE IF(ABS(EMAX).LT.100.) THEN
            SWITCH=MAX(SWITCH/100.*ABS(EMAX),ONE)
          ENDIF
        ENDIF
C
C  UPDATE POPULATIONS
C
        DO 830 K=1,NDEP
          DO 820 I=1,NK
            N(I,K)=N(I,K)*(1.0+E(I,K))
  820     CONTINUE
  830   CONTINUE
C
        IF(HSE .AND. EMAX.LT.ELIM1 .AND. SWITCH.EQ.1.) THEN
          CALL HSEINT(HSE)
          IF(ITHSE.EQ.1 .AND. NK+1.LE.MK) THEN
C
C  NG ACCELERATION FOR HSE ITERATIONS.
C  SWITCHED ON IF HSE STEP EACH ITERATION (ITHSE.EQ.1)
C  ZERO COUNT IN NG IF NGHSE=.FALSE.
C  ACCELERATION ON BOTH N AND NE THROUGH LOCAL ARRAY AN
C  TO HAVE ROOM FOR NE, NK+1.LE.MK IS REQUIRED
C
            IF(.NOT.NGHSE) CALL NG(N,NK,NDEP,NGHSE)
            NGHSE=.TRUE.
            DO 910 I=1,NK
              DO 900 K=1,NDEP
                AN(I,K)=N(I,K)
 900          CONTINUE
 910        CONTINUE
            DO 920 K=1,NDEP
              AN(NK+1,K)=NE(K)
 920        CONTINUE
            CALL NG(AN,NK+1,NDEP,NGHSE)
            DO 940 I=1,NK
              DO 930 K=1,NDEP
                N(I,K)=AN(I,K)
 930          CONTINUE
 940        CONTINUE
            DO 950 K=1,NDEP
              NE(K)=AN(NK+1,K)
 950        CONTINUE
          ELSE
            NGHSE=.FALSE.
          ENDIF
          EMAX=10.
          ITHSE=0
        ENDIF
        CALL MCTIME('ITER  ',0,2,2)
        IF(EMAX.LT.ELIM2 .AND. SWITCH.EQ.1.) THEN
          IF(EMAX.GT.0.0) ICONV=1
          GOTO 9100
        ENDIF
C
C
C  GIVE UP IF CORRECTION IS ABOVE 1000% AND SWITCHING
C
        IF(ICRSW.GT.0 .AND. IT.GE.1 .AND. EMAX.GT.10.) THEN
          ICONV=09
          GOTO 9100
        ENDIF
C
C  NG ACCELERATION STEP
C  START NG SCHEME AS SOON AS CORRECTIONS FALL BELOW 100 PERCENT
C  AND SWITCH=1.0. NOTE THAT NG COUNT IS RESET TO 1 AT EACH HSE
C  ITERATION
C
        IF(NGACC .AND. .NOT.NGHSE) THEN
          NGON=EMAX.LT.1.0 .AND. SWITCH.EQ.ONE
          CALL NG(N,NK,NDEP,NGON)
        ENDIF
C
        IF(SWITCH.EQ.1.) IT2=IT2+1
        IF(IT2.GT.ITMAX) GOTO 9100
 9000 CONTINUE
C
 9100 CONTINUE
      CALL MCTIME('SUM IT',0,2,3)
      RETURN
      END
C
C**********************************************************************
C
      SUBROUTINE ZMAT(W,N)
C
C  SETS ARRAY W(N) TO ZERO
C
C: ZMAT   90-07-31 NEW ROUTINE: (MATS CARLSSON)
C:
      INCLUDE 'PREC'
      DIMENSION W(N)
C
      DO 100 I=1,N
        W(I)=0.0
  100 CONTINUE
C
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE NG(AN,NK,NDEP,NGON)
C
C  PERFORMS NG ACCELERATION
C  SEE L.H. AUER , P. 101
C  IN KALKOFEN, ED., "NUMERICAL RADIATIVE TRANSFER", 
C  CAMBRIDGE UNIVERSITY PRESS 1987
C
C: NG     90-09-11  NEW ROUTINE: (MARTIN J. STIFT, MATS CARLSSON)
C:
C
      INCLUDE 'PREC'
      INCLUDE 'PARAM'
      INCLUDE 'CINPUT'
      INCLUDE 'CLU'
C
      DIMENSION AN(MK,MDEP)
      LOGICAL NGON
C
      DIMENSION YS(5,MK,MDEP),AA(3,3),BB(3)
      DATA ICALL/0/
      SAVE YS,ICALL
C
C  IF NGON=.FALSE., RESET ITERATION COUNT AND RETURN
C
      IF(.NOT.NGON) THEN
        ICALL=0
        RETURN
      ENDIF
      ICALL=ICALL+1
C
C STORE ITERATED LEVEL POPULATIONS FOR NG ACCELERATION
C
      IS0 = MOD(ICALL-1,5) + 1
      DO 910 I = 1,NK
        DO 900 K = 1,NDEP
          YS(IS0,I,K) = AN(I,K)
  900   CONTINUE
  910 CONTINUE
C
      IF (IS0.EQ.5) THEN
C
        IF(IWEMAX.NE.0) THEN
          WRITE(*,*) ' NG ACCELERATION'
          WRITE(LOUT,*) ' NG ACCELERATION'
        ENDIF
        DO 960 I = 1,NK
C
          DO 930 K1 = 1,3
            DO 920 K2 = 1,3
              AA(K1,K2) = 0.
  920       CONTINUE
            BB(K1) = 0.
  930     CONTINUE
C
          DO 940 K = 1,NDEP
            WT = 1. / YS(5,I,K)**2
            D0 = YS(5,I,K) - YS(4,I,K)
            D1 = D0 - YS(4,I,K) + YS(3,I,K)
            D2 = D0 - YS(3,I,K) + YS(2,I,K)
            D3 = D0 - YS(2,I,K) + YS(1,I,K)
            AA(1,1) = AA(1,1) + WT * D1 * D1
            AA(1,2) = AA(1,2) + WT * D1 * D2
            AA(1,3) = AA(1,3) + WT * D1 * D3
            AA(2,2) = AA(2,2) + WT * D2 * D2
            AA(2,3) = AA(2,3) + WT * D2 * D3
            AA(3,3) = AA(3,3) + WT * D3 * D3
            AA(2,1) = AA(1,2)
            AA(3,1) = AA(1,3)
            AA(3,2) = AA(2,3)
            BB(1)   = BB(1)   + WT * D0 * D1
            BB(2)   = BB(2)   + WT * D0 * D2
            BB(3)   = BB(3)   + WT * D0 * D3
  940     CONTINUE
C
          CALL LINEQ (AA,BB,3,3)

          DO 950 K = 1,NDEP
            AN(I,K) = (1. - BB(1) - BB(2) - BB(3)) * YS(5,I,K) +
     *                BB(1) * YS(4,I,K) + 
     *                BB(2) * YS(3,I,K) + 
     *                BB(3) * YS(2,I,K)
  950     CONTINUE
C
  960   CONTINUE
      ENDIF
C
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE LINEQ(A,B,N,M)
C
C  FINDS SOLUTION OF SYSTEM OF LINEAR EQUATIONS
C  WITH GAUSSIAN ELIMINATION WITH PIVOTING
C
C: LINEQ  90-06-05  NEW ROUTINE: (MARTIN J STIFT)
C:
C:        09-05-01  MODIFICATIONS: (MATS CARLSSON)
C:        MADE ROUTINE GENERAL AND NOT SPECIALIZED TO 3X3 MATRICES
C:
      INCLUDE 'PREC'
C
      PARAMETER (MDIM=10000)
      DIMENSION A(M,M),B(M),C(MDIM),ICOL(MDIM)
C
      IF(N.GT.MDIM) CALL STOP('LINEQ: N.GT.MDIM')
C
C INITIALIZE COLUMN COUNT AND STARTING POINTS
C
      DO 10 I = 1,N
        ICOL(I) = I
   10 CONTINUE
C
      IBEG = 1
      JBEG = 1
C
C DETERMINE PIVOT
C
   20 AMA = ABS(A(IBEG,JBEG))
      IMA = IBEG
      JMA = JBEG
      DO 30 I = IBEG,N
        DO 30 J = JBEG,N
          IF(ABS(A(I,J)).LE.AMA) GOTO 30
          AMA = ABS(A(I,J))
          IMA = I
          JMA = J
   30 CONTINUE
C
C ORDER MATRIX DEPENDING ON PIVOT
C
      DO 40 I = 1,N
        TEMP = A(I,JMA)
        A(I,JMA) = A(I,JBEG)
        A(I,JBEG) = TEMP
   40 CONTINUE
C
      DO 50 J = JBEG,N
        TEMP = A(IMA,J)
        A(IMA,J) = A(IBEG,J)
        A(IBEG,J) = TEMP
   50 CONTINUE
C
      TEMP = B(IMA)
      B(IMA) = B(IBEG)
      B(IBEG) = TEMP
      IT = ICOL(JBEG)
      ICOL(JBEG) = ICOL(JMA)
      ICOL(JMA)= IT
      IBEG = IBEG + 1
      JBEG = JBEG + 1
C
C ELIMINATE
C
      IMIN = IBEG - 1
      JMIN = JBEG - 1
      DO 70 I = IBEG,N
        QUOT = A(I,JMIN) / A(IMIN,JMIN)
        DO 60 J = JBEG,N
          A(I,J) = A(I,J) - QUOT * A(IMIN,J)
   60   CONTINUE
        B(I) = B(I) - QUOT * B(IMIN)
   70 CONTINUE
      IF (IBEG.LT.N) GOTO 20
C
C DETERMINE COEFFICIENTS
C
      B(N) = B(N) / A(N,N)
      N1 = N - 1
      DO 90 I = N1,1,-1
        I1 = I + 1
        DO 80 J = N,I1,-1
          B(I) = B(I) - A(I,J) * B(J)
   80   CONTINUE
        B(I) = B(I) / A(I,I)
   90 CONTINUE
C
C REORDER COEFFICIENTS
C
      DO 100 I = 1,N
        C(ICOL(I))=B(I)
  100 CONTINUE
C
      DO 110 I = 1,N
        B(I) = C(I)
  110 CONTINUE
C
      RETURN
      END
C
C**********************************************************************
C
      SUBROUTINE ENEQ(NK,NDEP,A,B,NEWMAT)
C
C  SOLVES THE EQUATION  SYSTEM  A*X=B.
C  WHEN NEWMAT=TRUE, THE SYSTEM IS REARRANGED INTO U*X=L*B, WHERE U
C  IS UPPER AND L IS LOWER TRIANGULAR. THESE ARE THEN REUSED IN LATER
C  CALLS WITH NEWMAT=FALSE AND NEW RIGHT HAND SIDES B. THE SOLUTION
C  VECTOR IS RETURNED IN B. NO PIVOTING, I.E. THE MATRIX A IS ASSUMED
C  TO HAVE NONZERO DIAGONAL ELEMENTS.
C
C  CODED BY: A. NORDLUND (FEB-1979)
C
C  THIS IS A MODIFIED VERSION OF EQSYST WHICH TESTS FOR ZERO ELEMENTS
C  BELOW THE DIAGONAL AND ALSO STOPS AT THE LAST NON-ZERO ELEMENT ABOVE
C  THE DIAGONAL. CONSIDERABLE SAVINGS ARE OBTAINED FOR LOOSE MATRICES.
C
C  THIS IS A COLUMN ORIENTED VERSION (M. CARLSSON JAN-1986)
C  TEMPORARY SCALARS BL, ALL, ALM AND BK ARE USED TO SHOW THE COMPILER
C  THAT THERE IS NO VECTOR DEPENDENCY IN THE INNERMOST DO-LOOP
C
C: ENEQ   90-06-08  NEW ROUTINE: (MARTIN J STIFT)
C:        LIKE EQSYST BUT FOR A BLOCK DIAGONAL GRAND MATRIX
C:
C:        03-08-12  MODIFICATIONS: (MATS CARLSSON)
C:        ADDED TEST OF DIMENSION MDIM
C:
      INCLUDE 'PREC'
      INCLUDE 'PARAM'
      PARAMETER (MDIM=1000)
C
      DIMENSION A(NK,NK,NDEP),B(NK,NDEP),LASTN(MDIM,MDIM)
      LOGICAL   NEWMAT
      SAVE      LASTN
C
      IF(NK.GT.MDIM) CALL STOP('ENEQ: NK.GT.MDIM')
      IF(NDEP.GT.MDIM) CALL STOP('ENEQ: NDEP.GT.MDIM')
      N = NK
C
      DO 120 MD = 1,NDEP
C
        IF (NEWMAT) THEN
C
C     FIND THE LAST NON-ZERO ELEMENT IN EACH COLUMN
C
          DO 30 L = 1,N
            DO 10 K = N,L+1,-1
              IF (A(K,L,MD).NE.0.0) GOTO 20
   10       CONTINUE
            K = L
   20       LASTN(L,MD) = K
   30     CONTINUE
C
C  COLUMN LOOP: ELIMINATE ELEMENTS BELOW THE DIAGONAL IN COLUMN L.
C
          DO 70 L = 1,N-1
C
C  STORE -A(K,L)/A(L,L) IN ELEMENT A(K,L)
C  MULTIPLY RIGHT HAND SIDE WITH -A(K,L)/A(L,L)
C
            ALL = A(L,L,MD)
            BL = B(L,MD)
            DO 40 K = L+1,LASTN(L,MD)
              A(K,L,MD) = -A(K,L,MD) / ALL
              B(K,MD) = B(K,MD) + A(K,L,MD) * BL
   40       CONTINUE
C
C  ADD FRACTION -A(K,L)/A(L,L) OF ROW L TO ROW K.
C  IN EACH COLUMN GO THROUGH ALL ROWS
C
            DO 60 M = L+1,N
              IF (A(L,M,MD).NE.0.0) THEN
                ALM = A(L,M,MD)
                LASTN(M,MD) = MAX(LASTN(L,MD),LASTN(M,MD))
                DO 50 K = L+1,LASTN(L,MD)
                  A(K,M,MD) = A(K,M,MD) + A(K,L,MD) * ALM
   50           CONTINUE
              END IF
   60       CONTINUE
C
   70     CONTINUE
C
        ELSE
C
          DO 90 L = 1,N-1
            BL = B(L,MD)
            DO 80 K = L+1,LASTN(L,MD)
              B(K,MD) = B(K,MD) + A(K,L,MD) * BL
   80       CONTINUE
   90     CONTINUE
C
        END IF
C
C BACKSUBSTITUTE
C
        DO 110 K = N,1,-1
          BK = B(K,MD)
          DO 100 L = K+1,N
            BK = BK - A(K,L,MD) * B(L,MD)
  100     CONTINUE
          B(K,MD) = BK / A(K,K,MD)
  110   CONTINUE
C
  120 CONTINUE
C
      RETURN
      END
C
C**********************************************************************
C
      SUBROUTINE WWMAT(W,E,NK1,NDEP1,EMAX)
C
C  WRITES MATRIX W TO UNFORMATTED FILE WMAT
C  THE MATRIX IS OUTPUT COLUMN BY COLUMN WITH ONE RECORD
C  PER COLUMN.
C  THE FIRST RECORD IS AN ID RECORD CONTAINING:
C  NDEP,NK,NLINE,NRAD,(IRAD(KR),KR=1,NRAD),(JRAD(KR),KR=1,NRAD),
C  ISUM,ATOMID,ATMOID,DPID
C
C  IF IWWMAT.GT.0 THE FILE IS WRITTEN ONLY WHEN EMAX.LT.ELIM1*100
C  TO AVOID TOO MANY WRITE OPERATIONS. IF EMAX GOES FROM A VALUE
C  LARGER THAN ELIM1*100 TO A VALUE SMALLER THAN ELIM1 IN ONE
C  ITERATION THIS MEANS THAT THERE WILL BE NO WMAT OUTPUT.
C:
C: WWMAT  88-02-26  MODIFICATIONS: (MATS CARLSSON)
C:        WRITES RIGHT HAND SIDE TO FILE IN LAST RECORD
C:
C:        90-06-05  MODIFICATIONS: (MARTIN J STIFT)
C:        LOCAL OPERATOR VERSION
C:
C:        90-08-14  MODIFICATIONS: (MATS CARLSSON)
C:        CHANGED TO CONTAIN DECLARATION OF W AND E
C:        THESE PASSED IN ARGUMENT LIST INSTEAD OF IN COMMON TO
C:        AVOID NK=MK1 AND NDEP=MDEP1 REQUIREMENT
C:
C:        90-10-19  MODIFICATIONS: (MATS CARLSSON)
C:        DOUBLE PRECISION VERSION. MATRIX OUTPUT IN SINGLE
C:        PRECISION
C:
      INCLUDE 'PREC'
      INCLUDE 'PARAM'
      INCLUDE 'CATOM'
      INCLUDE 'CATMOS'
      INCLUDE 'CINPUT'
      INCLUDE 'CLU'
C
      DIMENSION W(NK1,NK1,NDEP1),E(NK1,NDEP1)
C
      IF(IWWMAT.EQ.0) RETURN
C
      CALL MCTIME('WWMAT ',0,0,1)
      IF(IWWMAT.GT.0) THEN
        IF(EMAX.LT.ELIM1*100.) THEN
          CALL REWIND(LWMAT)
          WRITE(LWMAT) NDEP,NK,NLINE,NRAD,ISUM,(IRAD(KR),KR=1,NRAD),
     *    (JRAD(KR),KR=1,NRAD),ATOMID,ATMOID,DPID
          WRITE(LWMAT) ((( REAL(W(IK,JK,JL)),IK=1,NK),
     *     JK=1,NK),JL=1,NDEP)
          WRITE(LWMAT) (( REAL(E(IK,JL)),IK=1,NK),JL=1,NDEP)
        ENDIF
      ELSE
        IF(IT.EQ.1) THEN
          WRITE(LWMAT) NDEP,NK,NLINE,NRAD,ISUM,(IRAD(KR),KR=1,NRAD),
     *    (JRAD(KR),KR=1,NRAD),ATOMID,ATMOID,DPID
        ENDIF
        WRITE(LWMAT) ((( REAL(W(IK,JK,JL)),IK=1,NK),
     *   JK=1,NK),JL=1,NDEP)
        WRITE(LWMAT) (( REAL(E(IK,JL)),IK=1,NK),JL=1,NDEP)
      ENDIF
      CALL MCTIME('WWMAT ',0,1,1)
C
      RETURN
      END
C
C***********************************************************************
C
