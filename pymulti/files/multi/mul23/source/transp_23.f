      SUBROUTINE TRANSP
C
C  SOLVES THE RADIATIVE TRANSFER EQUATION WITH GIVEN SOURCE FUNCTION.
C  SOLVES FOR P-S WHERE P IS FEAUTRIERS P; I.E.
C
C    P = 0.5*(I(XMU)+I(-XMU))
C
C  IN THE PRESENCE OF VELOCITY FIELDS:
C
C    P = 0.5*(I(NY,XMU) + I(-NY,-XMU))
C
C  ITRAN DETERMINES THE MODE OF FORMAL SOLUTION
C  ITRAN  = 0  FEAUTRIER SOLUTION
C           1  FEAUTRIER SOLUTION TO CUBIC SPLINE ACCURACY, REF:
C              KUNASZ, HUMMER, 1974, MNRAS 166,19
C              MIHALAS, 1974, APJ SUPPL 28,343
C           2  FEAUTRIER SOLUTION HERMITE, REF:
C              AUER, 1976, JQSRT 16,931
C
C           3  INTEGRAL CUBIC SPLINE METHOD ACCORDING TO SCHARMER
C           4  INTEGRAL CUBIC SPLINE METHOD WITH LOCALLY DETERMINED
C              FIRST DERIVATIVE
C
C  VARIABLES:
C
C  XMU   : ANGULAR COSINE                                       (IN)
C  NDEP  : NUMBER OF DEPTH POINTS                               (IN)
C  TAU   : STANDARD OPTICAL DEPTH SCALE                         (IN)
C  X     : RATIO OF MONOCHROMATIC TO STANDARD OPACITY           (IN)
C  S     : MONOCHROMATIC SOURCE FUNCTION                        (IN)
C  ITRAN : DETERMINES MODE OF FORMAL SOLUTION                   (IN)
C
C  P     : MEAN BIDIRECTIONAL INTENSITY (CF. ABOVE)            (OUT)
C  PMS   : P-S                                                 (OUT)
C  IPLUS : RADIATION INTENSITY, OUTGOING RAYS                  (OUT)
C  IMINUS: RADIATION INTENSITY, INGOING RAYS                   (OUT)
C  TAUQ  : MONOCHROMATIC OPTICAL DEPTH                         (OUT)
C  DTAUQ : DTAUQ(K)=TAUQ(K)-TAUQ(K-1)                          (OUT)
C  A1    : 1/(DTAUQ(K+0.5)*DTAUQ(K))                           (OUT)
C  C1    : 1/(DTAUQ(K+0.5)*DTAUQ(K+1))                         (OUT)
C
C:
C: TRANSP 86-09-04  MODIFICATIONS: (GORAN SCHARMER, MATS CARLSSON)
C:        INCLUDES CALL TO CUBIC SPLINE INTEGRAL FORMAL SOLVER
C:
C:        92-10-08  MODIFICATIONS: (MATS CARLSSON)
C:        ITRAN=10-14 ALLOWED (ITRAN.GT.10 INCLUDES INCOMING RADIATION)
C:
C:        95-08-17  MODIFICATIONS: (MATS CARLSSON)
C:        INCOMING RADIATION SWITCHED ON BY INCRAD.NE.0 INSTEAD OF
C:        ITRAN=10-14
C:
      INCLUDE 'PREC'
      INCLUDE 'PARAM'
      INCLUDE 'CATMOS'
      INCLUDE 'CTRAN'
      INCLUDE 'CSLINE'
      INCLUDE 'CINPUT'
C
C K=1: UPPER BOUNDARY
C
      CMU=0.5/XMY
      DTAUQ(2)=(X(1)+X(2))*(TAU(2)-TAU(1))*CMU
      A1(1)=1./DTAUQ(2)
      T=TAU(1)*X(1)*2.0*CMU
      TAUQ(1)=T
      DTAUQ(1)=T
C
C  CALCULATE DTAUQ 
C
      IF(DPTYPE.EQ.'H') THEN
        DO 100 K=2,NDEP
          DTAUQ(K)=(X(K)*XNORM(K)+X(K-1)*XNORM(K-1))*
     *     (HEIGHT(K-1)-HEIGHT(K))*1.E5*CMU
  100   CONTINUE
C        DO 102 K=2,3
C          TAUQ(K)=TAUQ(K-1)+DTAUQ(K)
C  102   CONTINUE
C        TAUQ(1)=EXP(2*LOG(TAUQ(2))-LOG(TAUQ(3)))
      ELSE
        DO 105 K=2,NDEP
          DTAUQ(K)=(X(K)+X(K-1))*(TAU(K)-TAU(K-1))*CMU
  105   CONTINUE
      ENDIF
C
C  CALCULATE TAUQ
C
      DO 110 K=2,NDEP
        TAUQ(K)=TAUQ(K-1)+DTAUQ(K)
  110 CONTINUE
C
C  CALCULATE A1 AND C1
C
      DO 120 K=2,NDEP-1
        A1(K)=2./(DTAUQ(K)+DTAUQ(K+1))/DTAUQ(K)
        C1(K)=2./(DTAUQ(K)+DTAUQ(K+1))/DTAUQ(K+1)
  120 CONTINUE
C
C  CHOOSE THE METHOD OF FORMAL SOLUTION
C
      IF(ISCAT.EQ.1) THEN
        CALL TRANSC
      ELSE IF(ITRAN.LE.2) THEN
        CALL TRANF
      ELSE
        CALL TRANI
      ENDIF
C
C CALCULATE P(K)-S(K)
C
      PMS(1)=P(1)-S(1)
      DO 130 K=2,NDEP-1
        IF(ABS(A1(K)).GT.1.0) THEN
          PMS(K)=P(K)-S(K)
        ELSE
          PMS(K)=C1(K)*(P(K+1)-P(K))-A1(K)*(P(K)-P(K-1))
        END IF
  130 CONTINUE
      PMS(NDEP)=(P(NDEP-1)-P(NDEP)+S(NDEP)-S(NDEP-1))
     * /(DTAUQ(NDEP)+0.5*DTAUQ(NDEP)**2)
C
      RETURN
      END
C
C**********************************************************************
C
      SUBROUTINE TRANF
C
C  FORMAL SOLVER FOR THE TRANSFER EQUATION USING FEAUTRIER TECHNIQUE
C
C  ITRAN DETERMINES THE MODE OF FORMAL SOLUTION
C        =  0  FEAUTRIER SOLUTION
C           1  FEAUTRIER SOLUTION TO CUBIC SPLINE ACCURACY, REF:
C              AUER, 1976, JQSRT 16,931
C           2  FEAUTRIER SOLUTION HERMITE
C
C  THE  UPPER AND LOWER  BOUNDARY CONDITIONS ARE DERIVED FROM SECOND
C  ORDER TAYLOR EXPANSIONS, WITH ESTIMATES OF THE INCIDENT RADIATION
C  FIELDS:
C
C    P(2)    =  P(1) + DTAU*P'(1) + 0.5*DTAU**2*P''(1)
C    P'(1)   =  P(1) - I(-XMU)
C    I(-XMU) =  S(1)*(1.-EXP(-TAUQ(1)))
C
C    P(N-1)  =  P(N) - DTAU*P'(N) + 0.5*DTAU**2*P''(N)
C    P'(N)   =  I(XMU) - P(N)
C    I(XMU)  =  S(N) + (S(N)-S(N-1))/DTAU
C
C  STEINS  TRICK
C  (STORING THE SUM OF ELEMENTS INSTEAD OF THE DIAGONAL) IS USED  TO
C  AVOID  NUMERICAL  ROUND-OFF  PROBLEMS  AT  SMALL  OPTICAL DEPTHS.
C
C  CODED BY: A.NORDLUND (OCT-1981). REVISED (MAR-1982).
C
C: TRANF  92-10-08  MODIFICATIONS: (MATS CARLSSON)
C:        INCOMING RADIATION FIELD IS SPECIFIED IN IMINUS(0)
C:        THIS IS INDICATED BY ITRAN.GT.10
C:
C:        95-08-17  MODIFICATIONS: (MATS CARLSSON)
C:        INCOMING RADIATION FIELD INDICATED BY INCRAD.NE.0 
C:
C:        97-06-27  MODIFICATIONS: (MATS CARLSSON)
C:        SAVE COEFFICIENTS FOR LATER USE IN LOCAL OPERATOR
C:        FOLLOWS CODING BY MARTIN STIFT
C:        
      INCLUDE 'PREC'
      INCLUDE 'PARAM'
      INCLUDE 'CATMOS'
      INCLUDE 'CTRAN'
      INCLUDE 'CSLINE'
      INCLUDE 'CINPUT'
C
      COMMON /CRYBH/ SP21,SP31,SP1N,SP2N
C
      DIMENSION SP1(MDEP),SP2(MDEP),SP3(MDEP)
C
      T=TAUQ(1)
      SP1(1)=0.
      SP2(1)=1.+2.*A1(1)
      SP3(1)=-2.*A1(1)*A1(1)
      IF (T.LT.0.01) THEN
        EX1=T*(1.-T*(0.5-T*(0.1666667-T*0.041666667)))
      ELSE
        IF(T.LT.20.) THEN
          EX1=1.-EXP(-T)
        ELSE
          EX1=1.
        END IF
      END IF
      EX=1.-EX1
      FACT=1.+2.*A1(1)*EX1
      SP2(1)=SP2(1)/FACT
      SP3(1)=SP3(1)/FACT
      IF(INCRAD.EQ.0) THEN
        IMINUS(0)=0.0
        P(1)=S(1)
      ELSE
        P(1)=S(1)+2.*A1(1)/FACT*EXP(-T)*IMINUS(0)
      ENDIF
      SP31   = SP3(1)
      SP21   = SP2(1)
C
C  CALCULATE TRIDIAGONAL COEFFICIENTS
C
      IF(ITRAN.EQ.0) THEN
        DO 200 K=2,NDEP-1
          SP1(K)=-A1(K)
          SP2(K)=1.
          SP3(K)=-C1(K)
          P(K)=S(K)
  200   CONTINUE
      ELSE IF(ITRAN.EQ.1) THEN
        DO 201 K=2,NDEP-1
          AD=.166666666*DTAUQ(K)*2./(DTAUQ(K)+DTAUQ(K+1))
          CD=.166666666*DTAUQ(K+1)*2./(DTAUQ(K)+DTAUQ(K+1))
          SP1(K)=-A1(K)+AD
          SP2(K)=1.
          SP3(K)=-C1(K)+CD
          P(K)=S(K)+AD*(S(K-1)-S(K))+CD*(S(K+1)-S(K))
  201   CONTINUE
      ELSE IF(ITRAN.EQ.2) THEN
        DO 202 K=2,NDEP-1
          AD=.166666666-0.083333333*DTAUQ(K+1)**2/DTAUQ(K)*
     *     2./(DTAUQ(K)+DTAUQ(K+1))
          CD=.166666666-0.083333333*DTAUQ(K)**2/DTAUQ(K+1)*
     *     2./(DTAUQ(K)+DTAUQ(K+1))
          SP1(K)=-A1(K)+AD
          SP2(K)=1.
          SP3(K)=-C1(K)+CD
          P(K)=S(K)+AD*(S(K-1)-S(K))+CD*(S(K+1)-S(K))
  202   CONTINUE
      ELSE
        CALL STOP(' TRANF: ITRAN OUTSIDE RANGE')
      ENDIF
C
C K=NDEP: LOWER BOUNDARY
C
      SP1(NDEP)=-1.
      SP2(NDEP)=DTAUQ(NDEP)+0.5*DTAUQ(NDEP)**2
      SP3(NDEP)=0.
      P(NDEP)=S(NDEP)*(DTAUQ(NDEP)+0.5*DTAUQ(NDEP)**2)+
     *       (S(NDEP)-S(NDEP-1))
      SP1N=SP1(NDEP)/SP2(NDEP)
      SP2N=SP2(NDEP)/SP2(NDEP)
C
C ELIMINATE SUBDIAGONAL
C
      DO 300 K=1,NDEP-1
        F=-SP1(K+1)/(SP2(K)-SP3(K))
        P(K+1)=P(K+1)+F*P(K)
        SP2(K+1)=SP2(K+1)+F*SP2(K)
        SP2(K)=SP2(K)-SP3(K)
  300 CONTINUE
      SP2(NDEP)=SP2(NDEP)-SP3(NDEP)
C
C BACKSUBSTITUTE
C
      P(NDEP)=P(NDEP)/SP2(NDEP)
      DO 320 K=NDEP-1,1,-1
        P(K)=(P(K)-SP3(K)*P(K+1))/SP2(K)
  320 CONTINUE
C
C SURFACE INTENSITY
C
      IPLUS(0)=2.*(EX*P(1)+S(1)*0.5*EX1**2) - IMINUS(0)*EX*EX
C
C  INTENSITIES, CALCULATED FROM A WEIGHTED AVERAGE OF
C  DP/DTAUNY(K+0.5) AND DP/DTAUNY(K-0.5)
C
      IMINUS(1)=EX1*S(1) + EXP(-T)*IMINUS(0)
      IPLUS(1)=2.*P(1)-IMINUS(1)
      DPDT2=(P(2)-P(1))/DTAUQ(2)
      DO 340 K=2,NDEP-1
        DPDT1=DPDT2
        DPDT2=(P(K+1)-P(K))/DTAUQ(K+1)
        PPRIMK=(DTAUQ(K)*DPDT2+DTAUQ(K+1)*DPDT1)/
     *         (DTAUQ(K)+DTAUQ(K+1))
        IPLUS(K)=P(K)+PPRIMK
        IMINUS(K)=P(K)-PPRIMK
  340 CONTINUE
      IPLUS(NDEP)=S(NDEP)+(S(NDEP)-S(NDEP-1))/DTAUQ(NDEP)
      IMINUS(NDEP)=2.0*P(NDEP)-IPLUS(NDEP)
C
      RETURN
      END
C
C**********************************************************************
C
      SUBROUTINE TRANSC
C
C  FORMAL SOLVER FOR THE TRANSFER EQUATION USING FEAUTRIER TECHNIQUE
C  VERSION SOLVING FOR ALL ANGLES AT ONCE TO BE USED IF
C  SCATTERING DOMINATES
C  FOR THIS ROUTINE TO WORK IT MUST BE POSSIBLE TO SEPARATE THE
C  ABSORPTION AND SCATTERING PARTS OF THE SOURCE FUNCTION
C  THIS IS DONE BY WRITING THE SOURCE FUNCTION AS
C  S = SABS + SSCAT*JNY
C  SABS AND SSCAT ARE EXTRACTED FROM THE VARIABLES SBCK, SC, SCAT
C  AND RNY FROM:  
C  SSCAT(K)=RNY(K)*SCAT(K)
C  SABS(K)=S(K)-RNY(K)*SBCK(K)+RNY(K)*SC(K)
C  THESE VARIABLES THUS HAVE TO BE LOADED WITH THE CORRECT VALUES
C
C  ROUTINE HAS NOT BEEN OPTIMIZED, MAJOR SAVINGS PROBABLY POSSIBLE
C  ITRAN=0 IS THE ONLY FEAUTRIER OPTION IMPLEMENTED
C
C: TRANSC 95-02-22  NEW ROUTINE: (MATS CARLSSON)
C:
      INCLUDE 'PREC'
      INCLUDE 'PARAM'
      INCLUDE 'CATMOS'
      INCLUDE 'CTRAN'
      INCLUDE 'CSLINE'
      INCLUDE 'CINPUT'
      INCLUDE 'CGAUSI'
      INCLUDE 'CLU'
C
      DIMENSION SP1(MMU,MMU,MDEP),SP2(MMU,MMU,MDEP),SP3(MMU,MMU,MDEP)
      DIMENSION PSC(MMU,MDEP),SCJ(MDEP),SABS(MDEP),SSCAT(MDEP)
      DIMENSION WORK(MMU,MMU),FWORK(MMU,MMU),PWORK(MMU),PWORK2(MMU)
      LOGICAL LWARN
      DATA LWARN/.FALSE./
      SAVE PSC,SCJ,SABS,SSCAT,LWARN
      COMMON /CRYBH/ SP21,SP31,SP1N,SP2N
C
C  CALCULATE SP21,SP31,SP1N AND SP2N FOR RYBICKI-HUMMER OPERATOR
C
      T=TAUQ(1)
      SP21=1.+2.*A1(1)
      SP31=-2.*A1(1)*A1(1)
      IF (T.LT.0.01) THEN
        EX1=T*(1.-T*(0.5-T*(0.1666667-T*0.041666667)))
      ELSE
        IF(T.LT.20.) THEN
          EX1=1.-EXP(-T)
        ELSE
          EX1=1.
        END IF
      END IF
      EX=1.-EX1
      FACT=1.+2.*A1(1)*EX1
      SP21=SP21/FACT
      SP31=SP31/FACT
      SP1N=-1.
      SP2N=DTAUQ(NDEP)+0.5*DTAUQ(NDEP)**2
      SP1N=SP1N/SP2N
      SP2N=SP2N/SP2N
C
C  SOLVE FOR ALL ANGLES WHEN MU=1 AND STORE RESULTS IN LOCAL VARIABLE
C  PSC
C  FOR ALL MU-VALUES, COPY FROM LOCAL VARIABLE INTO MU-DEPENDENT 
C  P
C  QUANTITIES THAT INCLUDE XMY HAS TO BE ADJUSTED
C
C  IF ITRAN HAS A NON-IMPLEMENTED VALUE, ISSUE A WARNING BUT ONLY ONCE
C
      IF(ITRAN.NE.0 .AND. .NOT.LWARN) THEN
        WRITE(LJOBLO,*) 'TRANSC: WARNING - ONLY ITRAN=0 IMPLEMENTED'
        LWARN=.TRUE.
      ENDIF
      DO 100 MU=1,NMU
        IF(XMY.EQ.XMU(MU)) GOTO 110
  100 CONTINUE
      CALL STOP('TRANSC: MU NOT FOUND')
  110 CONTINUE
      MU0=MU
      IF(MU0.EQ.1) THEN
        DO 120 K=1,NDEP
          SSCAT(K)=RNY(K)*SCAT(K)
          SABS(K)=S(K)-RNY(K)*SBCK(K)+RNY(K)*SC(K)
  120   CONTINUE
        DO 280 MU=1,NMU
          DO 150 MU2=1,NMU
            DO 140 K=1,NDEP
              SP1(MU,MU2,K)=0.0
              SP2(MU,MU2,K)=0.0
              SP3(MU,MU2,K)=0.0
  140       CONTINUE
  150     CONTINUE
          XMU1=XMU(MU)/XMU(1)
          XMU2=XMU1*XMU1
          T=TAUQ(1)/XMU1
          SP1(MU,MU,1)=0.
          SP2(MU,MU,1)=1.+2.*XMU1*A1(1)
          SP3(MU,MU,1)=-2.*XMU2*A1(1)*A1(1)
          IF (T.LT.0.01) THEN
            EX1=T*(1.-T*(0.5-T*(0.1666667-T*0.041666667)))
          ELSE
            IF(T.LT.20.) THEN
              EX1=1.-EXP(-T)
            ELSE
              EX1=1.
            END IF
          END IF
          EX=1.-EX1
          FACT=1.+2.*XMU1*A1(1)*EX1
          SP2(MU,MU,1)=SP2(MU,MU,1)/FACT
          SP3(MU,MU,1)=SP3(MU,MU,1)/FACT
          IF(INCRAD.EQ.0) THEN
            IMINUS(0)=0.0
            PSC(MU,1)=SABS(1)
          ELSE
            PSC(MU,1)=SABS(1)+2.*XMU1*A1(1)/FACT*EXP(-T)*IMINUS(0)
          ENDIF
C
C  CALCULATE TRIDIAGONAL COEFFICIENTS
C
          DO 200 K=2,NDEP-1
            SP1(MU,MU,K)=-XMU2*A1(K)
            SP2(MU,MU,K)=1.
            SP3(MU,MU,K)=-XMU2*C1(K)
            PSC(MU,K)=SABS(K)
  200     CONTINUE
C
C K=NDEP: LOWER BOUNDARY
C
          SP1(MU,MU,NDEP)=-1.
          SP2(MU,MU,NDEP)=DTAUQ(NDEP)/XMU1+0.5*DTAUQ(NDEP)**2/XMU2
          SP3(MU,MU,NDEP)=0.
          SK=(DTAUQ(NDEP)/XMU1+0.5*DTAUQ(NDEP)**2/XMU2)+1.0
          PSC(MU,NDEP)=SABS(NDEP)*SK - SABS(NDEP-1)
C
C NON-DIAGONAL SCATTERING ELEMENTS
C
          DO 250 MU2=1,NMU
            SP2(MU,MU2,1)=SP2(MU,MU2,1)-SSCAT(1)*WMU(MU2)
            DO 220 K=2,NDEP-1
              SP2(MU,MU2,K)=SP2(MU,MU2,K)-SSCAT(K)*WMU(MU2)
  220       CONTINUE
            SP2(MU,MU2,NDEP)=SP2(MU,MU2,NDEP)-SSCAT(NDEP)*WMU(MU2)*SK+
     *       SSCAT(NDEP-1)*WMU(MU2)
            SP1(MU,MU2,NDEP)=SP1(MU,MU2,NDEP)+SSCAT(NDEP-2)*WMU(MU2)
  250     CONTINUE
  280   CONTINUE
C
C ELIMINATE SUBDIAGONAL
C
        DO 390 K=1,NDEP-1
          CALL MATSUB(SP2(1,1,K),SP3(1,1,K),WORK,NMU,NMU,MMU,MMU)
          CALL MATINV(WORK,NMU,MMU)
          CALL MATMUL(SP1(1,1,K+1),WORK,FWORK,NMU,NMU,NMU,MMU,MMU,MMU)
          CALL MATMUL(FWORK,PSC(1,K),PWORK,NMU,NMU,1,MMU,MMU,1)
          CALL MATMUL(FWORK,SP2(1,1,K),WORK,NMU,NMU,NMU,MMU,MMU,MMU)
          DO 380 MU=1,NMU
            PSC(MU,K+1)=PSC(MU,K+1)-PWORK(MU)
            DO 370 MU2=1,NMU
              SP2(MU,MU2,K+1)=SP2(MU,MU2,K+1)-WORK(MU,MU2)
              SP2(MU,MU2,K)=SP2(MU,MU2,K)-SP3(MU,MU2,K)
  370       CONTINUE
  380     CONTINUE
  390   CONTINUE
        DO 420 MU=1,NMU
          DO 410 MU2=1,NMU
            SP2(MU,MU2,NDEP)=SP2(MU,MU2,NDEP)-SP3(MU,MU2,NDEP)
  410     CONTINUE
  420   CONTINUE
C
C BACKSUBSTITUTE
C
        CALL MATINV(SP2(1,1,NDEP),NMU,MMU)
        CALL MATMUL(SP2(1,1,NDEP),PSC(1,NDEP),PWORK,NMU,NMU,1,MMU,MMU,1)
        SCJ(NDEP)=0.0
        DO 500 MU=1,NMU
          PSC(MU,NDEP)=PWORK(MU)
          SCJ(NDEP)=SCJ(NDEP)+WMU(MU)*PSC(MU,NDEP)
  500   CONTINUE
        DO 690 K=NDEP-1,1,-1
          CALL MATMUL(SP3(1,1,K),PSC(1,K+1),PWORK,NMU,NMU,1,MMU,MMU,1)
          CALL MATSUB(PSC(1,K),PWORK,PWORK2,NMU,1,MMU,1)
          CALL MATINV(SP2(1,1,K),NMU,MMU)
          CALL MATMUL(SP2(1,1,K),PWORK2,PWORK,NMU,NMU,1,MMU,MMU,1)
          SCJ(K)=0.0
          DO 650 MU=1,NMU
            PSC(MU,K)=PWORK(MU)
            SCJ(K)=SCJ(K)+WMU(MU)*PSC(MU,K)
  650     CONTINUE
  690   CONTINUE
      ENDIF
C
C  MOVE LOCAL VARIABLES INTO MULTI VARIABLES
C
      DO 700 K=1,NDEP
        P(K)=PSC(MU0,K)
        S(K)=SABS(K)+SSCAT(K)*SCJ(K)
  700 CONTINUE
C
      T=TAUQ(1)
      IF (T.LT.0.01) THEN
        EX1=T*(1.-T*(0.5-T*(0.1666667-T*0.041666667)))
      ELSE
        IF(T.LT.20.) THEN
          EX1=1.-EXP(-T)
        ELSE
          EX1=1.
        END IF
      END IF
      EX=1.-EX1
C
C SURFACE INTENSITY
C
      IPLUS(0)=2.*(EX*P(1)+S(1)*0.5*EX1**2) - IMINUS(0)*EX*EX
C
C  INTENSITIES, CALCULATED FROM A WEIGHTED AVERAGE OF
C  DP/DTAUNY(K+0.5) AND DP/DTAUNY(K-0.5)
C
      IMINUS(1)=EX1*S(1) + EXP(-T)*IMINUS(0)
      IPLUS(1)=2.*P(1)-IMINUS(1)
      DPDT2=(P(2)-P(1))/DTAUQ(2)
      DO 800 K=2,NDEP-1
        DPDT1=DPDT2
        DPDT2=(P(K+1)-P(K))/DTAUQ(K+1)
        PPRIMK=(DTAUQ(K)*DPDT2+DTAUQ(K+1)*DPDT1)/
     *         (DTAUQ(K)+DTAUQ(K+1))
        IPLUS(K)=P(K)+PPRIMK
        IMINUS(K)=P(K)-PPRIMK
  800 CONTINUE
      IPLUS(NDEP)=S(NDEP)+(S(NDEP)-S(NDEP-1))/DTAUQ(NDEP)
      IMINUS(NDEP)=2.0*P(NDEP)-IPLUS(NDEP)
C
      RETURN
      END
C
C**********************************************************************
C
      SUBROUTINE TRANI
C
C  FORMAL SOLVER FOR THE TRANSFER EQUATION USING DIRECT INTEGRATION
C  FOR CUBIC SPLINE SOURCES
C
C  ITRAN DETERMINES THE MODE OF FORMAL SOLUTION
C         = 3  INTEGRAL CUBIC SPLINE METHOD ACCORDING TO SCHARMER
C           4  INTEGRAL CUBIC SPLINE METHOD WITH LOCALLY DETERMINED
C              FIRST DERIVATIVE
C
C  REF G.B.SCHARMER 
C
C:
C: TRANI  86-09-04  NEW ROUTINE: (GORAN SCHARMER)
C:        INTEGRAL CUBIC SPLINE FORMAL SOLVER
C:
C:        92-10-08  MODIFICATIONS: (MATS CARLSSON)
C:        INCOMING RADIATION FIELD IS SPECIFIED IN IMINUS(0)
C:        THIS IS INDICATED BY ITRAN.GT.10
C:
C:        95-08-17  MODIFICATIONS: (MATS CARLSSON)
C:        INCOMING RADIATION FIELD INDICATED BY INCRAD.NE.0 
C:
C:        09-05-05  MODIFICATIONS: (MATS CARLSSON)
C:        ADDED ITRAN=5 (I+-0 AT BOTTOM BOUNDARY, OTHERWISE AS ITRAN=4)
C:
      INCLUDE 'PREC'
      INCLUDE 'PARAM'
      INCLUDE 'CATMOS'
      INCLUDE 'CTRAN'
      INCLUDE 'CSLINE'
      INCLUDE 'CINPUT'
C
      DIMENSION SPRIM(MDEP),SBISS(MDEP),STRISS(MDEP),
     * EXPD(MDEP),WW(MDEP)
      INTEGER KP(MDEP)
      DIMENSION AIP(0:MDEP),AIM(0:MDEP)
      DATA DT1/1.E-2/,DT2/1.E2/
      COMMON /CRYBH/ SP21,SP31,SP1N,SP2N
C
C  CALCULATE SP21,SP31,SP1N AND SP2N FOR RYBICKI-HUMMER OPERATOR
C
      T=TAUQ(1)
      SP21=1.+2.*A1(1)
      SP31=-2.*A1(1)*A1(1)
      IF (T.LT.0.01) THEN
        EX1=T*(1.-T*(0.5-T*(0.1666667-T*0.041666667)))
      ELSE
        IF(T.LT.20.) THEN
          EX1=1.-EXP(-T)
        ELSE
          EX1=1.
        END IF
      END IF
      EX=1.-EX1
      FACT=1.+2.*A1(1)*EX1
      SP21=SP21/FACT
      SP31=SP31/FACT
      SP1N=-1.
      SP2N=DTAUQ(NDEP)+0.5*DTAUQ(NDEP)**2
      SP1N=SP1N/SP2N
      SP2N=SP2N/SP2N
C
C  THE PARAMETER DT1 DETERMINES WHEN I HAS TO BE CALCULATED IN A
C  THIN WAY DUE TO NUMERICAL REASONS. THE PARAMETER SHOULD BE SET
C  TO ABOUT THE THIRD ROOT OF THE MACHINE ACCURACY.
C  
      IF(INCRAD.EQ.0) IMINUS(0)=0.0
      T=TAUQ(1)
C
      IF(ITRAN.EQ.3) THEN
        CALL SPLIN0(NDEP,DTAUQ,S,SPRIM,SBISS,STRISS,WW)
      ELSE
        CALL SPLIN1(NDEP,DTAUQ,S,SPRIM,SBISS,STRISS)
      ENDIF
C
C  CALCULATE ONLY NECESSARY EXPONENTIALS, AND DO IT IN A SEPARATE
C  LOOP TO ALLOW THE CRAY TO VECTORIZE THE CALLS
C
      NM=0
      DO 400 K=1,NDEP-1
        IF(DTAUQ(K+1).GT.DT1 .AND. DTAUQ(K+1).LE.DT2) THEN
          NM=NM+1
          WW(NM)=DTAUQ(K+1)
          KP(NM)=K
        ENDIF
  400 CONTINUE
      DO 410 M=1,NM
        WW(M)=EXP(-WW(M))
  410 CONTINUE
      DO 420 M=1,NM
        EXPD(KP(M))=WW(M)
  420 CONTINUE
C
C  CALCULATE OUTGOING INTENSITY IPLUS(K)
C
      AIP(NDEP)=S(NDEP)+SPRIM(NDEP)+SBISS(NDEP)+STRISS(NDEP)
      IPLUS(NDEP)=AIP(NDEP)
      IF(ITRAN.EQ.5) THEN
        IPLUS(NDEP)=0.
        AIP(NDEP)=0.
      ENDIF
      DO 500 K=NDEP-1,1,-1
        IF(DTAUQ(K+1).GT.DT2) THEN
          AIP(K)=S(K)+SPRIM(K)+SBISS(K)+STRISS(K)
        ELSE IF (DTAUQ(K+1).GT.DT1) THEN
          SBISS1=SBISS(K)+DTAUQ(K+1)*STRISS(K)
          AIP(K)=S(K)+SPRIM(K)+SBISS(K)+STRISS(K)+EXPD(K)*
     *     (AIP(K+1)-S(K+1)-SPRIM(K+1)-SBISS1-STRISS(K))
        ELSE
          SBISS1=SBISS(K)+DTAUQ(K+1)*STRISS(K)
          AIP(K)=AIP(K+1)-DTAUQ(K+1)*(AIP(K+1)-S(K+1)-
     *     0.5*DTAUQ(K+1)*(AIP(K+1)-S(K+1)-SPRIM(K+1)-
     *     0.333333333*DTAUQ(K+1)*(AIP(K+1)-S(K+1)-SPRIM(K+1)-
     *     SBISS1-0.25*DTAUQ(K+1)*(IPLUS(K+1)-SPRIM(K+1)-SBISS1-
     *     STRISS(K)))))
        ENDIF
        IPLUS(K)=AIP(K)
  500 CONTINUE
C
C  SURFACE INTENSITY
C
      IPLUS(0)=IPLUS(1)+T*(1.-T*(.5-T*(1./6.-T/24.)))*(S(1)-IPLUS(1))
C
C  CALCULATE INCOMING INTENSITY
C
      IF(T.LT.0.01) THEN
        AIM(1)=IMINUS(0)+T*(1.-T*(.5-T*(1./6.-T/24.)))*
     *          (S(1)-IMINUS(0))
      ELSE IF(T.LT.20.) THEN
        AIM(1)=IMINUS(0)+(1.-EXP(-T))*(S(1)-IMINUS(0))
      ELSE
        AIM(1)=S(1)
      ENDIF
      IMINUS(1)=AIM(1)
      DO 600 K=1,NDEP-1
        IF(DTAUQ(K+1).LE.DT1) THEN
          AIM(K+1)=AIM(K)-DTAUQ(K+1)*(AIM(K)-S(K)-
     *     .5*DTAUQ(K+1)*(AIM(K)-S(K)+SPRIM(K)-.333333333*
     *     DTAUQ(K+1)*(AIM(K)-S(K)+SPRIM(K)-SBISS(K)-.25*
     *     DTAUQ(K+1)*(AIM(K)+SPRIM(K)-SBISS(K)+STRISS(K)))))
        ELSE IF(DTAUQ(K+1).LE.DT2) THEN
          SBISS1=SBISS(K)+DTAUQ(K+1)*STRISS(K)
          AIM(K+1)=S(K+1)-SPRIM(K+1)+SBISS1-STRISS(K)+EXPD(K)*
     *     (AIM(K)-S(K)+SPRIM(K)-SBISS(K)+STRISS(K))
        ELSE
          SBISS1=SBISS(K)+DTAUQ(K+1)*STRISS(K)
          AIM(K+1)=S(K+1)-SPRIM(K+1)+SBISS1-STRISS(K)
        ENDIF
        IMINUS(K+1)=AIM(K+1)
  600 CONTINUE
C
C  FEAUTRIERS P
C
      DO 700 K=1,NDEP
        P(K)=0.5*(IPLUS(K)+IMINUS(K))
  700 CONTINUE
C
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE SPLIN0(N,DX,F,D,D2,D3,WW)
C
C  CALCULATES STANDARD CUBIC SPLINE WITH CONTINOUS SECOND DERIVATIVE
C
C  CODED BY G.SCHARMER, 1982
C
C  DX(K)=X(K)-X(K-1)
C  WW IS A WORKING ARRAY
C
      INCLUDE 'PREC'
      DIMENSION DX(N),F(N)
      DIMENSION D(N),D2(N),D3(N),WW(N)
C
      FAC=-DX(2)/DX(3)
      D(2)=(2.-FAC)*(DX(2)+DX(3))
      C2=DX(3)+FAC*DX(2)
      WW(2)=(F(2)-F(1))/DX(2)
      WW(3)=(F(3)-F(2))/DX(3)
      D3(2)=6.*(WW(3)-WW(2))
      DO 100 K=3,N-1
        D(K)=2.*(DX(K)+DX(K+1))
        WW(K+1)=(F(K+1)-F(K))/DX(K+1)
        FAC=-DX(K)/D(K-1)
        D(K)=D(K)+FAC*DX(K)
        IF(K.EQ.3) D(K)=D(K)+FAC*(C2-DX(K))
        D3(K)=6.*(WW(K+1)-WW(K))+FAC*D3(K-1)
  100 CONTINUE
      FAC=-DX(N)/D(N-2)
      AN=-DX(N-1)-DX(N)+FAC*DX(N-1)
      D3(N)=FAC*D3(N-2)
      FAC=-AN/D(N-1)
      D2(N)=(D3(N)+FAC*D3(N-1))/(DX(N-1)+FAC*DX(N))
      DO 150 K=N-1,3,-1
        D2(K)=(D3(K)-DX(K+1)*D2(K+1))/D(K)
  150 CONTINUE
      D2(2)=(D3(2)-C2*D2(3))/D(2)
      D2(1)=((DX(2)+DX(3))*D2(2)-DX(2)*D2(3))/DX(3)
      DO 180 K=1,N-1
        D(K)=WW(K+1)-(D2(K+1)+D2(K)+D2(K))*DX(K+1)/6.
        D3(K)=(D2(K+1)-D2(K))/DX(K+1)
  180 CONTINUE
      D3(N)=D3(N-1)
      D(N)=D(N-1)+DX(N)*(D2(N-1)+0.5*DX(N)*D3(N-1))
C
      RETURN
      END
C
C**********************************************************************
C
      SUBROUTINE SPLIN1(N,DX,F,D,D2,D3)
C
C  CALCULATES CUBIC SPLINES WITH LOCALLY DETERMINED FIRST DERIVATIVE
C  STABILITY BETTER THAN FOR STANDARD CUBIC SPLINE, AT COST OF
C  DISCONTINOUS SECOND DERIVATIVE.
C
C  CODED BY AA.NORDLUND, MAR-83
C
C  DX(K)=X(K)-X(K-1)
C
      INCLUDE 'PREC'
      DIMENSION DX(N),F(N)
      DIMENSION D(N),D2(N),D3(N)
C
C  FIRST DERIVATIVE BY CENTERED DIFFERENCE
C
      DO 100 K=2,N-1
        D(K)=(F(K+1)-F(K-1))/(DX(K+1)+DX(K))
  100 CONTINUE
      D(1)=(F(2)-F(1))/DX(2)
      D(N)=(F(N)-F(N-1))/DX(N)
C
C  SECOND AND THIRD DERIVATIVE FROM SPLINE CONDITIONS
C
      DO 110 K=1,N-1
        CX=1.0/DX(K+1)
        DFDX=(F(K+1)-F(K))*CX
        D2(K)=(6.*DFDX-4.*D(K)-2.*D(K+1))*CX
        D3(K)=6.*(D(K)+D(K+1)-2.*DFDX)*CX*CX
  110 CONTINUE
      CXN=1.0/DX(N)
      DFDXN=(F(N)-F(N-1))*CXN
      D2(N)=(4.*D(N)+2.*D(N-1)-6.*DFDXN)*CXN
      D3(N)=D3(N-1)
C
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE MATMUL(A,B,C,N1,N2,N3,M1,M2,M3)
C
C  GIVES A#B=C WHERE DIMENSIONS ARE
C  A:   N1 ROWS BY N2 COLUMNS
C  B:   N2 ROWS BY N3 COLUMNS
C  C:   N1 ROWS BY N3 COLUMNS 
C
C  DIMENSIONS OF ARRAYS ARE M1,M2,M3, ACTUAL USED DIMENSIONS N1,N2,N3
C
      INCLUDE 'PREC'
C
      DIMENSION A(M1,M2),B(M2,M3),C(M1,M3)
C
      DO 300 I3=1,N3
        DO 200 I1=1,N1
          C(I1,I3)=0.0
          DO 100 I2=1,N2
            C(I1,I3)=C(I1,I3)+A(I1,I2)*B(I2,I3)
  100     CONTINUE
  200   CONTINUE
  300 CONTINUE
C
      END
C
C****************************************************************
C
      SUBROUTINE MATADD(A,B,C,N1,N2,M1,M2)
C
C  GIVES A+B=C WHERE DIMENSIONS ARE
C  A:   N1 ROWS BY N2 COLUMNS
C  B:   N1 ROWS BY N2 COLUMNS
C  C:   N1 ROWS BY N2 COLUMNS 
C
C  DIMENSIONS OF ARRAYS ARE M1,M2 ACTUAL USED DIMENSIONS N1,N2
C
      INCLUDE 'PREC'
C
      DIMENSION A(M1,M2),B(M1,M2),C(M1,M2)
C
      DO 200 I2=1,N2
        DO 100 I1=1,N1
          C(I1,I2)=A(I1,I2)+B(I1,I2)
  100   CONTINUE
  200 CONTINUE
C
      END
C
C****************************************************************
C
      SUBROUTINE MATSUB(A,B,C,N1,N2,M1,M2)
C
C  GIVES A-B=C WHERE DIMENSIONS ARE
C  A:   N1 ROWS BY N2 COLUMNS
C  B:   N1 ROWS BY N2 COLUMNS
C  C:   N1 ROWS BY N2 COLUMNS 
C
C  DIMENSIONS OF ARRAYS ARE M1,M2 ACTUAL USED DIMENSIONS N1,N2
C
      INCLUDE 'PREC'
C
      DIMENSION A(M1,M2),B(M1,M2),C(M1,M2)
C
      DO 200 I2=1,N2
        DO 100 I1=1,N1
          C(I1,I2)=A(I1,I2)-B(I1,I2)
  100   CONTINUE
  200 CONTINUE
C
      END

C
C****************************************************************
C
      SUBROUTINE MATINV(A,N,NP)
C
C  ADAPTED FROM GAUSSJ OF NUMERICAL RECIPES
C
C: MATINV 07-12-26  MODIFICATIONS: (MATS CARLSSON)
C:        REPLACED PAUSE WITH CALL STOP
C:
      INCLUDE 'PREC'
C
      INTEGER N,NP,NMAX
      DIMENSION A(NP,NP)
      PARAMETER (NMAX=50)
      INTEGER I,ICOL,IROW,J,K,L,LL,INDXC(NMAX),INDXR(NMAX),IPIV(NMAX)
C
      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
        BIG=0.
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (ABS(A(J,K)).GE.BIG)THEN
                  BIG=ABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
                CALL STOP('SINGULAR MATRIX IN MATINV')
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
        IF (A(ICOL,ICOL).EQ.0.) CALL STOP('SINGULAR MATRIX IN MATINV')
        PIVINV=1./A(ICOL,ICOL)
        A(ICOL,ICOL)=1.
        DO 16 L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.
            DO 18 L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
23        CONTINUE
        ENDIF
24    CONTINUE
      RETURN
      END
C  (C) COPR. 1986-92 NUMERICAL RECIPES SOFTWARE VS1&V%1JW#<?4210(9P#.