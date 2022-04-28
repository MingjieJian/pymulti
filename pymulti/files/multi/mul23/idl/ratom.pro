;+
; NAME:
;	RATOM
;
; PURPOSE:
;	Reads atomic data from Multi atom file without having to run
;	Multi first.
;
; CATEGORY:
;	Multi.
;
; CALLING SEQUENCE:
;
;	RATOM, File
;
; INPUTS:
;	File:	Atomic data input file
;
; OUTPUTS:
;	Fills common blocks
;
; COMMON BLOCKS:
;	common_multi
;
; MODIFICATION HISTORY:
;	Written by:	Mats Carlsson, March 1988.	
;-
pro ratom,file ,qnorm0=qnorm0
;
@common_multi                                                                                
if (n_params(0) eq 0) then begin
  message,'ratom,file ,qnorm=qnorm',/info
  return
endif          

;
; strip comment lines from input file
;
dum=file_search(file,count=count)
if(count eq 0) then begin
  message,file+' not found',/info
  return
endif
text=' '

openr,lu1,file,/get_lun
openw,lu2,'dums.dat',/get_lun
while not eof(lu1) do begin
  readf,lu1,text
  if (strmid(text,0,1) ne '*') then begin
    k0=strpos(text,"'")
    if (k0 ne -1) then begin
      k1=strpos(text,"'",k0+1)
      k2=strlen(text)
      printf,lu2,strmid(text,0,k0)           ; text up to label
      printf,lu2,strmid(text,k0+1,k1-k0-1)   ; label
      printf,lu2,strmid(text,k1+1,k2-k1)     ; text after label
    endif else begin
      printf,lu2,text
    endelse
  endif
endwhile
free_lun,lu1
close,lu2
openr,lu2,'dums.dat'

if(n_elements(qnorm0) eq 0) then read,'Give QNORM: ',qnorm else qnorm=qnorm0

EE=1.602189E-12
HH=6.626176E-27
CC=2.99792458E10
EM=9.109534E-28
UU=1.6605655E-24
BK=1.380662E-16
PI=3.14159265359
HCE=HH*CC/EE*1.E8
HC2=2.*HH*CC *1.E24
HCK=HH*CC/BK*1.E8
EK=EE/BK
HNY4P=HH*CC/QNORM/4./PI*1.E-5
atomid=''

readf,lu2, ATOMID
atomid=strtrim(atomid,2)

READF,LU2, ABND,AWGT
AWGT=AWGT*UU

nk=0
nline=0
ncont=0
nrfix=0
READF,LU2, NK,NLINE,NCONT,NRFIX
NRAD=NLINE+NCONT
ev=dblarr(nk)
g=fltarr(nk)
label=strarr(nk)
ion=intarr(nk)
krad=intarr(nk,nk)

lab=''
ev2=0.d0
for I=0,NK-1 do begin
  READF,LU2, ev2,g2
  readf,lu2, lab
  readf,lu2, ion2
  G[I]=g2
  LABEL[I]=lab
  ION[I]=ion2
  EV[I]=EV2*CC *HH/EE
  for J=0,NK-1 do begin
    KRAD[I,J]=0.
  endfor
endfor
;C
;C  BOUND-BOUND TRANSITIONS IN DETAIL
;C  CALCULATE LAMBDA, A AND B
;C  IF QMAX OR Q0.LT.0 FREQUENCY POINTS IN DOPPLER UNITS ARE READ
;C
KT=-1
ktrm=-1
mq=4000
if (nrad gt 0) then begin
  f=fltarr(nrad)
  nq=fltarr(nrad)
  qmax=fltarr(nrad)
  q0=fltarr(nrad)
  iwide=intarr(nrad)
  ga=fltarr(nrad)
  gw=fltarr(nrad)
  gq=fltarr(nrad)
  irad=intarr(nrad)
  jrad=intarr(nrad)
  alamb=fltarr(nrad)
  a=fltarr(nrad)
  b=fltarr(nk,nk)
  ktrans=intarr(nrad)
  q=fltarr(mq,nrad)
  mtrm=20
  mterm=5
  if(nline gt 0) then begin
    nterm=intarr(nline)
    kterm=intarr(mterm,nline)
  endif
  determ=fltarr(mtrm)
  wterm=fltarr(mtrm)
  fterm=fltarr(mtrm)
  gaterm=fltarr(mtrm)
  gwterm=fltarr(mtrm)
  gqterm=fltarr(mtrm)
endif                   
if (nline gt 0) then begin
  for KR=0,NLINE-1 do begin
    readf,lu2, j,i,f2,nq2,qmax2,q02,io2,ga2,gw2,gq2
    i=i-1
    j=j-1
    f[kr]=f2
    nq[kr]=nq2
    qmax[kr]=qmax2
    q0[kr]=q02
    iwide[kr]=io2
    nterm[kr]=io2
    ga[kr]=ga2   
    gw[kr]=gw2
    gq[kr]=gq2
    IF(QMAX[KR] LT 0.0) OR (Q0[KR] LT 0.0) THEN begin
      dum=fltarr(nq2)
      READF,LU2, dum
      q[0:nq2-1,kr]=dum
    ENDIF
    if(nterm(kr) ge 2) then begin
      for itrm=0,nterm[kr]-1 do begin
        readf,lu2,de2,w2,f2,ga2,gw2,gq2
        ktrm=ktrm+1
        determ[ktrm]=de2
        wterm[ktrm]=w2
        fterm[ktrm]=f2
        gaterm[ktrm]=ga2
        gwterm[ktrm]=gw2
        gqterm[ktrm]=gq2
      endfor
    endif
    KRAD[I,J]=KR+1
    KRAD[J,I]=KR+1
    IRAD[KR]=I+1
    JRAD[KR]=J+1
    ALAMB[KR]=HCE/(EV[J]-EV[I])
    A[KR]=F[KR]*6.671E15*G[I]/G[J]/ALAMB[KR]/ALAMB[KR]
    B[J,I]=ALAMB(KR)^3/HC2*A(KR)
    B[I,J]=G[J]/G[I]*B[J,I]
  endfor
endif
;C
;C  BOUND-FREE TRANSITIONS IN DETAIL
;C  IF QMAX.LT.0.0 FREQUENCY POINTS IN ANGSTROM (STARTING WITH
;C  THRESHOLD AND DECREASING) AND CROSSECTIONS IN CM2 ARE READ.
;C  UNIT CONVERSION IN  ROUTINE  FREQC
;C
mq=4000
if (nrad gt nline) then begin
  alfac=fltarr(mq,nrad-nline)
  frq=fltarr(mq+1,nrad-nline)
  for KR=NLINE,NRAD-1 do begin
    READF,LU2, J,I,f2,nq2,qmax2
    i=i-1
    j=j-1
    kt=kt+1
    F[KR]=f2
    NQ[KR]=nq2
    QMAX[KR]=qmax2
    IF(QMAX[KR] LT 0.0) THEN begin
      tab=fltarr(2,nq2)
      READF,LU2, tab
      tab=transpose(tab)
      wledge=hh*cc/ee/(ev[j]-ev[i])*1.e8         ; wavelength at edge
      frq[0,kt]=cc*1.e8/wledge
      q[0:nq2-1,kr]=tab[*,0]                     ; wavelength table
      frq[1:nq2,kt]=cc*1.e8/q[0:nq2-1,kr]      ; frequencies
      alfac[0:nq2-1,kt]=tab[*,1]
    ENDIF else begin
      wledge=hh*cc/ee/(ev[j]-ev[i])*1.e8         ; wavelength at edge
      frq[0,kt]=cc*1.e8/wledge                   ; frequency at edge
      frq[nq[kr],kt]=cc*1.e8/qmax[kr]        ; frequency at shortest wavelength
      frq[0,kt]=findgen(nq[kr]+1)/nq[kr]*(frq[nq[kr],kt]-frq[0,kt])+frq[0,kt]
      q[0:nq[kr]-1,kr]=cc*1.e8/frq[1:nq[kr],kt]
      alfac[0:nq[kr]-1,kt]=alfac[0,kt]*(frq[1,kt]/frq[1:nq[kr],kt])^3
    endelse
    KTRANS[KR]=KR-NLINE+KT
    IRAD[KR]=I+1
    JRAD[KR]=J+1
    KRAD[I,J]=KR
    KRAD[J,I]=KR
    GA[KR]=0.
    GW[KR]=0.
    GQ[KR]=0.
    ALAMB[KR]=HCE/(EV[J]-EV[I])
  endfor
  mq=max(nq)+1
  alfac=alfac[0:mq,*]
  frq=frq[0:mq,*]
endif
;C
;C  INPUT FIXED TRANSITIONS PARAMETERS
;C
if (nrfix gt 0) then begin
  jfx=intarr(nrfix)
  ifx=intarr(nrfix)
  ipho=intarr(nrfix)
  a0=fltarr(nrfix)
  trad=fltarr(nrfix)
  itrad=intarr(nrfix)
  for KF=0,NRFIX-1 do begin
    readf,lu2, jfx2,ifx2,ipho2,a02,trad2,itrad2
    JFX[KF]=jfx2
    IFX[KF]=ifx2
    IPHO[KF]=ipho2
    A0[KF]=a02
;    TRAD[KF]=trad2
    ITRAD[KF]=itrad2
    IF(ITRAD[KF] EQ 4) THEN begin
      IF(IPHO[KF] EQ 0) THEN begin
        readf,lu2,j,i,ff2,nqf2,qmaxf2,q0f2,io2,gaf2,gwf2,gqf2
        IF(qmaxf2 LT 0.0) OR (q0f2 LT 0.0) THEN begin
          tab=fltarr(nqf2)
          READF,LU2, tab
        ENDIF
      endif ELSE begin
        READF,LU2, J,I,FF2,NQF2,QMAXF2
        IF(QMAXF2 LT 0.0) THEN begin
          tab=fltarr(2*nqf2)
          READF,LU2, tab
        ENDIF
      ENDelse
    endif
  endfor
endif
free_lun,lu2

q=q[0:max(nq),*]

RETURN
END
