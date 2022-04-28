pro watom,fileout ,luw=luw
;+
;   watom,fileout [,luw=luw]
;
;            writes multi-format atom file
;            output on fileout
;            if luw is given fileout is assumed open and there is no close
;            collisional data is not written
;
;-
@common_multi

if (n_params(0) lt 1) then begin
  print,'watom,fileout [,luw=luw]'
  return
endif

if(n_elements(luw) eq 0) then append=0 else append=1

; printout

if(append eq 0) then openw,luw,fileout,/get_lun

printf,luw, ATOMID
printf,luw,'*  ABUND    AWGT'
printf,luw, ABND,AWGT/uu,format='(2f8.2)'

ncont=nrad-nline
printf,luw,'*  NK NLIN NCNT NFIX'
printf,luw, NK,NLINE,NCONT,NRFIX,format='(4i5)'
blip="'"
for i=0,nk-1 do begin
  lab=label(i)
  for j=strlen(lab)+1,20 do lab=lab+' '            ; pad label with blanks
  EV2=EV(I)*EE/HH/CC
  printf,luw,ev2,g(i),blip,lab,blip,ion(i),$
   format="(f12.3,f8.2,2x,a1,a20,a1,i4)"
endfor
;
;  BOUND-BOUND TRANSITIONS IN DETAIL
;
if (nline gt 0) then begin
  printf,luw,$
   '* UP LO      F        NQ   QMAX   Q0 IW      GA          GVW        GS'
  for KR=0,NLINE-1 do begin
    j=jrad(kr)-1
    i=irad(kr)-1
    printf,luw, j+1,i+1,f(kr),nq(kr),qmax(kr),q0(kr),$
     iwide(kr),ga(kr),gw(kr),gq(kr),$
     format="(i4,i4,e11.3,i6,f7.1,f5.1,i3,e11.3,f10.1,e11.3)"
    IF(QMAX(KR) LT 0.0) OR (Q0(KR) LT 0.0) THEN begin
      printf,luw, q(0:nq(kr)-1,kr)
    endif
  endfor
endif
;
;  BOUND-FREE TRANSITIONS IN DETAIL
;
if (nrad gt nline) then begin
  printf,luw,'* UP LO     A0        NQ   QMAX'
  for KR=NLINE,NRAD-1 do begin
    printf,luw, jrad(kr),irad(kr),f(kr),nq(kr),qmax(kr),$
     format="(i4,i4,e11.3,i6,f7.1)"
    IF(QMAX(KR) LT 0.0) THEN begin
      tab=fltarr(2,nq(kr))
      tab(0,*)=q(0:nq(kr)-1,kr)
      tab(1,*)=alfac(0:nq(kr)-1,kr-nline)
      printf,luw, tab,format='(f13.3,e13.5)'
    ENDIF
  endfor
endif
;
;  FIXED TRANSITIONS PARAMETERS
;
if (nrfix gt 0) then begin
  printf,luw,$
   '* UP LO IPHO         A0       TRAD ITRAD'
  for KF=0,NRFIX-1 do begin
    printf,luw, jfx(kf),ifx(kf),ipho(KF),a0(KF),trad(KF),itrad(KF),$
     format='(i4,i4,i5,e11.3,f11.0,i6)'
    IF(ITRAD(KF) EQ 4) THEN begin
      print,'itrad=4 NOT IMPLEMENTED'
      return
    endif
  endfor
endif

if(append eq 0) then free_lun,luw

RETURN
END



