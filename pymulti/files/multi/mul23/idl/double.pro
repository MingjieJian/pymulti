;+
; NAME:
;	DOUBLE
;
; PURPOSE:
;	This procedure makes profile symmetric around x(0). Used for 
;	fluxes and intensities. For two-sided profiles, original profile 
;	is returned.
;
; CATEGORY:
;	Multi
;
; CALLING SEQUENCE:
;
;	DOUBLE, KR, Y, XX,YY
;
; INPUTS:
;	KR:	transition number
;
;	Y:	array to be made symmetric, either OUTINT or FLUX
;
; OUTPUTS:
;	XX:	delta lambda in Angstrom
;
;	YY:	symmetric Y values as function of wavelength
;
; COMMON BLOCKS:
;	common_multi
;
; EXAMPLE:
;	Typical call sequence:
;	double,0,flux,xx,yy
;	plot,xx,yy
;
; MODIFICATION HISTORY:
; 	Written by:	Mats Carlsson
;-
pro double,kr,y,xx,yy

@common_multi

if n_params(0) eq 0 then begin
  print,'double,kr,y,xx,yy'
  return
endif

ndum=size(y)
if ndum(0) eq 2 then begin  ; two dimensional y array: flux values in y
  y2=y(1:nq(kr),kr)
endif else begin
  y2=y(1:nq(kr),nmu-1,kr)
endelse
qn=qnorm*1.e5/cc

if(ind(kr) eq 2) then begin
  yy=y2
  qx=q(0:nq(kr)-1,kr)                  ; extract q
endif else begin
  yy=fltarr(nq(kr)*2-1)                ; arrays for symmetrizized variables
  xx=fltarr(nq(kr)*2-1)
  qx=fltarr(nq(kr)*2-1)
  q2=q(0:nq(kr)-1,kr)                  ; extract q
  qx(nq(kr)-1)=q2
  qx(0)=-reverse(q2)                   ; make symmetric q
  yy(nq(kr)-1)=y2
  yy(0)=reverse(y2)
endelse

xx=-alamb(kr)*qx*qn/(1.0+qx*qn)

end
