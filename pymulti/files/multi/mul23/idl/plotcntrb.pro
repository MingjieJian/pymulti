;+
; NAME:
;	PLOTCNTRB
;
; PURPOSE:
;	This procedure plots contribution function
;
; CATEGORY:
;	Multi
;
; CALLING SEQUENCE:
;
;	PLOTCNTRB, Kr, Xscale, Cntrb
;
; INPUTS:
;	Kr:	Transition number. 
;
;	Xscale:	X scale on plot (i.e. taulg).
;
;	Cntrb:	Contribution function to be plotted.
;
; COMMON BLOCKS:
;	common_multi
;
; SIDE EFFECTS:
;	Plots in current window
;
; EXAMPLE:
;	Typical calls are: plotcntrb,0,alog10(tau),cntrbi
;                          plotcntrb,0,alog10(tau),cntrbr
;
; MODIFICATION HISTORY:
; 	Written by:	Mats Carlsson, March 1988.
;-
pro plotcntrb,kr,xscale,cntrb,yrange=yrange,_extra=e

@common_multi

if n_params(0) eq 0 then begin
  print,'plotcntrb,kr,xscale,cntrb'
  return
endif

ymin=min(cntrb(*,*,kr))                ;find minimum and
ymax=max(cntrb(*,*,kr))                ;maximum values for plot scale
if(n_elements(yrange) eq 0) then yrange=[ymin,ymax]
plot,xscale,cntrb(*,0,kr),yrange=yrange,_extra=e     ;plot first frequency
for ny=1,nq(kr)-1 do oplot,xscale,cntrb(*,ny,kr) ;plot all other frequencies

return
end
