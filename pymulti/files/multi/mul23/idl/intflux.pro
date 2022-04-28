;+
; NAME:
;	INTFLUX
;
; PURPOSE:
;	This function calculates the integrated flux or intensity.
;
; CATEGORY:
;	Multi
;
; CALLING SEQUENCE:
; 
;	Result = INTFLUX(Kr [,Mu=Mu,/Nocont])
;
; INPUTS:
;	Kr:	The transition for which the integrated flux/intensity
;		is calculated. 
;
; KEYWORD PARAMETERS:
;	MU:	If mu-index is given, integrated intensity is returned.
;
;	NOCONT:	If /nocont  is given, continuum flux/intensity is not 
;		subtracted.
;
; OUTPUTS:
;	Integrated flux or intensity
;
; COMMON BLOCKS:
;	common_multi
;
; PROCEDURE:
;	Uses Trapez integration
;
; MODIFICATION HISTORY:
; 	Written by:	Mats Carlsson.
;-
function intflux,kr,mu=mu,nocont=nocont

@common_multi

if (n_params(0) eq 0) then begin
  print,'intflux(kr [,mu=mu,/nocont])'
  return,0
endif

if (n_elements(mu) ne 0) then ff=reform(outint(*,mu,*)) else ff=flux
double,kr,ff,xx,yy
x1=cc*1.e8/(xx+alamb(kr))
if(keyword_set(nocont)) then yy=yy else yy = yy-ff(0,kr)
return,trapez(x1,yy)

end
