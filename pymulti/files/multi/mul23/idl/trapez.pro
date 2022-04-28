;+
; NAME:
;	TRAPEZ
;
; PURPOSE:
;	This function performs trapezoidal integration.
;
; CATEGORY:
;	Multi
;
; CALLING SEQUENCE:
; 
;	Result = TRAPEZ(X, Y)
;
; INPUTS:
;	X:	X-array
;
;	Y:	Y-array
;
; OUTPUTS:
;	This function returns the integral Y*dx
;
; MODIFICATION HISTORY:
; 	Written by:	Mats Carlsson
;-
function trapez,x0,y0

n=n_elements(x0)
x=fltarr(n)+x0
y=fltarr(n)+y0
integrand=(y+shift(y,-1))*(shift(x,-1)-x)
return,total(integrand(0:n-2))*0.5

end
