;+
; NAME:
;	TRADB
;
; PURPOSE:
;	This function calculates trad from given (i,lambda) array. Trad is 
;	the radiation temperature.
;
; CATEGORY:
;	Multi
;
; CALLING SEQUENCE:
; 
;	Result = trad(I, Lambda)
;
; INPUTS:
;	I:	Intensity I_nu in cgs units
;
;	Lambda:	Wavelength in Angstrom
;
; OUTPUTS:
;	This function returns trad
;
; MODIFICATION HISTORY:
; 	Written by:	Mats Carlsson, April 1988.
;-
function tradb,i,lambda

c=2.99792458e10
h=6.626176e-27
k=1.380662e-16
l=lambda*1.e-8
trd=h*c/k/l/alog(2.*h*c/i/(l*l*l)+1.)
return,trd

end
