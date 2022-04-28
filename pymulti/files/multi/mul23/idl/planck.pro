;+
; NAME:
;	PLANCK
;
; PURPOSE:
;	This function calculates the Planck function, B_ny(lambda,t).
;
; CATEGORY:
;	Multi
;
; CALLING SEQUENCE:
;
;	Result = PLANCK(Lambda, T)
;
; INPUTS:
;	Lambda:	Wavelength in Angstroms.
;
;	T:	Temperature in Kelvins.
;
;		either input (but not both) can be an array
;
; OUTPUTS:
;	Returns B_ny in cgs units
;
; MODIFICATION HISTORY:
; 	Written by:	Mats Carlsson, April 1988.
;-
function planck,lambda,t

c=2.99792458e10
h=6.626176e-27
k=1.380662e-16
l=lambda*1.e-8

return, 2.*h*c/l/l/l/(exp(h*c/l/k/t)-1.)
end
