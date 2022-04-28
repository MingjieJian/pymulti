;+
; NAME:
;	DLAMB
;
; PURPOSE:
;	This function calculates and returns delta lambda for a given 
;	Q array.
;
; CATEGORY:
;	Multi
;
; CALLING SEQUENCE:
; 
;	Result = DLAMB(QQ, Lambda)
;
; INPUTS:
;	QQ:	Frequency parameter in typical Doppler units
;
;	Lambda:	Central wavelength in Angstrom
;
; OUTPUTS:
;	Delta lambda from line center in Angstrom
;
; COMMON BLOCKS:
;	common_multi
;
; MODIFICATION HISTORY:
; 	Written by:	Mats Carlsson
;-
function dlamb,qq,lambda

@common_multi

qn=qnorm*1.e5/cc
return,-lambda*qq*qn/(1.0+qq*qn)
end
