;+
; NAME:
;	CONVL
;
; PURPOSE:
;	This function converts vacuum wavelengths to air for wavelengths 
;	greater than 2000 Angstroms. 
;
; CATEGORY:
;	Multi
;
; CALLING SEQUENCE:
; 
;	Result = CONVL(Lambda)
;
; INPUTS:
;	Lambda:	Vacuum wavelength
;
; OUTPUTS:
;	This function returns the air wavelength for wavelengths
;       greater than 2000 Angstrom, else the vacuum wavelength
;
; PROCEDURE:
;	Algorithm from Starlink program IUEDR
;	You might not need this section for your routine.
;
; MODIFICATION HISTORY:
; 	Written by:	Mats Carlsson
;-
function convl,lambda

array=size(lambda)
if array(0) eq 0 then begin
  if (lambda lt 2000.) then $
    lambda_ut = lambda $
  else $
    lambda_ut = lambda/(1.0+2.735182e-4+131.4182/lambda/lambda+ $
    2.76249e8/lambda/lambda/lambda/lambda)
  return,lambda_ut
endif else begin
  lambda_ut=lambda
  if(max(lambda) gt 2000.) then begin
    iw=where (lambda gt 2000.)
    lambda_ut[iw] = lambda[iw]/(1.0+2.735182e-4+ $
     131.4182/lambda[iw]/lambda[iw]+ $
     2.76249e8/lambda[iw]/lambda[iw]/lambda[iw]/lambda[iw])
  endif
  return,lambda_ut
endelse
end
