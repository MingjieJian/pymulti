;+
; NAME:
;	GACALC
;
; PURPOSE:
;	This procedure calculates ga values from line list
;
; CATEGORY:
;	Multi
;
; CALLING SEQUENCE:
;
;	GACALC
;
; INPUTS:
;	From common
;
; OUTPUTS:
;	In common:
;	GA:	GA(kr) Summed A values for all transitions from upper level
;		and lower level of transition kr
;
; COMMON BLOCKS:
;	common_multi
;
; MODIFICATION HISTORY:
; 	Written by:	Mats Carlsson.
;-
pro gacalc,dummy

@common_multi

asum_lev=fltarr(nk)                 ; inverse lifetime of level
for kr=0,nline-1 do begin           ; add up inverse lifetimes
  j=jrad(kr)-1
  asum_lev(j)=asum_lev(j)+a(kr)
endfor

for kr=0,nline-1 do begin           ; set ga to gamma(i)+gamma(j)
  i=irad(kr)-1
  j=jrad(kr)-1
  ga(kr)=asum_lev(i)+asum_lev(j)
endfor

return
end
