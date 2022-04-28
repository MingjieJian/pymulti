;+
; NAME:
;	COOL_SUM
;
; PURPOSE:
;	This procedure adds up all cooling contributions
;	and puts the sum into cool_total
;
; CATEGORY:
;	Multi
;
; CALLING SEQUENCE:
;
;	COOL_SUM, Cool, Cool_total
;
; INPUTS:
;	Cool:		Cool(k,kr) cooling function at depth k, transition kr
;
; OUTPUTS:
;	Cool_total:	Cool_total(k) total cooling function at depth k
;
; EXAMPLE:
;	Typical call: cool_sum,cool,cool_total
;
; MODIFICATION HISTORY:
; 	Written by:	Mats Carlsson
;-
pro cool_sum,cool,cool_total

if n_params(0) eq 0 then begin
  print,'cool_sum,cool,cool_total'
  return
endif

dum=size(cool)
ndep=dum(1)
nrad=dum(2)
cool_total=fltarr(ndep)
for kr=0,nrad-1 do begin
  cool_total=cool_total+cool(*,kr)
endfor

return
end
