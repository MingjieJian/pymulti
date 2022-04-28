;+
; NAME:
;	OPCPLOT
;
; PURPOSE:
;	Plots opacity contributions as function of depth.
;
; CATEGORY:
;	Multi
;
; CALLING SEQUENCE:
;
;	OPCPLOT, Xscale, Il, Min_cont [, /TOTAL] 
;
; INPUTS:
;	Xscale:		X scale on plot (i.e. taulg).
;
;	Il:		Wavelength index.
;
;	Min_cont:	Minimum contribution to be plotted (default 0.02)
;
; KEYWORD PARAMETERS:
;	TOTAL:	Contributions relative to total opacity and not relative 
;		to background opacity.
;
; COMMON BLOCKS:
;	common_multi
;
; SIDE EFFECTS:
;	Plots in current window
;
; EXAMPLE:
;	opcrd		reads opacity file
;	print,xla	find wavelength number of interest
;	opcplot,taulg,4,/total
;
; MODIFICATION HISTORY:
; 	Written by:	Mats Carlsson
;-
pro opcplot,xscale,il ,min_cont,total=total

@common_multi

if(n_params(0) lt 2) then begin
  print,'opcplot,xscale,il [,min_cont,/total]'
  return
endif

if(n_elements(min_cont) eq 0) then min_cont=0.02
if(n_elements(total) eq 0) then total=0

iw=indgen(ndep)

y=xscale(iw)*0.+1.                 ; set up y array to get plot limits right
plot,xscale(iw),y,title=strtrim(atmoid,2)+'   !4k!3='+strtrim(string(xla(il)),2),$
 ytitle='Opacity contribution',charsize=1.5,/nodata
linestyle=0
for i=0,19 do begin
  if(total) and (il ne 0) then begin
    nyrd,kr_xla(il),ny_xla(il)
    yref=x*xnorm/rho
  endif else yref=chi(*,il)
  y=prov(i,iw,il)/yref(iw)
  if(max(y) gt min_cont) then begin
    oplot,xscale(iw),y,linestyle=linestyle
    linestyle=linestyle+1
    if(linestyle gt 5) then linestyle=0
    if(linestyle eq 1) then linestyle=2
    y0=max(y)
    x0=xscale(iw(!c))
    xyouts,x0,y0,provid(i),size=1.5
  endif
endfor

return
end
