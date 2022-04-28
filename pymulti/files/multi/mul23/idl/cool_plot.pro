;+
; NAME:
;	COOL_PLOT
;
; PURPOSE:
;	This procedure plots the cooling function. If line type is given, 
;	overlay plot.
;
; CATEGORY:
;	Multi
;
; CALLING SEQUENCE:
;
;	COOL_PLOT, X, Cool [,Line]
;
; INPUTS:
;	X:	Atmospheric height/depth scale to be used for x-axis
;
;	Cool:	Cooling function.
;
; OPTIONAL INPUTS:
;	Line:	Line type, gives overlay plot
;	
;
; COMMON BLOCKS:
;	CCOOL_PLOT:	saving local variables for overlay plot
;
; SIDE EFFECTS:
;	Plots in current window
;
; EXAMPLE:
;	Typical call sequences:
;	cool_plot,height,cool_total
;       for kr=0,nrad-1 do cool_plot,height,cool(*,kr),kr+1
;
; MODIFICATION HISTORY:
; 	Written by:	Mats Carlsson
;-
pro cool_plot,x,cool,line

common ccool_plot,yc_min,yc_max,yh_min,yh_max ; saves local variables
                                                                                
if n_params(0) eq 0 then begin
  print,'cool_plot,x,cool [,line]'
  return
endif
                                                                                
npar=n_params(0)
dum=size(x)
ndep=dum(1)
yc=fltarr(ndep)                     ; declare arrays
yh=fltarr(ndep)
ic=where(cool gt 0.0)               ; find where cooling
ih=where(cool lt 0.0)               ; find where heating
yc(ic)=alog10(cool(ic))             ; set logarithmic cooling value
yc(ih)=-100.                        ; set heating points outside cooling plot
yh(ih)=alog10(-cool(ih))            ; repeat for heating
yh(ic)=-100.
if npar le 2 then begin
  yc_max=float(fix(max(yc(ic))))    ; find maximum cooling and set plot limit
  if yc_max gt 0. then yc_max=yc_max+1.
  yc_min=yc_max-5.0                 ; set plot range to 5 orders of magnitude
  yh_max=float(fix(max(yh(ih))))
  if yh_max gt 0. then yh_max=yh_max+1.
  yh_min=yh_max-5.0
endif
                                                                                
if max(x) gt 10. then begin         ; turn height scale around
  x_range=[max(x),min(x)]
endif else begin
  x_range=[min(x),max(x)]
endelse  				
                                                                                
if npar eq 3 then begin
  linetype=line
  no_eras=1
endif else begin
  linetype=0
  no_eras=0
endelse
if npar le 2 then erase
plot,x,yh,xrange=x_range,yrange=[yh_max,yh_min],noeras=no_eras,$ ; plot heating
  linestyle=linetype,position=[0.1,0.1,0.9,0.45]
plot,x,yc,xrange=x_range,yrange=[yc_min,yc_max],/noeras,$        ; plot cooling
  linestyle=linetype,position=[0.1,0.5,0.9,0.95]
return
end
