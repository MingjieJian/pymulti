pro scale2,x,x2,ypos0, x2step,csize ,dy=dy,color=color,x2text=x2text
;+
;   scale2,x,x2,ypos0, [,x2step,csize,dy=dy,color=color,x2text=x2text]
;
;            plots a secondary scale x2 at screen pos ypos0 (0-1)
;            inputs:
;              x      plot x-scale eg. cmass, 
;              x2     secondary scale, eg. alog10(tauq)
;              x2step step between labels on secondary axis
;                     if given as array: values where to put secondary ticks
;              csize  character size
;              color  color for secondary scale
;              x2text text for secondary labels if given as x2step array
;                     defaults to strtrim(x2step,2)
;
;-
;COMMON CATMOS,GRAV,CMASS,TEMP,NNE,VEL,TAU,XCONT,SCAT,XNORM,HEIGHT
;COMMON CATMID,ATMOID,DPID,DPTYPE

if (n_params(0) lt 3) then begin 
  print,'scale2,x,x2,ypos0 [,x2step,csize,dy=dy,color=color,x2text=x2text]'
  return
endif

if(n_elements(dy) eq 0) then dy=0.015  ; length of tick marks in normalized
                                       ; coordinates
if(!p.multi(2) ne 0) then dy=dy/!p.multi(2)

ypos=!y.window(0)+ypos0*(!y.window(1)-!y.window(0))  ; transform to norm. coor.

mint=min(x2)
maxt=max(x2)
!c=0
minlog=fix(mint)
if (minlog gt 0) then minlog=minlog+1
maxlog=fix(maxt)
if (maxlog lt 0) then maxlog=maxlog-1
if (n_elements(x2step) eq 0) then x2step=1
if (n_elements(csize) eq 0) then csize=1.
if (n_elements(color) eq 0) then color=!p.color

dyc=!d.y_ch_size*csize/!d.y_vsize      ; character size in normalized coor.

!psym=0
linetype=!p.linestyle
!p.linestyle=0
if(n_elements(x2step) eq 1) then begin
  for i=minlog,maxlog,x2step do begin
    intep,x2,x,float(i),iposd       ; iposd = tick position in data coordinates
    if(iposd ge min(!x.crange)) and (iposd le max(!x.crange)) then begin
      if(!x.type eq 1) then iposd=alog10(iposd)
      ipos=!x.s(0)+!x.s(1)*iposd   ; ipos = tick position in normalized coor.
      if (n_elements(left) eq 0) then left=ipos         ; full extent of scale
      plots,[ipos,ipos],[ypos,ypos-dy],/norm,linestyle=0,color=color  ; plot tick mark  
      if(csize gt 0) then xyouts,ipos,ypos-dy-dyc*1.1,strtrim(string(i),2),$
       size=csize,/norm,alignment=0.5,color=color
    endif
  endfor
endif else begin
  if(n_elements(x2text) ne n_elements(x2step)) then begin
    x2text=strtrim(x2step,2)
  endif
  for i=0,n_elements(x2step)-1 do begin
    intep,x2,x,x2step(i),iposd   ; iposd = tick position in data coordinates
    if(iposd ge min(!x.crange)) and (iposd le max(!x.crange)) then begin
      if(!x.type eq 1) then iposd=alog10(iposd)
      ipos=!x.s(0)+!x.s(1)*iposd ; ipos = tick position in normalized coor.
      if (n_elements(left) eq 0) then left=ipos         ; full extent of scale
      plots,[ipos,ipos],[ypos,ypos-dy],/norm,linestyle=0,color=color ; plot tick mark  
      if(csize gt 0) then xyouts,ipos,ypos-dy-dyc*1.1,x2text(i),$
       size=csize,/norm,alignment=0.5,color=color
    endif
  endfor
endelse
right=ipos
!p.linestyle=linetype
plots,[left,right],[ypos,ypos],/norm,color=color

return
end
