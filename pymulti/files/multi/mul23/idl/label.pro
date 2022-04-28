pro label,l_sty,text ,x0=x0,y0=y0,thick=thick,psym=psym,dy=dy,size=size,$
                      color=color,dx=dx,normalized=normalized
;+
;   label,l_sty,text [,x0=x0,y0=y0,thick=thick,psym=psym,dy=dy,size=size]
;
;            produces line-style explanation in figure
;            l_sty   integer array with line-styles
;            text    string array with texts
;            x0,y0   optional position for upper right corner of
;                    explanation box in coordinates normalized
;                    to corners of plotting box if /normalized is given the coordinates
;                    are relative to plotting window instead of plotting box
;            thick   integer array with line thickness
;            psym    integer array with psym
;            color   integer array with color
;            size    text size
;-
if (n_params(0) lt 2) then begin
  print,'label,l_sty,text [,x0=x0,y0=y0,thick=thick,psym=psym,dy=dy,size=size]'
  return
endif

n=n_elements(l_sty)

if(n_elements(x0) eq 0) then x0=0.15
if(n_elements(y0) eq 0) then y0=0.8
if(n_elements(thick) eq 0) then thick=intarr(n)+1
if(n_elements(psym) eq 0) then psym=intarr(n)
if(n_elements(dy) eq 0) then dy=0.04
if(n_elements(size) eq 0) then size=1.0
if(n_elements(color) eq 0) then color=intarr(n)+!p.color

if(keyword_set(normalized)) then begin
  x0_n=x0
  y0_n=y0
endif else begin
  x0_n=!x.window(0)+x0*(!x.window(1)-!x.window(0))  ; transform to norm. coor.
  y0_n=!y.window(0)+y0*(!y.window(1)-!y.window(0))  ; transform to norm. coor.
endelse

if(n_elements(dx) eq 0) then dx=0.08
dx0=0.01
dy0=0.01
for i=0,n-1 do begin
  y00=y0_n-i*dy+dy0
  x=x0_n+findgen(3)/2*dx
  y=x-x+y00
  plots,x,y,/normal,linestyle=l_sty(i),thick=thick(i),$
   psym=psym(i),color=color(i)
  xyouts,x0_n+dx+dx0,y00,text(i),/normal,size=size,color=color(i)
endfor

return
end
