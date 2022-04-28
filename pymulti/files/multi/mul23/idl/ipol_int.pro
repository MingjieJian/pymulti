pro ipol_int,x1,y1,x2,y2
;+
;   ipol_int,x1,y1,x2,y2
;
;            interpolate x1,y1 to x2,y2 preserving trapez(1./x1,y1)
;            x2 is generated with dx=10
;
;-
if(n_params() lt 4) then begin
  message,'ipol_int,x1,y1,x2,y2',/continue,/informational
  return
endif

; sort x in ascending order

iw=sort(x1)
x=x1(iw)
y=y1(iw)
n=n_elements(x)

; add new points with spacing 10

nxnew=fix(x(n-1)/10.)-fix(x(0)/10.)+2
xnew=findgen(nxnew)*10.+fix(x(0)/10.)*10.+5.

; only include inner points

iw=where((xnew gt x(0)) and (xnew lt x(n-1)))
xnew=xnew(iw)

; interpolate values at midpoints

intep,x,y,xnew,ynew

; add endpoints

xnew=[x(0),xnew,x(n-1)]
ynew=[y(0),ynew,y(n-1)]

; integrate in each interval

nxnew=n_elements(xnew)
x2=fltarr(nxnew+1)
y2=x2

for i=0,nxnew-2 do begin
  iw=where((x gt xnew(i)) and (x lt xnew(i+1)),count)
  if(count gt 0) then begin
    xx=1./[xnew(i),x(iw),xnew(i+1)]
    yy=[ynew(i),y(iw),ynew(i+1)]
  endif else begin
    xx=1./[xnew(i),xnew(i+1)]
    yy=[ynew(i),ynew(i+1)]
  endelse
  int=trapez(xx,yy)
  x2(i+1)=1./(0.5*(1./xnew(i+1)+1./xnew(i)))
  y2(i+1)=int/(1./xnew(i+1)-1./xnew(i))
endfor

; add end-points

x2(0)=x(0)
y2(0)=y(0)
x2(nxnew)=x(n-1)
y2(nxnew)=y(n-1)

; reverse if input is descending

if(x1(n-1) lt x1(0)) then begin
  x2=reverse(x2)
  y2=reverse(y2)
endif

end
