pro remove_xl,x1,y1,x2,y2 ,maxdx=maxdx,maxdy=maxdy,maxerr=maxerr,$
 npoint=npoint,verbose=verbose,si30=si30,tref=tref,iref=iref,xlref=xlref,osmet=osmet
;+
;   remove_xl,x1,y1,x2,y2
;
;            remove wavelength points that minimizes the change
;            in trapez(1/x,y*bny(x,t)) where t goes from 5000K to 40000K
;            t is set 20% higher shortward of 912 Å
;            point is not removed if dx becomes larger than maxdx or
;            error gets larger than maxerr
;
;            if osmet is set: does not remove points in
;            range [1500,20000]
;
;-
if(n_params() lt 4) then begin
  print,'remove_xl,x1,y1,x2,y2 ,maxdx=maxdx,maxdy=maxdy,maxerr=maxerr,$'
  print,' npoint=npoint,/verbose,/si28,/si30,tref=tref,iref=iref,xlref=xlref,/osmet'
  return
endif

if(n_elements(maxerr) eq 0) then maxerr=0.01
if(n_elements(maxdx) eq 0) then maxdx=0.01
if(n_elements(maxdy) eq 0) then maxdy=5.0
if(n_elements(npoint) eq 0) then npoint=1

if(n_elements(si30) ne 0) then begin
  edges=[911.8,1100.16,1100.93,1200.,1238.81,1239.78,$
   1514.49,1516.26,1519.62,1521.10,1522.89,1526.29,$
   1674.20,1682.29,$
   1768,1974.94,1986.20,2070.,2513.5,3647.0,3756,8205.8]
endif else begin
  edges=[911.8,1100.16,1100.93,1200.,1238.81,1239.78,$
   1517.93,1524.58,$
   1674.20,1682.29,$
   1768,1974.94,1986.20,2070.,2513.5,3647.0,3756,8205.8]
endelse
x=x1
y=y1

; reference integral

n=n_elements(x)
if(n_elements(tref) eq 0) then begin
  tref=[5000.,10000.,15000.,20000.,30000.,40000.]
endif
if(n_elements(iref) eq 0) then begin
  nt=n_elements(tref)
  int_ref=fltarr(nt)
  for it=0,nt-1 do begin
    t=fltarr(n)+tref(it)
    iw=where(x lt 912.,count)
    if(count gt 0) then t(iw)=t(iw)*1.2
    int_ref(it)=trapez(1./x,y*planck(x,t))
  endfor
endif else begin
  nt=1
  intep,xlref,iref,x,inu
  int_ref=fltarr(nt)
  int_ref[0]=trapez(1./x,y*inu)
endelse

; remove loop

loop:
  n=n_elements(x)
  err_cum=fltarr(n)
  for it=0,nt-1 do begin
    t=tref(it)
    err=fltarr(n)
    if(n_elements(iref) ne 0) then int0=trapez(1./x,y*inu)
    for i=1,n-1 do begin             ; remove point i
      case i of
        0:   iw=indgen(n-1)+1
        n-1: iw=indgen(n-1)
        else:iw=[indgen(i),indgen(n-i-1)+i+1]
      endcase
      if((i eq 0) or (i eq n-1) or (x[i] gt 20000.) or $
        (x[i] lt 1500.) or (n_elements(osmet) eq 0)) then begin
        if(n_elements(iref) eq 0) then begin
          tt=fltarr(n-1)+t
          iw2=where(x(iw) lt 912,count)
          if(count gt 0) then tt(iw2)=t*1.2
          int=trapez(1./x(iw),y(iw)*planck(x(iw),tt))
        endif else begin
          case i of
            0:   int=int0-trapez(1./x[0:1],y[0:1]*inu[0:1])
            n-1: int=int0-trapez(1./x[n-2:n-1],y[n-2:n-1]*inu[n-2:n-1])
            else:int=int0-trapez(1./x[i-1:i+1],y[i-1:i+1]*inu[i-1:i+1])+$
             trapez([1./x[i-1],1./x[i+1]],[y[i-1]*inu[i-1],y[i+1]*inu[i+1]])
          endcase
        endelse
        err(i)=abs(int/int_ref(it)-1.)
      endif else begin
        err[i]=4.0                                               ; OS
      endelse
      if(i gt 0) and (i lt n-1) then $
       if(abs(y(i)/(y(i-1)<y(i))-1.) gt maxdy) then err(i)=1.0   ; fast change
      if(i gt 0) and (i lt n-1) then $
       if(abs(y(i)/(y(i+1)<y(i))-1.) gt maxdy) then err(i)=1.0   ; fast change
      if((i gt 0) and (i lt n-1)) then $
       if(abs(x(i+1)-x(i-1))/x[i] gt maxdx) then err(i)=2.0 
      if(n_elements(iref) eq 0) then if(min(abs(x(i)-edges)) lt 1.1) then err(i)=3.0  ; don't remove edges
    endfor
    err_cum=err > err_cum
  endfor

; remove point minimizing error

  minerr=min(err_cum(1:n-1),imin)
  i=imin+1
  if(minerr gt maxerr) then goto,finish
  case i of
    0:   iw=indgen(n-1)+1
    n-1: iw=indgen(n-1)
    else:iw=[indgen(i),indgen(n-i-1)+i+1]
  endcase
  x=x(iw)
  y=y(iw)
  if(n_elements(iref) ne 0) then intep,xlref,iref,x,inu
  err=err(iw)
  errmax=minerr
  if(n-1 le npoint) then goto,finish
  if(n_elements(verbose) ne 0) then begin
    if(verbose ne 0) then begin
      plot,x1,y1,lin=2
      oplot,x,y,psym=-1
      iw=where(err eq 1.,count)
      if(count gt 0) then oplot,x(iw),y(iw),psym=4
      iw=where(err eq 2.,count)
      if(count gt 0) then oplot,x(iw),y(iw),psym=5
      iw=where(err eq 3.,count)
      if(count gt 0) then oplot,x(iw),y(iw),psym=2
    endif
  endif
goto,loop

finish:
x2=x
y2=y
if(n_elements(verbose) ne 0) then begin
  if(verbose ne 0) then begin
    print,'max error=',errmax,' n=',n-1
    label,[0,0,0],['maxdx','maxdy','edge'],psym=[5,4,2],x0=0.8,y0=0.9
;    text=''
;    read,text
;    plot,err_cum,yran=[1.e-5,10],/ylog
;    oplot,[i,i],[err_cum[i],err_cum[i]],psym=2
;    oplot,[0,n-1],[1,1]*maxerr,lin=2
  endif
endif

end
