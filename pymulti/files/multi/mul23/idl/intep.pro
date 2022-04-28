;+
; NAME:
;	INTEP
;
; PURPOSE:
;	This procedure is used for interpolation
;       uses a Hermite spline interpolation that often avoids over-shoot
;	found with cubic splines
;
; CATEGORY:
;	Multi
;
; CALLING SEQUENCE:
;
;	INTEP,X,Y,Xpin,Ypout [, Nowarning=Nowarning]
;
; INPUTS:
;	X:	X-array to interpolate in
;
;	Y:	Y-array to interpolate in
;
;	Xpin:	X-array where Y-values are wanted
;
; KEYWORD PARAMETERS:
;	NOWARNING:	If /NOWARNING there will be no warning messages
;			if extrapolation occurs
;
; OUTPUTS:
;	Ypout:	Y-array with returned interpolated values
;
; OPTIONAL OUTPUTS:
;	If extrapolation is attempted, there is a warning message
;
; RESTRICTIONS:
;	Instead of extrapolation, the nearest end value is used
;
; PROCEDURE:
;	ref: publications of the dominion astrophysical observatory, 
;	xvi,6,67 graham hill: intep, an effective interpolation subroutine
;
; MODIFICATION HISTORY:
; 	Written by:	Mats Carlsson
;-
pro intep,x,y,xpin,ypout ,nowarning=nowarning

if (n_params(0) ne 4) then begin
  print,'intep,x,y,xp,yp'
  return
endif
if(n_elements(nowarning) eq 0) then nowarning=0
n=n_elements(x)
dum=size(xpin)
if(dum(0) eq 0) then begin
  np=1
  xp=fltarr(1)
  yp=fltarr(1)
  xp(0)=xpin
endif else begin
  np=dum(1)
  xp=xpin
  yp=fltarr(np)
endelse

      ier=1
      io=0
      iup=0
      if(x(1) lt x(0)) then iup=1
      n1=n-2
      for it=0,np-1 do begin
        if ((xp(it) ge x(n-1)) and (iup eq 0)) or $
           ((xp(it) le x(n-1)) and (iup eq 1)) then begin
   g5:    yp(it)=y(n-1)
          goto,g6
        endif else begin
          if ((xp(it) le x(0)) and (iup eq 0)) or $
             ((xp(it) ge x(0)) and (iup eq 1)) then begin
            yp(it)=y(0)
   g6:      if (xp(it) ne x(0)) and (xp(it) ne x(n-1)) then ier=2
            goto,loop
          endif
        endelse
   g8:
        for i=io,n-1 do begin
          if(xp(it) lt x(i)) and (iup eq 0) then goto,g2
          if(xp(it) gt x(i)) and (iup eq 1) then goto,g2
        endfor
        goto,g5
  g2:   i=i-1
        if(i eq io-1) then goto,g4
        io=i+1
        lp1=1./(x(i)-x(i+1))
        lp2=1./(x(i+1)-x(i))
        if(i eq 0) then fp1=(y(1)-y(0))/(x(1)-x(0))
        if(i eq 0) then goto,g3
        fp1=(y(i+1)-y(i-1))/(x(i+1)-x(i-1))
   g3:  if(i ge n1) then fp2=(y(n-1)-y(n-2))/(x(n-1)-x(n-2))
        if(i ge n1) then goto,g4
        fp2=(y(i+2)-y(i))/(x(i+2)-x(i))
   g4:  xpi1=xp(it)-x(i+1)
        xpi=xp(it)-x(i)
        l1=xpi1*lp1
        l2=xpi*lp2
        yp(it)=y(i)*(1.-2.*lp1*xpi)*l1*l1+y(i+1)*(1.-2.*lp2*xpi1)* $
         l2*l2+fp2*xpi1*l2*l2+fp1*xpi*l1*l1
        goto,loop
      loop:
      endfor
result:
if(dum(0) eq 0) then begin
  ypout=yp(0)
endif else begin
  ypout=yp
endelse
if (ier eq 2) and (nowarning eq 0) then $
  print,' warning: xp outside range, yp set to end-point value'
return
end
