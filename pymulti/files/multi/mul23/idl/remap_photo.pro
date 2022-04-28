pro remap_photo,npoint=npoint,maxerr=maxerr,maxdx=maxdx,osmet=osmet,$
  verbose=verbose,si30=si30,edge=edge,tref=tref,iref=iref,xlref=xlref
;+
;   remap_photo
;
;            remap photo-ionization cross-sections to a given wavelength set
;
;-
@common_multi
common copctab_xl,xl0

if(n_elements(maxerr) eq 0) then maxerr=0.001
if(n_elements(npoint) eq 0) then npoint=50
if(n_elements(maxdx) eq 0) then maxdx=0.01

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
dx=50.0

if(n_elements(xlref) eq 0) then begin
  opctab_xl      ; find grid used by opctab
  xl0=reverse(xl0)
endif else begin
  iw=sort(xlref)
  xl0=reverse(xlref[iw])
endelse

xl=cc*1.e8/frq

openw,luw,'remap_photo.tmp',/get_lun

if(n_elements(verbose) ne 0) then begin
  window,0
  window,2
endif

for kr=nline,nrad-1 do begin
  kt=kr-nline
  x1=xl(1:nq(kr),kt)
  y1=alfac(0:nq(kr)-1,kt)
  ipol_int,x1,y1,x2,y2                   ; remap to dx=10 preserving integral
  if(n_elements(edge) ne 0) then begin
    xl2=[xl0,xl(0)]                      ; add edge wavelength
    xl2=reverse(xl2(sort(xl2)))
  endif else begin
    xl2=reverse(xl0(sort(xl0)))
  endelse
  iw=where((xl2 ge min(x2)) and (xl2 le max(x2)))
  x=xl2(iw)
  intep,x2,y2,x,y                        ; interpolate to opctab grid
  y=y > 0.                               ; prevent negative values
  remove_xl,x,y,x2,y2,npoint=npoint,maxerr=maxerr,maxdx=maxdx,$
   verbose=verbose,si30=si30,tref=tref,iref=iref,xlref=xlref,osmet=osmet
  nq2=n_elements(x2)
  printf,luw,label[irad[kr]-1],label[jrad[kr]-1],format="('*',a,' -- ',a)"
  printf,luw,'* up  lo alphac    nq qmax  q0'
  printf,luw,jrad(kr),irad(kr),f(kr),nq2,$
   format="(2i4,e9.2,i4,' -1.0 0.0')"
  for ny=0,nq2-1 do begin
    printf,luw,x2(ny),y2(ny),format="(f13.2,e10.2)"
  endfor
  if(n_elements(verbose) ne 0) then wset,2-!d.window
endfor
free_lun,luw

print,'new photoionization grid written to remap_photo.tmp'
end

