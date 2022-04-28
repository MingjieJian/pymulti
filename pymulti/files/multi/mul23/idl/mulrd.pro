;+
; NAME:
;	MULRD
;
; PURPOSE:
;	Reads data from file. If no file name is given, 
;	the name idl1.def_ext is assumed where def_ext is set with
;	procedure default. After execution of this
;	procedure, most common block variables are accessible.
;
; CATEGORY:
;	Multi.
;
; CALLING SEQUENCE:
;
;	MULRD, File
;
; INPUTS:
;	File:	Input file containing all the multi-data.
;               Default file name is IDL1 if def_ext='none', else
;               idl1.def_ext
;
; OUTPUTS:
;	Fills most common blocks
;
; COMMON BLOCKS:
;	common_multi
;
; MODIFICATION HISTORY:
;	Written by:	Mats Carlsson, March 1988.	
;       040916  added automatic detection of endianness 
;-
pro mulrd,file

@common_multi

if n_params(0) eq 0 then begin
  if(n_elements(def_ext) eq 0) then begin
    message,'you need to set default extension with procedure default first',/info
    return
  endif
  if(def_ext eq 'none') then file='IDL1' else file='idl1.'+def_ext
endif
dum=file_search(file,count=count)
if(count eq 0) then begin
  message,file+' not found',/info
  return
endif

hh=6.626176E-27
cc=2.99792458E+10
if n_elements(lu1) ne 0 then free_lun,lu1
get_lun,lu1
openr,lu1,file,/f77_unformatted
rec=assoc(lu1,lonarr(1))
test=rec[0]
free_lun,lu1
if(test eq 32) then swap_endian=0 else swap_endian=1
openr,lu1,file,/f77_unformatted,swap_endian=swap_endian

ndep=0l
nk=0l
nline=0l
nwide=0l
nrad=0l
nrfix=0l
nmu=0l
mq=0l

iformat=0L          ; standard format
print,'reading values from '+file
forrd,lu1, ndep,nk,nline,nwide,nrad,nrfix,nmu,mq
if(ndep eq 0) then begin
  forrd,lu1, iformat
  if(iformat eq 2) then begin
    nsel=0L
    forrd,lu1, nsel
    krsel=lonarr(nsel)
    forrd,lu1, krsel
  endif
  forrd,lu1, ndep,nk,nline,nwide,nrad,nrfix,nmu,mq
endif

if(nrad gt 0) then begin
  nq=lonarr(nrad)
  forrd,lu1, nq
endif

ev=fltarr(nk)
g=fltarr(nk)
ion=lonarr(nk)
if(nrad gt 0) then begin
  ktrans=lonarr(nrad)
  jrad=lonarr(nrad)
  irad=lonarr(nrad)
  f=fltarr(nrad)
  iwide=lonarr(nrad)
  ga=fltarr(nrad)
  gw=fltarr(nrad)
  gq=fltarr(nrad)
endif
krad=lonarr(nk,nk)
z=fltarr(ndep)
gij=fltarr(ndep)
if (nwide gt 0) then begin
  alfac=fltarr(mq,nwide)
  frq=fltarr(mq+1,nwide)
endif
if(nrad gt 0) then alamb=fltarr(nrad)
if(nline gt 0) then a=fltarr(nline)
b=fltarr(nk,nk)
totn=fltarr(ndep)
if(nrad gt 0) then bp=fltarr(ndep,nrad)
nstar=fltarr(nk,ndep)
n=fltarr(nk,ndep)
c=fltarr(nk,nk,ndep)
if nrfix gt 0 then begin
  jfx=lonarr(nrfix)
  ifx=lonarr(nrfix)
  ipho=lonarr(nrfix)
  a0=fltarr(nrfix)
  trad=fltarr(nrfix)
  itrad=lonarr(nrfix)
endif
dnyd=fltarr(ndep)
if(nline gt 0) then adamp=fltarr(ndep,nline)
label=strarr(nk)
for i=0,nk-1 do label[i]=string(0,'(i20)')
atomid=label(0)
crout=string(0,'(i6)')
cmass=fltarr(ndep)
temp=fltarr(ndep)
nne=fltarr(ndep)
vel=fltarr(ndep)
tau=fltarr(ndep)
taulg=fltarr(ndep)
xnorm=fltarr(ndep)
height=fltarr(ndep)
atmoid=string(0,'(i72)')
dpid=atmoid
dptype=' '
vturb=fltarr(ndep)
bh=fltarr(5,ndep)
nh=fltarr(6,ndep)
rho=fltarr(ndep)
if(nrad gt 0) then begin
  qmax=fltarr(nrad)
  q0=fltarr(nrad)
  ind=lonarr(nrad)
  q=fltarr(mq,nrad)
  wq=fltarr(mq,nrad)
  wphi=fltarr(ndep,nrad)
  sl=fltarr(ndep,nrad)
endif
sbck=fltarr(ndep)
s=fltarr(ndep)
rny=fltarr(ndep)
if(nline gt 0) then begin
  weqlte=fltarr(nline)
  weq=fltarr(nline)
endif
if(nrad gt 0) then begin
  rij=fltarr(ndep,nrad)
  rji=fltarr(ndep,nrad)
  flux=fltarr(mq+1,nrad)
  outint=fltarr(mq+1,nmu,nrad)
  cool=fltarr(ndep,nrad)
endif
xmu=fltarr(nmu)
wmu=fltarr(nmu)

forrd,lu1, qnorm
forrd,lu1, abnd,awgt
forrd,lu1, ev
forrd,lu1, g
forrd,lu1, ion
forrd,lu1, hn3c2
if(nrad gt 0) then begin
  forrd,lu1, ktrans
  forrd,lu1, jrad
  forrd,lu1, irad
  forrd,lu1, f
  forrd,lu1, iwide
  forrd,lu1, ga
  forrd,lu1, gw
  forrd,lu1, gq
endif
if(iformat eq 0) then begin
  forrd,lu1, krad
endif else begin
  for kr=0,nrad-1 do begin
    i=irad(kr)-1
    j=jrad(kr)-1
    krad(i,j)=kr+1
    krad(j,i)=kr+1
  endfor
endelse
if(iformat eq 0) then forrd,lu1, z
if (nwide gt 0) then forrd,lu1, alfac
forrd,lu1, hny4p
if(nrad gt 0) then forrd,lu1, alamb
if(nline gt 0) then forrd,lu1, a
if(iformat eq 0) then begin
  forrd,lu1, b
endif else begin
  for kr=0,nline-1 do begin $
    i=irad(kr)-1
    j=jrad(kr)-1
    hn3c2_line=2.*hh*cc/(alamb(kr)*1.e-8)^3
    b(j,i)=a(kr)/hn3c2_line
    b(i,j)=g(j)/g(i)*b(j,i)
  endfor
endelse

forrd,lu1, totn
if(nrad gt 0) and (iformat eq 0) then begin
  forrd,lu1, bp
endif
forrd,lu1, nstar
forrd,lu1, n
if(iformat eq 0) then forrd,lu1, c
if nrfix gt 0 then begin
  forrd,lu1, jfx
  forrd,lu1, ifx
  forrd,lu1, ipho
  forrd,lu1, a0
  forrd,lu1, trad
  forrd,lu1, itrad
endif
forrd,lu1, dnyd
if(nline gt 0) and (iformat eq 0) then forrd,lu1, adamp
forrd,lu1, label
forrd,lu1, atomid
forrd,lu1, crout
;
;  common block catmos
;
forrd,lu1, grav
forrd,lu1, cmass
forrd,lu1, temp
forrd,lu1, nne
forrd,lu1, vel
forrd,lu1, tau
taulg=alog10(tau)
forrd,lu1, xnorm
forrd,lu1, height
forrd,lu1, atmoid,dpid,dptype
;
;  common block catmo2
;
forrd,lu1, vturb
if(iformat eq 0) then forrd,lu1, bh
forrd,lu1, nh
if(iformat eq 0) then forrd,lu1, rho
;
;  common block csline
;
if(nrad gt 0) then begin
  forrd,lu1, qmax
  forrd,lu1, q0
  forrd,lu1, ind
  forrd,lu1, diff
  forrd,lu1, q
  if(iformat eq 0) then forrd,lu1, wq
endif
forrd,lu1, wqmu
if (nwide gt 0) then forrd,lu1, frq
if(nrad gt 0) then begin
  if(iformat eq 0) then forrd,lu1, wphi
  forrd,lu1, sl
endif
if(nline gt 0) then begin
  forrd,lu1, weqlte
  forrd,lu1, weq
endif
if(nrad gt 0) then begin
  if(iformat eq 0) then forrd,lu1, rij
  if(iformat eq 0) then forrd,lu1, rji
  if(iformat ne 2) then begin
    forrd,lu1, flux
  endif else begin
    dum=fltarr(mq+1)
    for kr0=0,nsel-1 do begin
      forrd,lu1, dum
      flux(*,krsel(kr0)-1)=dum
    endfor
  endelse
  if(iformat ne 2) then begin
    forrd,lu1, outint
  endif else begin
    dum=fltarr(mq+1,nmu)
    for kr0=0,nsel-1 do begin
      forrd,lu1, dum
      outint(*,*,krsel(kr0)-1)=dum
    endfor
  endelse
  if(iformat eq 0) then forrd,lu1, cool
endif
;
;  common block cgausi
;
forrd,lu1, xmu
forrd,lu1, wmu
;
;  common block cconst
;
forrd,lu1, ee,hh,cc,bk,em,uu,hce,hc2,hck,ek,pi
free_lun,lu1

if(iformat ne 0) then  begin
  for kr=0,nline-1 do begin $
    bp(*,kr)=planck(alamb(kr),temp)
  endfor
endif

return
end

@common_multi
end
