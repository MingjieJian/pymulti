pro calc_sotnfi_na,xl,inu,nab,nar,na0 ,dx0=dx0,verbose=verbose,xrange=xrange
;+
;   calc_sotnfi_na,xl,inu,nab,nar,na0 ,dx0=dx0,/verbose
;
;            calculate Hinode SOT NFI intensities
;            wavelength array in xl (Angstrom) in vacuum
;            intensities in inu[0:nq-1,nx,ny] (erg cm-2 ster-1 s-1 Hz-1)
;            output in nab, nar, na0 (photons s-1) for blue wing, red
;            wing and line-center
;            pixelsize is 0.0796 arcseconds and 15 microns
;
;-
common ccalc_sotnfi_na,nl,wl,tl

If(n_params() lt 5) then begin
  print,'calc_sotnfi_na,xl,inu,nab,nar,na0 ,dx0,/verbose'
  return
endif

if(n_elements(dx0) eq 0) then dx0=0.0

cc=2.99792d+10
hh=6.62618d-27
siz=size(xl)
nq=siz[1]
xl0=589.756          ; line-center wavelength in nm in vacuum
dxl=(xl*0.1-xl0)     ; wavelengths in nm relative to line center
siz=size(inu)
if(siz[1] ne nq) then begin
  message,'nq not the same for xl and inu',/info
  return
endif
nx=1
ny=1
if(siz[0] ge 2) then nx=siz[2]
if(siz[0] eq 3) then ny=siz[3] 
restore,'../idl/tf_doppler_crude.idlsave'
nu=cc*1.d7/(wnad+xl0)    ; frequencies
Const=(0.053/3600.*!pi/180.)^2*(!pi*50.^2)    ; Solid angle*area of telescope

; go through x and y

nab=fltarr(nx,ny)
nar=fltarr(nx,ny)
na0=fltarr(nx,ny)
if(keyword_set(verbose)) then mtimer,0,ny,'calc_sotnfi_na',/start
for j=0,ny-1 do begin
  for i=0,nx-1 do begin
    i0=interpol(inu[*,i,j],dxl,wnad+dx0*0.1)   ; interpolate intensity
    na0[i,j]=trapez(nu,tfmean*i0/hh/nu)  
    nar[i,j]=trapez(nu,tfred*i0/hh/nu)  
    nab[i,j]=trapez(nu,tfblue*i0/hh/nu)
    if(keyword_set(verbose) and i+j eq 0) then begin
      text=''
      spec_ct
      if(n_elements(xrange) eq 0) then xrange=[-200,200]
      plot,wnad*1.e3,i0/max(i0),xrange=xrange,$
       ytitle='I/I!dc!n',col=0,back=255,xstyle=4,ymargin=[4,4]
      axis,color=0,xtitle='!4Dk!3 [pm]'
      axis,xaxis=1,xrange=!x.crange/xl0*3e2,xstyle=1,xtitle='Velocity [km/s]',color=0
      inorm=!y.crange[1]*0.6/max([tfred,tfblue,tfmean])
      oplot,wnad*1.e3,tfred*inorm,col=1
      oplot,wnad*1.e3,tfmean*inorm,col=100
      oplot,wnad*1.e3,tfblue*inorm,col=200
      loadct,0
      read,'q to quit: ',text
      if(strtrim(text,2) eq 'q') then return
    endif
    if(keyword_set(verbose)) then mtimer,'calc_sotnfi_na',j,ny,/remain
  endfor
endfor
if(keyword_set(verbose)) then mtimer,0,ny,'calc_sotnfi_na',/start
na0=na0*const
nab=nab*const
nar=nar*const

end
