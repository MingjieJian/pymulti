pro calc_sotbfi,xl,nq,inu,sotbfi ,verbose=verbose
;+
;   calc_sotbfi,xl,nq,inu,sotbfi ,/verbose
;
;            calculate Hinode SOT BFI intensities
;            wavelength array in xl[0:nq[kr]-1,kr] (Angstrom)
;            intensities in inu[0:nq[kr]-1,kr] (erg cm-2 ster-1 s-1 Hz-1)
;            output in sotbfi[i],i=0,5 (photons s-1)
;            pixelsize is 0.053 arcseconds and 15 microns
;
;-
common ccalc_sotbfi,nl,wl,tl

If(n_params() lt 4) then begin
  print,'calc_sotbfi,xl,nq,inu,sotbfi ,/verbose'
  return
endif


cc=2.99792d+10
hh=6.62618d-27
Siz=size(xl)
if(siz[0] eq 1) then nrad=1 else nrad=siz[2]
if(n_elements(nl) eq 0) then restore,'../idl/bfi.idlsave',/verb
nu=cc*1.d8/xl    ; frequencies
Const=(0.053/3600.*!pi/180.)^2*(!pi*50.^2)    ; Solid angle*area of telescope

; go through the SOT filters

sotbfi=fltarr(6)
for i=0,5 do begin
  for kr=0,nrad-1 do begin
    iw=where((xl[0:nq[kr]-1,kr] ge min(wl[0:nl[i]-1,i])*10) and $
             (xl[0:nq[kr]-1,kr] le max(wl[0:nl[i]-1,i])*10),count) ; find line
    if(count gt 0) then begin
      if(keyword_set(verbose)) then begin
        print,'match, i,kr,count=',i,kr,count
      endif
      tf=interpol(tl[0:nl[i]-1,i],wl[0:nl[i]-1,i]*10,xl[iw,kr])
      sotbfi[i]=trapez(nu[iw,kr],tf*inu[iw,kr]/hh/nu[iw,kr])
    endif
  endfor
endfor
sotbfi=sotbfi*const

end
