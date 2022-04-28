pro watmos,atmos_file,dscale_file,rstrt_file
;+
;   watmos,atmos_file,dscale_file,rstrt_file
;
;            writes multi input files from idl variables
;-
@common_multi

if (n_params(0) lt 1) then begin
  message,'watmos,atmos_file ,dscale_file,rstrt_file',/info
  return
endif

dptype=strupcase(dptype)

tab1=fltarr(5,ndep)

openw,luatm,atmos_file,/get_lun
printf,luatm,strtrim(atmoid,2)              ; file heading

if (dptype eq 'M') then begin
  printf,luatm,'MASS SCALE'
  tab1(0,0)=transpose(alog10(cmass))
endif else if (dptype eq 'T') then begin
  printf,luatm,'TAU(5000) SCALE'
  tab1(0,0)=transpose(alog10(tau))
endif else begin
  printf,luatm,'Height scale'
  tab1(0,0)=transpose(height)
  if(n_elements(radyn) ne 0) then tab1(0,0)=transpose(height)*1.e5
endelse
printf,luatm,alog10(grav),format='(f7.4)'
printf,luatm,ndep

tab1(1,0)=transpose(temp)
tab1(2,0)=transpose(nne)
tab1(3,0)=transpose(vel)*qnorm
tab1(4,0)=transpose(vturb)*qnorm
tab2=nh

if(dptype ne 'H') and (n_elements(radyn) eq 0) then begin
;  printf,luatm,tab1,format='(1x,f16.8,f13.2,e13.5,2f13.5)'
  printf,luatm,tab1
endif else begin
;  printf,luatm,tab1,format='(1x,e16.8,f13.2,e13.5,2f13.5)'
  printf,luatm,tab1
endelse
;printf,luatm,tab2,format='(1x,6e13.5)'
printf,luatm,tab2,format='(6e12.4)'
free_lun,luatm


if(n_elements(dscale_file) ne 0) then begin
  openw,ludsc,dscale_file,/get_lun
  printf,ludsc,strtrim(atmoid,2)              ; file heading
  if (dptype eq 'M') then begin
    printf,ludsc,' MASS SCALE'
    printf,ludsc,ndep,alog10(tau(0))
  endif else if (dptype eq 'T') then begin
    printf,ludsc,' TAU(5000) SCALE'
    printf,ludsc,ndep,alog10(cmass(0))
  endif else begin
    printf,ludsc,' HEIGHT SCALE'
    printf,ludsc,ndep,alog10(tau(0))
  endelse
  if(dptype ne 'H') and (n_elements(radyn) eq 0) then begin
    printf,ludsc,tab1(0,*),format='(f16.8)'
  endif else begin
    printf,ludsc,tab1(0,*),format='(e16.8)'
  endelse
  free_lun,ludsc
endif

if(n_elements(rstrt_file) ne 0) then begin
  openw,lurstrt,rstrt_file,/get_lun
  printf,lurstrt,strtrim(atmoid,2)              ; file heading
  if (dptype eq 'M') then begin
    printf,lurstrt,' MASS SCALE'
  endif else if (dptype eq 'T') then begin
    printf,lurstrt,' TAU(5000) SCALE'
  endif else begin
    printf,lurstrt,' HEIGHT SCALE'
  endelse
  printf,lurstrt,alog10(grav)
  printf,lurstrt,ndep
  if(dptype ne 'H') and (n_elements(radyn) eq 0) then begin
    printf,lurstrt,tab1,format='(1x,f16.8,f13.2,e13.5,2f13.5)'
  endif else begin
    printf,lurstrt,tab1,format='(1x,e13.5,f13.2,e13.5,2f13.5)'
  endelse
  for k=0,ndep-1 do begin
    printf,lurstrt,n[*,k],$
   format='(1x,6e13.5)'
  endfor
  free_lun,lurstrt
endif

end
