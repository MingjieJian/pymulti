;+
; NAME:
;	OPCRD
;
; PURPOSE:
;	This procedure reads opacity data from file.
;
; CATEGORY:
;	Multi
;
; CALLING SEQUENCE:
;
;	OPCRD, File
;
; INPUTS:
;	File:	Name of file containing opacity data.
;               Default file name is IDLOPC if def_ext='none', else
;               idlopc.def_ext
;
; OUTPUTS:
;	In common block
;
; COMMON BLOCKS:
;	common_multi
;
; MODIFICATION HISTORY:
; 	Written by:	Mats Carlsson
;       04-09-16        added swap endian
;-
pro opcrd,file

@common_multi

if n_params(0) eq 0 then begin
  if(def_ext eq 'none') then file='IDLOPC' else file='idlopc.'+def_ext
endif

get_lun,lu1
dum=file_search(file,count=count)
if(count eq 0) then begin
  message,file+' not found, remember to set IDLOPC=1 in INPUT',/info
  return
endif
openr,lu1,file,/f77_unformatted,swap_endian=swap_endian

ncont=nrad-nline                            ; find number of wavelengths
nfreq=nline+1                               ; with background opacities
for kr=nline,nrad-1 do nfreq=nfreq+nq[kr]

xla=fltarr(nfreq)
chi=fltarr(ndep,nfreq)
sumabs=fltarr(ndep,nfreq)
sumsca=fltarr(ndep,nfreq)
prov=fltarr(20,ndep,nfreq)
provid0=strarr(20)
provid0[00]='H bf '
provid0[01]='H ff '
provid0[02]='H-bf '
provid0[03]='H-ff '
provid0[04]='H2-  '
provid0[05]='H2+  '
provid0[06]='H+H  '
provid0[07]='He-  '
provid0[08]='He   '
provid0[09]='Si   '
provid0[10]='Mg   '
provid0[11]='Al   '
provid0[12]='Fe   '
provid0[13]='C    '
provid0[14]='N    '
provid0[15]='O    '
provid0[16]='elec '
provid0[17]='R[H) '
provid0[18]='R(H2)'
provid0[19]='R(He)'

provid=strarr(20)
provid[00]='H!dbf!n'
provid[01]='H!dff!n'
provid[02]='H!u-!dbf!n'
provid[03]='H!u-!dff!n'
provid[04]='H!d2!u-!n'
provid[05]='H!d2!u+!n'
provid[06]='H+H'
provid[07]='He!u-!n'
provid[08]='He   '
provid[09]='Si   '
provid[10]='Mg   '
provid[11]='Al   '
provid[12]='Fe   '
provid[13]='C    '
provid[14]='N    '
provid[15]='O    '
provid[16]='!4r!3(e!u-!n)'
provid[17]='R(H) '
provid[18]='R(H!d2!n)'
provid[19]='R(He)'

dum=fltarr(23)
ntp=0l

for kr=0,nfreq-1 do begin
  forrd,lu1,dum1
  xla[kr]=dum1
  for k=0,ndep-1 do begin
    forrd,lu1,ntp,dum
    chi[k,kr]=dum[0]
    prov[0,k,kr]=dum[1:20]
    sumabs[k,kr]=dum[21]
    sumsca[k,kr]=dum[22]
  endfor
endfor
                                ; set up mapping between transition
nla=n_elements(xla)             ; numbers, frequencies and opacity numbers
kr_xla=intarr(nla)
ny_xla=intarr(nla)
for kr=0,nline-1 do begin
  kr_xla[kr+1]=kr
  ny_xla[kr+1]=0
endfor
ii=nline
for kr=nline,nrad-1 do begin
  for ny=0,nq[kr]-1 do begin
    ii=ii+1
    kr_xla[ii]=kr
    ny_xla[ii]=ny
  endfor
endfor

free_lun,lu1
return
end
