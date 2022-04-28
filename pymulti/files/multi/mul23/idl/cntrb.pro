;+
; NAME:
;	CNTRB
;
; PURPOSE:
;	This procedure reads contribution functions from file and places 
;	contribution function to intensity in cntrbi(k,ny,kr), to relative 
;	absorption in cntrbr(k,ny,kr) and average height of formation in 
;	xmeani(ny,kr) and xmeanr(ny,kr).
;
; CATEGORY:
;	Multi
;
; CALLING SEQUENCE:
;
;	CNTRB, File
;
; INPUTS:
;	File:	File containing contribution function data. 
;               Default file name is IDLCNT if def_ext='none', else
;               idlcnt.def_ext
;
; OUTPUTS:
;	in common block:
;	cntrbi(k,ny,kr)  contribution function to intensity
;	cntrbr(k,ny,kr)  contribution function to relative intensity
;	cntrbf(k,ny,kr)  contribution function to flux
;       xmeani(ny,kr)    average taulg of formation for intensity
;	xmeanr(ny,kr)    average taulg of formation for relative intensity
;	xmeanf(ny,kr)    average taulg of formation for flux
;
;       to screen:
;       'reading contribution functions for kr=',kr
;
; COMMON BLOCKS:
;	common_multi
;
; PROCEDURE:
;	cntrbr as defined by Magain
;
; MODIFICATION HISTORY:
; 	Written by:	Mats Carlsson.
;       04-09-16        added swap endian
;-
pro cntrb,file

@common_multi

if n_params(0) eq 0 then begin
  if(def_ext eq 'none') then file='IDLCNT' else file='idlcnt.'+def_ext
endif

get_lun,lcnt
openr,lcnt,file,/f77_unformatted,swap_endian=swap_endian
ndep=0l
nline=0l
nrad=0l
mq=0l
forrd,lcnt,ndep,nline,nrad,mq
nq=lonarr(nrad)
forrd,lcnt,nq
cntrbi=fltarr(ndep,mq,nrad)              ;declare arrays
cntrbf=fltarr(ndep,mq,nrad)              ;declare arrays
if (nline gt 0) then cntrbr=fltarr(ndep,mq,nline)             ;declare arrays
xmeani=fltarr(mq,nrad)
xmeanf=fltarr(mq,nrad)
if (nline gt 0) then xmeanr=fltarr(mq,nline)
dummy=fltarr(ndep)

for kr=0,nrad-1 do begin
  print,'reading contribution functions for kr=',kr
  for ny=0,nq(kr)-1 do begin
    forrd,lcnt,dummy
    cntrbi(0,ny,kr)=dummy
    xmeani(ny,kr)=trapez(taulg,taulg*dummy)/trapez(taulg,dummy)
    forrd,lcnt,dummy
    cntrbf(0,ny,kr)=dummy
    xmeanf(ny,kr)=trapez(taulg,taulg*dummy)/trapez(taulg,dummy)
  endfor
  if kr le nline-1 then begin
    for ny=0,nq(kr)-1 do begin
      forrd,lcnt,dummy
      cntrbr(0,ny,kr)=dummy
      xmeanr(ny,kr)=trapez(taulg,taulg*dummy)/trapez(taulg,dummy)
    endfor
  endif
endfor

free_lun,lcnt
return
end
