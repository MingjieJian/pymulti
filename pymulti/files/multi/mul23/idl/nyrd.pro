;+
; NAME:
;	NYRD
;
; PURPOSE:
;	This procedure reads ny dependent variables from files 
;	file_idlny and file_jny for line kr, frequency ny. if no 
;	filename is given, the defaults IDLNY and JNY are used. Variables 
;	are: pms, iplus, iminus, p, s, tauq, dtauq, xcont, sc, scat, x, 
;	jny, sbck, rny.
;
; CATEGORY:
;	Multi
;
; CALLING SEQUENCE:
;
;	NYRD, Kr, Ny [, File_idlny, File_jny, MDEP=MDEP, /DP, /JNYDP]
;
; INPUTS:
;	Kr:		The transition number for which variables are read
;
;	Ny:		Frequency number for which variables are read
;
;	File_idlny:	idlny file name. 
;			Default file name is IDLNY if def_ext='none', else
;			idlny.def_ext
;
;	File_jny:	jny file name
;			Default file name is JNY if def_ext='none', else
;			jny.def_ext
;
;	
; KEYWORD PARAMETERS:
;
;	MDEP:		Dimension MDEP. Defaults to NDEP. If dimension is
;			different from NDEP this keyword has to be given
;			in order to get correct JNY
;
;	DP:		signals that program was compiled in double precision
;			this means that the record length for odd values of
;			NDEP is different than in SP which this keyword
;			takes care of. Note that the file JNY is still
;                       assumed to be in single precision (see JNYDP keyword)
;
;       JNYDP:          JNY file is assumed to be in double precision.
;                       sets /DP as well
;
;
; OUTPUTS:
;	In common:
;	pms		P-S
;	iplus		IPLUS
;	iminus		IMINUS
;	p		Feautrier mean intensity
;	s		Source function
;	tauq		Monochromatic optical depth
;	dtauq		dtauq(k)=tauq(k)-tauq(k-1)
;	xcont		continuum opacity relative to standard opacity
;	sc		absorption part of source function
;	scat		scattering part of source function
;	x		total opacity relative to standard opacity
;	jny		mean intensity
;	sbck		background source function, SBCK=SC+SCAT*JNY
;	rny		xcont/x
;
; COMMON BLOCKS:
;	common_multi
;
; RESTRICTIONS:
;	The jny file has to be converted to single precision - this is
;	NOT taken care of by the keyword /DP - or you have to use /JNYDP
;
; MODIFICATION HISTORY:
; 	Written by:	Mats Carlsson
;       95-11-30        JNYDP keyword added
;       04-09-16        added swap endian
;-
pro nyrd,kr,ny,file_idlny,file_jny,mdep=mdep,dp=dp,jnydp=jnydp

@common_multi
common cnyrd,openfile,lu2,ljny
common cnyrd2,nrec1,nrec2
common cnyrd3,count


if (n_params(0) lt 2) or (n_params(0) eq 3) then begin
  print,'nyrd,kr,ny [,file_idlny,file_jny,mdep=mdep,/dp,/jnydp]'
  return
endif

if(n_elements(dp) eq 0) then dp=0
if(n_elements(jnydp) ne 0) then dp=1
;
; new files should be opened if first call or if file parameter given
;
if (n_elements(openfile) eq 0) then openfile=0  ; if first call
if (n_params(0) gt 2) then openfile=0           ; if file parameter given:

if (openfile eq 0) then begin
  if (n_elements(lu2) ne 0) then if(lu2 ne 0) then free_lun,lu2 ; close open fil
  if (n_elements(ljny) ne 0) then if(ljny ne 0) then free_lun,ljny
  get_lun,lu2                                   ; get unit numbers
  get_lun,ljny
  if (n_elements(file_idlny) eq 0) then begin
    if(def_ext eq 'none') then file_idlny='IDLNY' else $
     file_idlny='idlny.'+def_ext
  endif
  if (n_elements(file_jny) eq 0) then begin
    if(def_ext eq 'none') then file_jny='JNY' else $
     file_jny='jny.'+def_ext
  endif
  openr,lu2,file_idlny,swap_endian=swap_endian   ; open files
  dum=findfile(file_jny,count=count)
  if(count eq 1) then begin
    openr,ljny,file_jny,swap_endian=swap_endian
  endif else begin
    print,'file: '+file_jny+' does not exist, setting jny=0'
  endelse
  if(n_elements(mdep) eq 0) then mdep=ndep
  nrec2=(ndep*11+1)/2
  nrec2=nrec2*2
  nrec1=ndep*11
  if(nrec2 eq nrec1) or (dp eq 0) then begin
    data_idl2=assoc(lu2,fltarr(ndep,11))
  endif else begin
    data_idl2=assoc(lu2,fltarr(nrec2))
  endelse
  if(count eq 1) then begin
    if(n_elements(jnydp) eq 0) then data_jny=assoc(ljny,fltarr(mdep)) $
     else data_jny=assoc(ljny,dblarr(mdep))
  endif
  openfile=1
endif
;
if (kr gt 0) then  irec=total(nq(0:kr-1))+ny  else irec=ny

data=data_idl2(irec)
if(nrec2 ne nrec1) then begin
  data=reform(data(0:nrec1-1),ndep,11,/overwrite)
endif
pms=data(0:ndep-1,0)
iplus=data(0:ndep-1,1)
iminus=data(0:ndep-1,2)
p=data(0:ndep-1,3)
s=data(0:ndep-1,4)
tauq=data(0:ndep-1,5)
dtauq=data(0:ndep-1,6)
xcont=data(0:ndep-1,7)
sc=data(0:ndep-1,8)
scat=data(0:ndep-1,9)
x=data(0:ndep-1,10)

if(count eq 1) then begin
  jny=data_jny(irec) 
  jny=jny(0:ndep-1)
endif else jny=fltarr(ndep)

if (kr gt nline-1) then begin
  kt=ktrans(kr)-1
  i=irad(kr)-1
  j=jrad(kr)-1
  hn3c2=2.*hh*frq(ny+1,kt)/cc*frq(ny+1,kt)/cc*frq(ny+1,kt)
  gij=transpose(nstar(i,*)/nstar(j,*))*exp(-hh*frq(ny+1,kt)/bk/temp)
  z(0)=n(i,*)-gij*n(j,*)
  sl(0,kr)=hn3c2*transpose(n(j,*))*gij/z
  bp(0,kr)=hn3c2/(exp(hh*frq(ny+1,kt)/bk/temp)+1.)
endif

sbck=sc+scat*jny
rny=xcont/x

end
