pro atom_renum,filein,index, fileout,remove=remove
;+
;   atom_renum,filein,index, fileout,remove=remove
;
;            reads multi-format atom file and renumbers levels
;            index(i) gives the new level number for level i
;            (both in IDL convention with 0 as lowest index)
;            negative value indicates that the level and all
;            transitions involving the level should be removed
;            output on fileout if filename given
;            with the remove keyword you can give the numbers of the
;            levels that are to be removed and index then becomes an output
;
;-
@common_multi

if (n_params(0) lt 1) then begin
  print,'atom_renum,filein ,index,fileout,remove=remove'
  return
endif

; read input file

ratom,filein

if(keyword_set(remove)) then begin
  index=intarr(nk)+1
  index[0]=0
  index[remove]=0
  index=round(total(index,/cumulative))
  index[remove]=-1
endif

; printout

if(n_elements(fileout) eq 0) then printout=0 else printout=1

if(printout eq 1) then begin
  openw,luw,fileout,/get_lun
endif

; copy comment lines at beginning of file
; position atomic file at GENCOL data

nt=0
openr,lu1,filein,/get_lun
text=''
readf,lu1,text
if(printout eq 1) then begin
  while (strmid(text,0,1) eq '*') do begin
    printf,luw,text
    readf,lu1,text
  endwhile
endif
while (strmid(strtrim(text,2),0,6) ne 'GENCOL') do readf,lu1,text

; renumber

nk2=max(index)+1
g2=g
label2=label
ion2=ion
ev2=ev
g=fltarr(nk2)
label=strarr(nk2)
ion=intarr(nk2)
ev=dblarr(nk2)
for i=0,nk-1 do begin
  ii=index(i)
  if(ii ge 0) then begin
    g(ii)=g2(i)
    label(ii)=label2(i)
    ion(ii)=ion2(i)
    ev(ii)=ev2(i)
  endif
endfor
if(nline gt 0) then begin
  f2=f
  nq2=nq
  qmax2=qmax
  q02=q0
  iwide2=iwide
  ga2=ga
  gw2=gw
  gq2=gq
  q2=q
  irad2=irad
  jrad2=jrad
  krad2=intarr(nk2,nk2)            ; kr corresponding to new i,j
  for KR=0,NLINE-1 do begin
    jj=index(jrad(kr)-1)
    ii=index(irad(kr)-1)
    j=ii > jj
    i=ii < jj
    if(i ge 0) and (j ge 0) then begin
      krad2(j,i)=kr+1
      krad2(i,j)=kr+1
    endif
  endfor
  kr2=-1
  for i=0,nk2-2 do begin
    for j=i+1,nk2-1 do begin
      kr=krad2(i,j)-1
      if(kr ge 0) then begin
        kr2=kr2+1
        ii=index(irad2(kr)-1)+1
        jj=index(jrad2(kr)-1)+1
        irad(kr2)=ii < jj
        jrad(kr2)=ii > jj
        f(kr2)=f2(kr)
        nq(kr2)=nq2(kr)
        qmax(kr2)=qmax2(kr)
        q0(kr2)=q02(kr)
        iwide(kr2)=iwide2(kr)
        ga(kr2)=ga2(kr)
        gw(kr2)=gw2(kr)
        gq(kr2)=gq2(kr)
        IF(QMAX2(KR) LT 0.0) OR (Q02(KR) LT 0.0) THEN begin
          q(0:nq(kr2)-1,kr2)=q2(0:nq2(kr)-1,kr)
        endif
      endif
    endfor
  endfor
endif
nline2=kr2+1              ; new number of lines

;  BOUND-FREE TRANSITIONS IN DETAIL

if (nrad gt nline) then begin
  alfac2=alfac
  i=index(irad2(nline:nrad-1)-1)+1
  iw=sort(i)+nline
  for KRR=NLINE,NRAD-1 do begin
    kr=iw(krr-nline)
    if(kr ge 0) then begin
      i=index(irad2(kr)-1)+1
      j=index(jrad2(kr)-1)+1
      if(i gt 0) and (j gt 0) then begin
        kr2=kr2+1
        irad(kr2)=i < j
        jrad(kr2)=i > j
        f(kr2)=f2(kr)
        nq(kr2)=nq2(kr)
        qmax(kr2)=qmax2(kr)
        IF(QMAX2(KR) LT 0.0) THEN begin
          q(0:nq(kr2)-1,kr2)=q2(0:nq2(kr)-1,kr)
          alfac(0:nq(kr2)-1,kr2-nline2)=alfac2(0:nq2(kr)-1,kr-nline)
        ENDIF
      endif
    endif
  endfor
endif
nrad2=kr2+1
;
;  FIXED TRANSITIONS PARAMETERS
;
kf=-1
if (nrfix gt 0) then begin
  jfx2=jfx
  ifx2=ifx
  for KFF=0,NRFIX-1 do begin
    j=index(jfx2(kff)-1)+1
    i=index(ifx2(kff)-1)+1
    if(i gt 0) and (j gt 0) then begin
      kf=kf+1
      jfx(kf)=i > j
      ifx(kf)=i < j
      IF(ITRAD(KF) EQ 4) THEN begin
        print,'itrad=4 NOT IMPLEMENTED'
        free_lun,lu1
        if(printout eq 1) then free_lun,luw
        return
      endif
    endif
  endfor
endif
nrfix2=kf+1

if(printout eq 1) then begin

; read GENCOL data and write renumbered data to gencol.tmp
  ckey=['OHM','CE','CP','CH','CALP','RECO','CH0','CH+','CI']
  nkey=n_elements(ckey)
  openw,lu2,'gencol.tmp',/get_lun
  printf,lu2,' GENCOL'
  comm=''
  gencol_loop:
    readf,lu1,text
    key=strtrim(strmid(text,0,20),2)
    if(strmid(key,0,1) eq '*') then begin
      comm=text
      goto,gencol_loop
    endif
; special keys
    if(strmid(key,0,4) eq 'TEMP') then begin
      if(comm ne '') then printf,lu2,comm
      printf,lu2,text
      readf,lu1,text
      dum=strsplit(text,/extract)
      nt=fix(dum[0])
      nt1=n_elements(dum)-1
      temp=fltarr(nt)
      if(nt1 gt 0) then temp[0:nt1-1]=dum[1:nt1]
      if(nt1 lt nt) then begin
        dum=fltarr(nt-nt1)
        readf,lu1,dum
        temp[nt1:nt-1]=dum
      endif
      printf,lu2,nt,temp[0:5<(nt-1)],format='(i4,6e14.5)'
      if(nt gt 6) then printf,lu2,temp[6:nt-1],format='(4x,6e14.5)'
      comm=''
      goto,gencol_loop
    endif else if (strmid(key,0,4) eq 'SEMI') then begin
      dum=fltarr(1)
      readf,lu1,ihi,ilo,dum
      if(index(ihi-1) ge 0) and (index(ilo-1) ge 0) then begin
        if(comm ne '') then printf,lu2,comm
        printf,lu2,text
        printf,lu2,index(ihi-1)+1,index(ilo-1)+1,dum,format='(2i4,2x,7e9.2)'
      endif
      comm=''
      goto,gencol_loop
    endif else if (strmid(key,0,7) eq 'BURGESS') then begin
      dum=fltarr(1)
      readf,lu1,ihi,ilo,dum
      if(index(ihi-1) ge 0) and (index(ilo-1) ge 0) then begin
        if(comm ne '') then printf,lu2,comm
        printf,lu2,text
        printf,lu2,index(ihi-1)+1,index(ilo-1)+1,dum,format='(2i4,2x,7e9.2)'
      endif
      comm=''
      goto,gencol_loop
    endif else if (strmid(key,0,7) eq 'SHULL82') then begin
      dum=fltarr(8)
      readf,lu1,ihi,ilo,dum
      if(index(ihi-1) ge 0) and (index(ilo-1) ge 0) then begin
        if(comm ne '') then printf,lu2,comm
        printf,lu2,text
        printf,lu2,index(ihi-1)+1,index(ilo-1)+1,dum,format='(2i4,2x,8e9.2)'
      endif
      comm=''
      goto,gencol_loop
    endif else if (strmid(key,0,8) eq 'AR85-CDI') then begin
      ndum=0L
      dum=fltarr(5)
      readf,lu1,ihi,ilo,ndum
      if(index(ihi-1) ge 0) and (index(ilo-1) ge 0) then begin
        if(comm ne '') then printf,lu2,comm
        printf,lu2,text
        printf,lu2,index(ihi-1)+1,index(ilo-1)+1,format='(2i4)'
        printf,lu2,ndum,format='(i4)'
        for idum=0,ndum-1 do begin
          readf,lu1,dum
          printf,lu2,dum,format='(5f9.2)'
        endfor
      endif else begin
        for idum=0,ndum-1 do begin
          readf,lu1,dum
        endfor
      endelse
      comm=''
      goto,gencol_loop
    endif else if (strmid(key,0,8) eq 'AR85-CEA') then begin
      dum=fltarr(1)
      readf,lu1,ihi,ilo,dum
      if(index(ihi-1) ge 0) and (index(ilo-1) ge 0) then begin
        if(comm ne '') then printf,lu2,comm
        printf,lu2,text
        printf,lu2,index(ihi-1)+1,index(ilo-1)+1,format='(2i4)'
        printf,lu2,dum,format='(7e9.2)'
      endif
      comm=''
      goto,gencol_loop
    endif else if (strmid(key,0,7) eq 'AR85-CH') then begin
      dum=fltarr(6)
      readf,lu1,ihi,ilo,dum
      if(index(ihi-1) ge 0) and (index(ilo-1) ge 0) then begin
        if(comm ne '') then printf,lu2,comm
        printf,lu2,text
        printf,lu2,index(ihi-1)+1,index(ilo-1)+1,format='(2i4)'
        printf,lu2,dum,format='(3e9.2,3F9.2)'
      endif
      comm=''
      goto,gencol_loop
    endif else if (strmid(key,0,4) eq 'LTDR') then begin
      dum=fltarr(5)
      readf,lu1,ihi,ilo,dum
      if(index(ihi-1) ge 0) and (index(ilo-1) ge 0) then begin
        if(comm ne '') then printf,lu2,comm
        printf,lu2,text
        printf,lu2,index(ihi-1)+1,index(ilo-1)+1,dum,format='(2i4,2x,7e9.2)'
      endif
      comm=''
      goto,gencol_loop
    endif else if (strmid(key,0,6) eq 'CORONA') then begin
      readf,lu1,ihi,ilo
      if(index(ihi-1) ge 0) and (index(ilo-1) ge 0) then begin
        if(comm ne '') then printf,lu2,comm
        printf,lu2,text
        printf,lu2,index(ihi-1)+1,index(ilo-1)+1,format='(2i4)'
      endif
      comm=''
      goto,gencol_loop
    endif
; keys with nt data points
    for i=0,nkey-1 do begin
      key0=strmid(key,0,strlen(ckey(i)))
      if(key0 eq ckey(i)) then begin
        dum=fltarr(nt)
        readf,lu1,ihi,ilo,dum
        if(index(ihi-1) ge 0) and (index(ilo-1) ge 0) then begin
          if(comm ne '') then printf,lu2,comm
          printf,lu2,text
          printf,lu2,index(ihi-1)+1,index(ilo-1)+1,dum[0:7<(nt-1)],format='(2i4,2x,8e10.2)'
          if(nt gt 8) then printf,lu2,dum[8:nt-1],format='(10x,8e10.2)'
        endif
        comm=''
        goto,gencol_loop
      endif
    endfor
    if(strmid(key,0,3) ne 'END') then begin
      print,'error in gencol data, text=',text
      free_lun,lu1
      free_lun,lu2
      free_lun,luw
      return
    endif else printf,lu2,text
  free_lun,lu2
endif
      
free_lun,lu1

nk=nk2             ; new number of levels
nline=nline2
nrad=nrad2
nrfix=nrfix2

; write atomic data

if(printout eq 1) then begin
  watom,fileout,luw=luw

; copy over gencol.tmp to fileout

  openr,lu1,'gencol.tmp',/get_lun
  while (not eof(lu1)) do begin
    readf,lu1,text
    printf,luw,text
  endwhile
  free_lun,lu1
  free_lun,luw
endif

RETURN
END



