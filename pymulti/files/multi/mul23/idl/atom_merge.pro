pro atom_merge,filein,merge_index,fileout ,new_label=new_label,$
 index_arr=index_arr
;+
;   atom_merge,filein,merge_index,fileout ,new_label=new_label, 
;              index_arr=index_arr
;
;            merges two levels (IDL index given in merge_index)
;            index_arr[i,0] gives number of levels in original atom
;            that have been merged into new level i, index_arr[i,1:n]
;            gives indices in original atom. If not given on input it
;            will be initialized to non-merged atom
;
;-
@common_multi

if (n_params(0) lt 3) then begin
  print,'atom_merge,filein,merge_index,fileout ,new_label=new_label,$'
  print,' index_arr=index_arr'
  return
endif

if(n_elements(merge_index) ne 2) then begin
  print,'index array has to have two elements'
  return
endif

ii=merge_index(0) < merge_index(1)
jj=merge_index(0) > merge_index(1)

if(n_elements(fileout) eq 0) then printout=0 else printout=1

if(printout eq 1) then begin
  openw,luw,fileout,/get_lun
endif

; copy comment lines at beginning of file if fileout given
; find number of temperature points in GENCOL data

if(n_elements(filein) ne 0) or (printout eq 1) then begin
  ratom,filein

  if(ii lt 0) or (jj gt nk-1) then begin
    print,'indices have to be in range [0,nk-1]'
    return
  endif

; weighting factors

  gsum=g(ii)+g(jj)
  wii=g(ii)/gsum
  wjj=g(jj)/gsum

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
  while (strmid(strtrim(text,2),0,4) ne 'TEMP') do readf,lu1,text
  readf,lu1,text
  nt=fix(text)
  text=strcompress(strtrim(text,2))
  t=fltarr(nt)
  ipos=strpos(text,' ')+1
  ilen=strlen(text)-ipos+1
  text=strmid(text,ipos,ilen)
  for i=0,nt-1 do begin
    t(i)=float(text)
    ipos=strpos(text,' ')+1
    ilen=strlen(text)-ipos+1
    text=strmid(text,ipos,ilen)
  endfor
  free_lun,lu1
; read GENCOL data
  c_colrd,nk,nt,filein,coll,clabi ,clab
endif

; merge

if(n_elements(index_arr) eq 0) then begin
  index_arr=intarr(nk,2)
  index_arr[*,0]=1
  index_arr[*,1]=indgen(nk)
endif

index_arr_old=index_arr
n_merge=index_arr_old[ii,0]+index_arr_old[jj,0]
old_dim=n_elements(index_arr_old[0,*])
new_dim=old_dim > n_merge+1
index_arr=intarr(nk,new_dim)
index_arr[0:nk-1,0:old_dim-1]=index_arr_old[0:nk-1,0:old_dim-1]
index_arr[ii,0]=n_merge
i_merge=[reform(index_arr_old[ii,1:index_arr_old[ii,0]]),$
         reform(index_arr_old[jj,1:index_arr_old[jj,0]])]
iw=sort(i_merge)
index_arr[ii,1:n_merge]=i_merge[iw]

ev(ii)=wii*ev(ii)+wjj*ev(jj)
g(ii)=g(ii)+g(jj)
;l=strlen(label(ii)) < strlen(label(jj))
;k3=0
;label_loop:
;  if(strmid(label(ii),k3,1) eq strmid(label(jj),k3,1)) then begin
;    k3=k3+1
;    goto,label_loop
;  endif
;label(ii)=strmid(label(ii),0,k3)          ; set label to common part
l=strlen(strtrim(label[ii]))
new_label=label[ii]
if(strupcase(strmid(label[ii],l-1,1)) ne 'M') then strput,new_label,'m',l
label[ii]=new_label

; merge bb transitions

if(nline gt 0) then begin
  f2=f
  irad2=irad
  jrad2=jrad
  for KR=0,NLINE-1 do begin
    i=irad2(kr)-1
    j=jrad2(kr)-1
    if(j eq ii) then begin
      if(krad(i,jj) ne 0) then begin
        kr2=krad(i,jj)-1
        f(kr)=f2(kr)+f2(kr2)
      endif
    endif else if(j eq jj) then begin
      if(krad(i,ii) eq 0) then begin
        jrad(kr)=ii+1
      endif
    endif else if(i eq ii) then begin
      if(krad(jj,j) ne 0) then begin
        kr2=krad(jj,j)-1
        f(kr)=wii*f2(kr)+wjj*f2(kr2)
      endif else begin
        f(kr)=wii*f2(kr)
      endelse
    endif else if(i eq jj) then begin
      if(krad(ii,j) eq 0) then begin
        f(kr)=wjj*f2(kr)
        irad(kr)=ii+1
      endif
    endif
  endfor
endif

; merge photoionization continua
; NB does not necessarily conserve integral, simple interpolation of
; one set of wavelength points to the other

krii=krad[ii,nk-1]-1      ; photoionization from lower level of merge
krjj=krad[jj,nk-1]-1      ; photoionization from upper level of merge

if(krii ge 0) then begin
  frqii=frq[1:nq[krii],krii-nline]
  if(krjj ge 0) then begin
    frqjj=frq[1:nq[krjj],krjj-nline]
    intep,frqjj,alfac[0:nq[krjj]-1,krjj-nline],frqii,alfacjj
    alfac[0:nq[krii]-1,krii-nline]=(g[ii]*alfac[0:nq[krii]-1,krii-nline]+$
     g[jj]*alfacjj)/(g[ii]+g[jj])
  endif else begin
    print,'no photoionization data for upper level, assuming same as lower'
  endelse
endif else if(krjj ge 0) then begin
  print,'no photoionization data for lower level, assuming same as upper'
  irad[krjj]=ii
  krad[ii,nk-1]=krjj+1
  krad[nk-1,ii]=krjj+1
endif

; merge collisions

coll2=coll
clab2=clab
coll(*,ii,*)=wii*coll2(*,ii,*)+wjj*coll2(*,jj,*)
coll(*,*,ii)=coll2(*,*,ii)+coll2(*,*,jj)

; check labels

for i=0,nk-1 do begin
  if(clab[ii,i] eq '' and clab[jj,i] ne '') then clab[ii,i]=clab[jj,i]
  if(clab[i,ii] eq '' and clab[i,jj] ne '') then clab[i,ii]=clab[i,jj]
endfor

ev2=ev
ev2(jj)=1.e30
iw=sort(ev2)
index=intarr(nk)
index(iw)=indgen(nk)
index(jj)=-1              ; delete jj level

index_arr_unsort=index_arr  ; unsorted copy of index_arr
for i=0,nk-1 do begin
  if(index[i] ge 0) then index_arr[index[i],*]=index_arr_unsort[i,*]
endfor
index_arr=index_arr[0:nk-2,*]

; write data to temporary file atom.tmp

openw,luw,'atom.tmp',/get_lun
watom,'atom.tmp',luw=luw

c_gencolw,t,coll,clab
printf,luw,' GENCOL'
openr,lu1,'gencol.dat',/get_lun
text=''
while (not eof(lu1)) do begin
  readf,lu1,text
  printf,luw,text
endwhile
printf,luw,' END'
free_lun,lu1
free_lun,luw

; renumber levels

atom_renum,'atom.tmp',index,fileout

return
end
