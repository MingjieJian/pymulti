pro c_colrd,nk,nt,file,coll,clabi ,clab
;+
;   c_colrd,nk,nt,file,coll,clabi [,clab]
;
;            reads collisonal data from file
;
;            clabi is collisional option (integer*1)
;            clab  is collisional label
;
;-
if (n_params(0) lt 5) then begin
  print,'c_colrd,nk,nt,file,coll,clabi [,clab]'
  return
endif
coll=fltarr(nt,nk,nk)
clabi=intarr(nk,nk)
clab=strarr(nk,nk)
dum=fltarr(nt)
text=''
text2=''
get_lun,lu
get_lun,lu2
openr,lu,file
openw,lu2,'dums.dat'
while not eof(lu) do begin                        ; strip away comment lines
  readf,lu,text
  if(strmid(text,0,1) ne '*') then printf,lu2,text
endwhile
close,lu
close,lu2
openr,lu,'dums.dat'

loop:
  readf,lu,text
if(strtrim(text,2) ne 'GENCOL') then goto,loop
readf,lu,text
readf,lu,text
nn=long(nk)*(nk-1)/2
ii=0l
text1=' '
text2=' '
on_ioerror,end_read
for ix=1L,nn do begin
  readf,lu,text
  if(strpos(text,'TEMP') ge 0) then begin
    readf,lu,text
    readf,lu,text
  endif
  text1='                    '
  strput,text1,strmid(text,0,20),0
  ii=fix(strmid(text,20,1))
  text2=strmid(text,21,60)
  if(strtrim(text1,2) eq 'END') then goto,end_read
  readf,lu,j,i,dum
  j=j-1
  i=i-1
  clabi(j,i)=ii
  clabi(i,j)=ii
  clab(i,j)=text1+strtrim(string(ii),2)+text2
  clab(j,i)=clab(i,j)
  coll(0,j,i)=dum  
  coll(0,i,j)=dum  
endfor
ix=nn+1L

end_read:
print,' number of rates read: ',ix-1L
free_lun,lu
free_lun,lu2
return
end
