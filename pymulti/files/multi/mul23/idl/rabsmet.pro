pro rabsmet,wnum ,file=file
;+
;   rabsmet,wnum ,file=file
;
;            reads osmet file ABSMET
;
;-
if(n_params() lt 1) then begin
  print,'rabsmet,wnum ,file=file'
  return
endif

if(n_elements(file) eq 0) then file='ABSMET'

openr,lur,file,/get_lun,/f77_unformatted

textatom=''
readu,lur,textatom
readu,lur,textatom
ispec=lonarr(2*92)
readu,lur,ispec
xiturb=0.
readu,lur,xiturb
maxp6=0L
maxtem=0L
maxnu=0L
mxi=1L
readu,lur,maxp6
p6pt=fltarr(maxp6)
readu,lur,p6pt
readu,lur,maxtem
tempt=fltarr(maxtem)
readu,lur,tempt
text=''
for i=1,4 do readu,lur,text
mnu=1000000L
wnum=dblarr(mnu)
wnu=0.d0
i=0L
while(not eof(lur)) do begin
  readu,lur,wnu
  wnum[i]=wnu
  i=i+1L
endwhile
nnu=i

wnum=wnum[0:nnu-1]
free_lun,lur
end
