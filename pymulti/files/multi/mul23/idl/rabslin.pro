pro rabslin,lines ,file=file
;+
;   rabslin,lines ,file=file
;
;            reads abslin file and put data into structure array lines
;
;-
if(n_params() lt 1) then begin
  message,'rabslin,lines ,file=file',/info
  return
endif 

if(n_elements(file) eq 0) then file='ABSLIN'
files=findfile(file,count=count)
if(count ne 1) then begin
  message,file+' not found',/info
  return
endif

cstrip,file,'~/dum.tmp'     ; strip abslin file from comment lines
openr,lur,'~/dum.tmp',/get_lun
nll=0L
readf,lur,nll

lines_structure={lines_str,srcell:'',abnill:'',ionill:0,chill:0.,gill:0.,$
                                     abnjll:'',ionjll:0,chjll:0.,gjll:0.,$
                 alamll:0.,bluell:0.,redll:0.,fll:0.,gall:0.,gwll:0.,gqll:0.}
lines=replicate(lines_structure,nll)
text=''
ii=0L
chi=0.
g=0.
alam=0.
blue=0.
red=0.
f=0.
ga=0.
gw=0.
gq=0.
for i=0,nll-1 do begin
  readf,lur,text
  lines[i].srcell=text
  readf,lur,text
  lines[i].abnill=text
  readf,lur,ii,chi,g
  lines[i].ionill=ii
  lines[i].chill=chi
  lines[i].gill=g
  readf,lur,text
  lines[i].abnjll=text
  readf,lur,ii,chi,g
  lines[i].ionjll=ii
  lines[i].chjll=chi
  lines[i].gjll=g
  readf,lur,alam,blue,red
  lines[i].alamll=alam
  lines[i].bluell=blue
  lines[i].redll=red
  readf,lur,f,ga,gw,gq
  lines[i].fll=f
  lines[i].gall=ga
  lines[i].gwll=gw
  lines[i].gqll=gq
endfor
free_lun,lur

end
