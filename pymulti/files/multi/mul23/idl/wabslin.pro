pro wabslin,lines ,file=file
;+
;   wabslin,lines ,file=file
;
;            writes abslin file 
;
;-
if(n_params() lt 1) then begin
  message,'wabslin,lines ,file=file',/info
  return
endif 

if(n_elements(file) eq 0) then file='ABSLIN'
files=findfile(file,count=count)
if(count eq 1) then begin
  text=''
  read,file+' exists, overwrite (y/n)? ',text
  if(strtrim(strlowcase(text),2) ne 'y') then return
endif
if(n_elements(lines) eq 0) then begin
  message,'lines does not exist',/info
  return
endif

nll=n_elements(lines)
openw,luw,file,/get_lun
printf,luw,nll,format='(i8)'
for i=0,nll-1 do begin
  printf,luw,lines[i].srcell
  printf,luw,lines[i].abnill
  printf,luw,lines[i].ionill,lines[i].chill,lines[i].gill,format='(I4,F18.3,f5.1)'
  printf,luw,lines[i].abnjll
  printf,luw,lines[i].ionjll,lines[i].chjll,lines[i].gjll,format='(I4,F18.3,f5.1)'
  printf,luw,lines[i].alamll,lines[i].bluell,lines[i].redll,format='(3f10.3)'
  printf,luw,lines[i].fll,lines[i].gall,lines[i].gwll,lines[i].gqll,format='(e13.4,e10.3,f8.2,e11.3)'
endfor
free_lun,luw

end
