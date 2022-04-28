pro cstrip,filein,fileout
;+
;   cstrip,filein,fileout
;
;            strip file from comment lines
;
;-
if(n_params() lt 2) then begin
  message,'cstrip,filein,fileout',/info
  return
endif
files=findfile(filein,count=count)
if(count ne 1) then begin
  message,filein+' not found',/info
  return
endif

openr,lur,filein,/get_lun
openw,luw,fileout,/get_lun
text=''
while (not eof(lur)) do begin
  readf,lur,text
  if(strpos(text,'*') ne 0) then printf,luw,text
endwhile
free_lun,lur
free_lun,luw

end
