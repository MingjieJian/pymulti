pro c_gencolw,t,coll,clab ,file=file
;+
;   c_gencolw,t,coll,clab ,file=file
;
;            writes gencol data to file gencol.dat
;
;-
COMMON CLABEL,LABEL,ATOMID,CROUT
                                                                                
COMMON CLIMIT,NDEP,NK,NLINE,NWIDE,NRAD,NRFIX,NMU,MQ,NQ

if (n_params(0) ne 3) then begin
  print,'c_gencolw,t,coll,clab ,file=file'
  return
endif

if(n_elements(file) eq 0) then file='gencol.dat'

nt=n_elements(t)
get_lun,lu
openw,lu,file
printf,lu,' TEMP'
printf,lu,nt,t,format='(i3,(20f10.0))'

for i=0,nk-3 do begin
  for j=i+1,nk-2 do begin
    if(max(coll(*,j,i)) gt 0.0) then begin
      printf,lu,label(j),label(i),format="('* ',a,' -> ',a)"
      printf,lu,clab(j,i)
      printf,lu,j+1,i+1,coll(*,j,i),format='(2i4,(20e10.2))'
    endif
  endfor
endfor

j=nk-1
for i=0,nk-2 do begin
  if(max(coll(*,j,i)) gt 0.0) then begin
    printf,lu,label(j),label(i),format="('* ',a,' -> ',a)"
    printf,lu,clab(j,i)
    printf,lu,j+1,i+1,coll(*,j,i),format='(2i4,(20e10.2))'
  endif
endfor

free_lun,lu

return
end
