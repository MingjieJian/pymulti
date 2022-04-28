pro c_l,label,l
;+
;   c_l,label,l
;
;            finds l from label, 1st character before multiplicity
;
;-
if (n_params(0) ne 2) then begin
  print,'c_l,label,l'
  return
endif

nk=n_elements(label)
l=intarr(nk-1)
for i=0,nk-2 do begin
  labtrim=strtrim(label(i),2)
  ipos=strpos(labtrim,' ')
  ipos2=strpos(labtrim,'E',ipos)                ; look for parity designation
  if(ipos2 lt 0) then ipos2=strpos(labtrim,'O',ipos)
  if(ipos2 lt 0) then begin
    print,' c_l: parity not given in label:'
    print,labtrim
    return
  endif
  ipos3=ipos2-1
  loop:
    ipos3=ipos3-1
    if(ipos3 lt 0) then begin
      print,'c_l: angular quantum number not found, label:'
      print,labtrim
      return
    endif
    lc=strmid(labtrim,ipos3,1)   
    l(i)=strpos('SPDFGHIKLMNOQRTUVXYZ',lc)
  if(l[i] lt 0) then goto,loop
  if(l(i) lt 0) then begin
    print,'c_l: angular quantum number not found, label:'
    print,labtrim
    return
  endif
endfor

return
end
