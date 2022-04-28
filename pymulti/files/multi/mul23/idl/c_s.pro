pro c_s,label,s
;+
;   c_s,label,s
;
;            finds s from label, 2nd character from parity
;
;-
if (n_params(0) ne 2) then begin
  print,'c_s,label,s'
  return
endif

nk=n_elements(label)
s=strarr(nk-1)
for i=0,nk-2 do begin
  labtrim=strtrim(label(i),2)
  ipos=strpos(labtrim,' ')
  ipos2=strpos(labtrim,'E',ipos)                ; look for parity designation
  if(ipos2 lt 0) then ipos2=strpos(labtrim,'O',ipos)
  if(ipos2 lt 0) then begin
    print,' c_s: parity not given in label:'
    print,labtrim
    return
  endif
  s(i)=strmid(labtrim,ipos2-2,1)
endfor
s=fix(s)

return
end
