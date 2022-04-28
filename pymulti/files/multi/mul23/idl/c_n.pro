pro c_n,label,n
;+
;   c_n,label,n
;
;            finds n from label, leading number in group before parity
;
;-
if (n_params(0) ne 2) then begin
  print,'c_n,label,n'
  return
endif

digits='0123456789'
nk=n_elements(label)
n=intarr(nk-1)
for i=0,nk-2 do begin
  labtrim=strtrim(label(i),2)
  group_pos=intarr(10)           ; position for group of characters
  jgroup=-1                       ; index for group of characters
  cjm=' '                             ; previous character
  for j=0,strlen(labtrim)-1 do begin  ; go through characters one by one
    cj=strmid(labtrim,j,1)            ; extract one character
    if(cj ne ' ') and (cjm eq ' ') then begin   ; new word begins here
      jgroup=jgroup+1
      group_pos(jgroup)=j
    endif
    cjm=cj
  endfor
; find group where parity is given
  find_parity:
  last_group=strmid(labtrim,group_pos(jgroup),10)
  if(strpos(last_group,'E') lt 0) and (strpos(last_group,'O') lt 0) then begin
    jgroup=jgroup-1
    goto,find_parity
  endif
; check number of digits in n
  ipos=group_pos(jgroup-1)            ; next to last group of characters
  if(strpos(digits,strmid(labtrim,ipos+1,1)) ge 0) then il=2 else il=1 
  n(i)=long(strmid(labtrim,ipos,il))
endfor

return
end
