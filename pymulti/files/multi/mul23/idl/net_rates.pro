pro net_rates
;+
;   net_rates
;
;            calculates net rates per particle
;            pr(k,i,j)=net radiative rate, positive=i->j, 
;            pc(k,j,i)=net collisional rate
;
;-
@common_multi 

pr=fltarr(ndep,nk,nk)
pc=fltarr(ndep,nk,nk)
 
for i=0,nk-2 do begin
  for j=i+1,nk-1 do begin
    pc(*,i,j)=pc(*,i,j)+(n(i,*)*c(i,j,*)- $
                         n(j,*)*c(j,i,*))/totn(*)
    pc(*,j,i)=-pc(*,i,j)
  endfor
endfor

for kr=0,nrad-1 do begin
  i=irad(kr)-1
  j=jrad(kr)-1
  pr(*,i,j)=pr(*,i,j)+(n(i,*)*rij(*,kr)-$
                       n(j,*)*rji(*,kr))/totn(*)
  pr(*,j,i)=-pr(*,i,j)
endfor

return 
end
