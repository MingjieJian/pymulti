pro hpop,temp,nne,rho,nh

if(n_params() lt 4) then begin
  print,'hpop,temp,nne,rho,nh
  return
endif

em=9.109534E-28
bk=1.380662E-16
hh=6.626176E-27
ee=1.602189E-12
nh=fltarr(6,n_elements(temp))
for k=0,n_elements(temp)-1 do begin
    phit=(2.d0*!pi*em*bk*temp[k]/hh/hh)^1.5*$
      exp(-13.6d0*ee/bk/temp[k])
    ratio=phit/nne[k]
    totnh=rho[k]/2.38049d-24
    nh[5,k]=ratio/(1.d0+ratio)*totnh
    ratio=nne[k]/phit
    nh[0,k]=ratio/(1.d0+ratio)*totnh
    for i=2,5 do begin
        nh[i-1,k]=i*i*exp(-13.6d0*(1.d0-1.d0/i/i)*$
                          ee/bk/temp[k])*nh[0,k]
    endfor
endfor

end
