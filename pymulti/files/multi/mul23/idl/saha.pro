pro saha,t,nne,xi,u0,u1,n1_ntot,n0_ntot
;+
;   saha,t,nne,xi,u0,u1,n1_ntot,n0_ntot
;
;            calculates saha ionization
;
;            parameters:
;
;            t      temperature K
;            nne    electron density
;            xi     ionization potential (eV)
;            u0     partition function of neutral element
;            u1     partition function of ionized element
;            n1_ntot ionization fraction
;            n0_ntot neutral fraction
;
;            all atoms are supposed to be either neutral or once ionized
;
;-
if(n_params(0) lt 7) then begin
  print,'saha,t,nne,xi,u0,u1,n1_ntot,n0_ntot'
  return
endif

bk= 1.38066E-16

theta=5039.77/t
phit=0.6665*u1/u0*t^2.5*10.^(-theta*xi)       ; function phi of Grey
n1_n0=phit/(nne*bk*t)                         ; N1/N0
n1_ntot=n1_n0/(1.+n1_n0)
n0_ntot=1/(1.+n1_n0)

end
