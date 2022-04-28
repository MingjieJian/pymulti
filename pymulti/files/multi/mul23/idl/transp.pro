pro transp,xmy,tau,x,s,p ,iplus,iminus,pms,tauq,itran=itran,iplus_0=iplus_0
;
;  solves the radiative transfer equation with given source function.
;  solves for p-s where p is feautriers p, i.e.
;
;    p = 0.5*(i(xmu)+i(-xmu))
;
;  in the presence of velocity fields:
;
;    p = 0.5*(i(ny,xmu) + i(-ny,-xmu))
;
;  itran determines the mode of formal solution
;  itran  = 0  feautrier solution
;           1  feautrier solution to cubic spline accuracy, ref:
;              kunasz, hummer, 1974, mnras 166,19
;              mihalas, 1974, apj suppl 28,343
;           2  feautrier solution hermite, ref:
;              auer, 1976, jqsrt 16,931
;
;  variables:
;
;  xmu   : angular cosine                                       (in)
;  tau   : standard optical depth scale                         (in)
;  x     : ratio of monochromatic to standard opacity           (in)
;  s     : monochromatic source function                        (in)
;
;  p     : mean bidirectional intensity (cf. above)            (out)
;  iplus : radiation intensity, outgoing rays                  (out)
;  iminus: radiation intensity, ingoing rays                   (out)
;  pms   : p-s                                                 (out)
;  tauq  : monochromatic optical depth                         (out)
;
if(n_params(0) lt 5) then begin
  print,'transp,xmy,tau,x,s,p ,iplus,iminus,pms,tauq'
  return
endif

if(n_elements(itran) eq 0) then itran=0
ndep=n_elements(tau)
sp1=fltarr(ndep)
sp2=fltarr(ndep)
sp3=fltarr(ndep)
a1=fltarr(ndep)
c1=fltarr(ndep)
p=fltarr(ndep)
iplus=fltarr(ndep)
iminus=fltarr(ndep)
pms=fltarr(ndep)
tauq=fltarr(ndep)
dtauq=fltarr(ndep)
;
; k=1: upper boundary
;
cmu=0.5/xmy
dtauq(1)=(x(0)+x(1))*(tau(1)-tau(0))*cmu
a1(0)=1./dtauq(1)
t=tau(0)*x(0)*2.0*cmu

tauq(0)=t
dtauq(0)=t
;
;  calculate dtauq
;
dtauq(1)=(x(1:ndep-1)+x(0:ndep-2))*(tau(1:ndep-1)-tau(0:ndep-2))*cmu
;
;  calculate tauq
;
for k=1,ndep-1 do begin
  tauq(k)=tauq(k-1)+dtauq(k)
endfor
;
;  calculate a1 and c1
;
a1(1)=2./(dtauq(1:ndep-2)+dtauq(2:ndep-1))/dtauq(1:ndep-2)
c1(1)=2./(dtauq(1:ndep-2)+dtauq(2:ndep-1))/dtauq(2:ndep-1)
;
;  call formal solver
;
t=tauq(0)
sp1(0)=0.
sp2(0)=1.+2.*a1(0)
sp3(0)=-2.*a1(0)*a1(0)
if (t lt 0.01) then begin
  ex1=t*(1.-t*(0.5-t*(0.1666667-t*0.041666667)))
endif else begin
  if(t lt 20.) then ex1=1.-exp(-t) else ex1=1.
endelse
ex=1.-ex1
fact=1.+2.*a1(0)*ex1
sp2(0)=sp2(0)/fact
sp3(0)=sp3(0)/fact
p(0)=s(0)
;
;  calculate tridiagonal coefficients
;
if(itran eq 0) then begin
  sp1(1)=-a1(1:ndep-2)
  sp2(1:ndep-2)=1.
  sp3(1)=-c1(1:ndep-2)
  p(1)=s(1:ndep-2)
endif else if(itran eq 1) then begin
  ad=.166666666*dtauq(1:ndep-2)*2./(dtauq(1:ndep-2)+dtauq(2:ndep-1))
  cd=.166666666*dtauq(2:ndep-1)*2./(dtauq(1:ndep-2)+dtauq(2:ndep-1))
  sp1(1)=-a1(1:ndep-2)+ad
  sp2(1:ndep-2)=1.
  sp3(1)=-c1(1:ndep-2)+cd
  p(1)=s(1:ndep-2)+ad*(s(0:ndep-3)-s(1:ndep-2))+cd*(s(2:ndep-1)-s(1:ndep-2))
endif else if(itran eq 2) then begin
  ad=.166666666-0.083333333*dtauq(2:ndep-1)^2/dtauq(1:ndep-2)* $
   2./(dtauq(1:ndep-2)+dtauq(2:ndep-1))
  cd=.166666666-0.083333333*dtauq(1:ndep-2)^2/dtauq(2:ndep-1)* $
   2./(dtauq(1:ndep-2)+dtauq(2:ndep-1))
  sp1(1)=-a1(1:ndep-2)+ad
  sp2(1:ndep-2)=1.
  sp3(1)=-c1(1:ndep-2)+cd
  p(1)=s(1:ndep-2)+ad*(s(0:ndep-3)-s(1:ndep-2))+cd*(s(2:ndep-1)-s(1:ndep-2))
endif else begin
  print,' tranf: itran outside range'
  return
endelse
;
; k=ndep: lower boundary
;
sp1(ndep-1)=-1.                        
sp2(ndep-1)=dtauq(ndep-1)+0.5*dtauq(ndep-1)^2
sp3(ndep-1)=0.
p(ndep-1)=s(ndep-1)*(dtauq(ndep-1)+0.5*dtauq(ndep-1)^2)+ $
 (s(ndep-1)-s(ndep-2))
;
; eliminate subdiagonal
;
for k=0,ndep-2 do begin
  f=-sp1(k+1)/(sp2(k)-sp3(k))
  p(k+1)=p(k+1)+f*p(k)
  sp2(k+1)=sp2(k+1)+f*sp2(k)
  sp2(k)=sp2(k)-sp3(k)
endfor
sp2(ndep-1)=sp2(ndep-1)-sp3(ndep-1)
;
; backsubstitute
;
p(ndep-1)=p(ndep-1)/sp2(ndep-1)
for k=ndep-2,0,-1 do begin
  p(k)=(p(k)-sp3(k)*p(k+1))/sp2(k)
endfor
;
; surface intensity
;
iplus_0=2.*(ex*p(0)+s(0)*0.5*ex1^2)
;
;  intensities, calculated from a weighted average of
;  dp/dtauny(k+0.5) and dp/dtauny(k-0.5)
;
iplus(0)=2.*p(0)-ex1*s(0)
iminus(0)=ex1*s(0)
dpdt2=(p(1)-p(0))/dtauq(1)
for k=1,ndep-2 do begin
  dpdt1=dpdt2
  dpdt2=(p(k+1)-p(k))/dtauq(k+1)
  pprimk=(dtauq(k)*dpdt2+dtauq(k+1)*dpdt1)/(dtauq(k)+dtauq(k+1))
  iplus(k)=p(k)+pprimk
  iminus(k)=p(k)-pprimk
endfor
iplus(ndep-1)=s(ndep-1)+(s(ndep-1)-s(ndep-2))/dtauq(ndep-1)
iminus(ndep-1)=2.0*p(ndep-1)-iplus(ndep-1)
;
; calculate p(k)-s(k)
;
pms(0)=p(0)-s(0)
for k=1,ndep-2 do begin
  if(abs(a1(k)) gt 1.0) then pms(k)=p(k)-s(k) else $
    pms(k)=c1(k)*(p(k+1)-p(k))-a1(k)*(p(k)-p(k-1))
endfor
pms(ndep-1)=(p(ndep-2)-p(ndep-1)+s(ndep-1)-s(ndep-2)) $
 /(dtauq(ndep-1)+0.5*dtauq(ndep-1)^2)
;
return
end
