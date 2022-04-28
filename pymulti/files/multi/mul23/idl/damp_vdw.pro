function damp_vdw,gvlg,zz,evi,evj,id
;+
;  gw=damp_vdw(gvlg,zz,evi,evj,id)
;
;            returns MULTI format van der Waals broadening parameter
;            given VALD lg(g.vdw)
;-
if(n_params() lt 5) then begin
  print,'gw=damp_vdw(gvlg,zz,evi,evj,id)'
  return,0
endif

idel=['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si',$
      'P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni']
awel=[1.0,4.0,6.9,9.0,10.8,12.0,14.0,16.0,19.0,20.2,23.0,24.3,27.0,28.1,$
 31.0,32.1,35.5,39.9,39.1,40.1,45.0,47.9,50.9,52.0,54.9,55.8,58.9,58.7]
evcel=fltarr(2,28)
evcel[0,*]=[13.6,24.6,5.4,9.3,8.3,11.3,14.5,13.6,17.4,21.6,5.1,7.6,6.0,8.1,10.5,$
 10.4,13.0,15.8,4.3,6.1,6.5,6.8,6.7,6.8,7.4,7.9,7.9,7.6]
evcel[1,*]=[100.,54.4,75.6,18.2,25.1,24.4,29.6,35.1,35.0,41.1,47.3,15.0,18.8,16.3,$
 19.7,23.4,23.8,27.6,31.8,11.9,12.8,13.6,14.7,16.5,15.6,16.2,17.1,18.2]
iel=where(strtrim(strupcase(id),2) eq strtrim(strupcase(idel),2),count)
if(count ne 1) then begin
  print,id,' not found among ',idel
  return,0
endif
evc=evcel[zz-1,iel]
awgt=awel[iel]

uu=1.66057e-24
bk=1.38066e-16

temp=10000.

c625=1.283984e-12*zz^0.8*(1./(evc-evj)^2-1./(evc-evi)^2)^.4
gv=8.411*(8.*bk*temp/!pi*(1./(1.008*uu)+1./(awgt*uu)))^0.3*c625
gw=10.^gvlg/gv

return,gw
end
