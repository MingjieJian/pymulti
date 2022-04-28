pro watmos_multi3d,dx,dy,z,tg,nne,rho,vx,vy,vz,file=file ,dp=dp
;+
;   watmos_multi3d,dx,dy,z,tg,nne,rho,vx,vy,vz,file=file ,dp=dp
;
;            writes atmos3d file
;            dx   grid step size in x in cm
;            dy   grid step size in y in cm
;            z    grid           in z in cm
;            tg   temperature in K
;            nne  Ne in cm^-3
;            rho  Mass density in g cm^-3
;            vx   Vx in cm/s [defaults to 0]
;            vy   Vy in cm/s [defaults to 0]
;            vz   Vz in cm/s [defaults to 0], positive is upward velocity
;
;            input cubes should not repeat first and last plane, this
;            is done by this routine such that nx and ny increase by
;            one and tg2[nx-1,*,*]=tg2[0,*,*] and tg2[*,ny-1,*]=tg2[*,0,*]
;
;            z index 0 is at the top of the atmosphere
;
;-
if(n_params() lt 6) then begin
  print,'watmos_multi3d,dx,dy,z,tg,nne,rho ,vx,vy,vz,file=file ,/dp'
  return
endif
if(n_elements(file) eq 0) then file='atmos3d'

; parameter checking

dim=size(tg)
if(dim[0] ne 3) then begin
  message,'dim(tg) ne 3',/informational
  return
endif
nx=long(dim[1])+1
ny=long(dim[2])+1
nz=long(dim[3])

check_dim,nne,tg,'nne','tg'
check_dim,rho,tg,'rho','tg'
if(n_elements(vx) ne 0) then begin
  check_dim,vx,tg,'vx','tg'
  if(keyword_set(dp)) then vx2=dblarr(nx,ny,nz) else vx2=fltarr(nx,ny,nz)
  vx2[0:nx-2,0:ny-2,*]=vx*1.d-5        ; in km/s in atmos3d
  vx2[nx-1,*,*]=vx2[0,*,*]
  vx2[*,ny-1,*]=vx2[*,0,*]
endif else begin
  if(keyword_set(dp)) then vx2=dblarr(nx,ny,nz) else vx2=fltarr(nx,ny,nz)
endelse
if(n_elements(vy) ne 0) then begin
  check_dim,vy,tg,'vy','tg'
  if(keyword_set(dp)) then vy2=dblarr(nx,ny,nz) else vy2=fltarr(nx,ny,nz)
  vy2[0:nx-2,0:ny-2,*]=vy*1.d-5
  vy2[nx-1,*,*]=vy2[0,*,*]
  vy2[*,ny-1,*]=vy2[*,0,*]
endif else begin
  if(keyword_set(dp)) then vy2=dblarr(nx,ny,nz) else vy2=fltarr(nx,ny,nz)
endelse
if(n_elements(vz) ne 0) then begin
  check_dim,vz,tg,'vz','tg'
  if(keyword_set(dp)) then vz2=dblarr(nx,ny,nz) else vz2=fltarr(nx,ny,nz)
  vz2[0:nx-2,0:ny-2,*]=vz*1.d-5
  vz2[nx-1,*,*]=vz2[0,*,*]
  vz2[*,ny-1,*]=vz2[*,0,*]
endif else begin
  if(keyword_set(dp)) then vz2=dblarr(nx,ny,nz) else vz2=fltarr(nx,ny,nz)
endelse

openw,luw,file,/get_lun,/f77_unformatted
itype=3L
isize=3L
cname='dim     '
writeu,luw,itype,isize,cname
writeu,luw,nx,ny,nz

if(keyword_set(dp)) then itype=5L else itype=4L
isize=nx
cname='x grid  '
writeu,luw,itype,isize,cname
if(keyword_set(dp)) then widthx=dindgen(nx)*dx else widthx=findgen(nx)*float(dx)
writeu,luw,widthx

isize=ny
cname='y grid  '
writeu,luw,itype,isize,cname
if(keyword_set(dp)) then widthy=dindgen(ny)*dy else widthy=findgen(ny)*float(dy)
writeu,luw,widthy

isize=nz
cname='z grid  '
writeu,luw,itype,isize,cname
if(keyword_set(dp)) then writeu,luw,double(z) else writeu,luw,float(z) 

isize=nx*ny*nz
cname='nne     '
writeu,luw,itype,isize,cname
if(keyword_set(dp)) then nne2=dblarr(nx,ny,nz) else nne2=fltarr(nx,ny,nz)
nne2[0:nx-2,0:ny-2,*]=nne
nne2[nx-1,*,*]=nne2[0,*,*]
nne2[*,ny-1,*]=nne2[*,0,*]
writeu,luw,nne2

isize=nx*ny*nz
cname='tg      '
writeu,luw,itype,isize,cname
if(keyword_set(dp)) then tg2=dblarr(nx,ny,nz) else tg2=fltarr(nx,ny,nz)
tg2[0:nx-2,0:ny-2,*]=tg
tg2[nx-1,*,*]=tg2[0,*,*]
tg2[*,ny-1,*]=tg2[*,0,*]
writeu,luw,tg2

isize=nx*ny*nz
cname='vx      '
writeu,luw,itype,isize,cname
writeu,luw,vx2

isize=nx*ny*nz
cname='vy      '
writeu,luw,itype,isize,cname
writeu,luw,vy2

isize=nx*ny*nz
cname='vz      '
writeu,luw,itype,isize,cname
writeu,luw,vz2

isize=nx*ny*nz
cname='rho     '
writeu,luw,itype,isize,cname
if(keyword_set(dp)) then rho2=dblarr(nx,ny,nz) else rho2=fltarr(nx,ny,nz)
rho2[0:nx-2,0:ny-2,*]=rho
rho2[nx-1,*,*]=rho2[0,*,*]
rho2[*,ny-1,*]=rho2[*,0,*]
writeu,luw,rho2

free_lun,luw
rho2=0 & tg2=0 & nne2=0 & vx2=0 & vy2=0 & vz2=0

end

pro check_dim,var1,var2,id1,id2
;+
;   check_dim,var1,var2,id1,id2
;
;            check dimensions of var1 and var2
;
;-
if(n_params() lt 4) then begin
  message,'check_dim,var1,var2,id1,id2'
  return
endif
dim1=size(var1)
dim2=size(var2)
if(max(abs(dim1-dim2)) gt 0) then begin
  message,'variable dimension mismatch'
  print,'size('+id1+'=',dim1
  print,'size('+id2+'=',dim2
  return
endif

end
