pro rn3d,file ,verbose=verbose
;+
;   rn3d,file ,/verbose
;
;            reads N3D file from mul23guc_3d
;-
@common_multi3d
common cn3d,n3d

if(n_params() lt 1) then begin
  print,'rn3d,file ,/verbose'
  return
endif

openr,lur,file,/get_lun

stat=fstat(lur)   ; get statistics of file
nk=stat.size/nx/ny/ndep/4
if(keyword_set(verbose)) then print,'allocating n3d array'
n3d=fltarr(nx,ny,ndep,nk)
rec=assoc(lur,fltarr(nx,ny,ndep))
itype=4L
isize=long(nx)*long(ny)*long(ndep)
for i=0,nk-1 do begin
  cname='nk i='+string(i+1,format='(i2)')
  if(keyword_set(verbose)) then print,itype,isize,' ',cname
  n3d[*,*,*,i]=rec[i]
endfor

free_lun,lur
end

common cn3d,n3d
end
