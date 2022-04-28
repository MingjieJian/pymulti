pro watmos_dt,file ,dt=dt,dk=dk
;+
;   watmos_dt,file ,dt=dt,dk=dk
;
;            produces atmos files with temperature perturbation
;            file  input atmosphere file is file
;            dt relative perturbation, defaults to 0.01
;            output written to files file_txxx
;            dk step in depth, defaults to 1
;
;-
@common_multi

if(n_params() lt 1) then begin
  print,'watmos_dt,file ,dt=dt,dk=dk'
  return
endif

if(n_elements(dt) eq 0) then dt=0.01
if(n_elements(dk) eq 0) then dk=1

watmos,file+'_t'+string3(ndep)
for k=ndep-1,0,-dk do begin
  temp[k]=temp[k]*(1.+dt)
  watmos,file+'_t'+string3(k)
endfor

end

