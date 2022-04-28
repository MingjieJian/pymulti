pro intt_dt,id
;+
;   intt_dt,id
;
;             read idl1 files for response function calculations
;             files are idl1.id_t000...idl1.id_tndep
;             results saved in intt_id_dt.idl
;             assumes all transitions are iwide=.true.
;-
@common_multi

if(n_params() lt 1) then begin
  print,'intt_dt,id'
  return
endif

f=file_search('idl1.'+id+'_t*',count=nfiles)
if(nfiles eq 0) then begin
  print,'files '+'idl1.'+id+'_t* not found'
  return
endif

mulrd,'idl1.'+id+'_t000'
if(ndep+1 ne nfiles) then begin
  print,'nfiles ne ndep+1',nfiles,ndep+1
  return
endif
nqq=max(nq)
intt=fltarr(nqq,ndep+1,nrad)
for i=0,nfiles-1 do begin
  mulrd,'idl1.'+id+'_t'+string3(i)
  intt[*,i,*]=outint[1:nqq,nmu-1,0:nrad-1]
endfor
xl=cc*1.d8/frq[1:nqq,*]
xx=[height,2*height[ndep-1]-height[ndep-2]]
save_file='intt_'+id+'_dt.idlsave'
save,xx,rho,xl,nq,intt,file=save_file
Print,'saved in '+save_file

end
