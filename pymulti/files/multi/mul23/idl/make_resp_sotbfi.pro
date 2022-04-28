pro make_resp_sotbfi,id
;+
;   make_resp_sotbfi,id
;
;            make response functions for SOT WB filters
;+
@common_multi

if(n_params() lt 1) then begin
  print,'make_resp_sotbfi,id'
  return
endif

restore,'intt_'+id+'_dt.idlsave',/verb

n=n_elements(xx)
sotbfii=fltarr(6,n)
for k=0,n-1 do begin
  calc_sotbfi,xl,nq,reform(intt[*,k,*]),yy
  sotbfii[*,k]=yy
endfor

sotbfir=fltarr(n,6)
for i=1,5 do begin
  yy=reform(sotbfii[i,*]/sotbfii[i,n-1]-1.)
  sotbfir[*,i]=deriv(xx,yy)*100.  ; 1% perturbation in delta(t)/t
endfor

save_file='resp_'+id+'.idlsave'
save,xx,rho,sotbfir,file=save_file
print,'saved in '+save_file

end
