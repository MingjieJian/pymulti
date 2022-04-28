pro watmos_3d21d,ext ,ixx=ixx,iyy=iyy,dscale=dscale,nh_in=nh_in,rstrt=rstrt,intep=intep
;+
;   watmos_3d21d,ext ,dscale=dscale,nh_in=nh_in
;
;            writes multi input files from idl variables
;            output on atmos.ext_ix_iy, dscale.ext_ix_iy
;            where ext is extension (string)
;            optional dscale in km in dscale variable
;            nh is written assuming LTE if nh is not given explicitly
;            if ix,iy not given, whole range is written
;            if /rstrt given, there will also be a rstrt.ext_ix_iy
;            file written
;            /intep gives interpolation to dscale the same way as in multi_3d
;-
@common_multi3d
common cout3d,id3d,taulg3d,cmasslg3d,xnorm3d,icont,x3d,height3d
common cn3d,n3d

if (n_params(0) lt 1) then begin
  print,'watmos_3d21d,ext ,ix=ix,iy=iy,dscale=dscale,nh_in=nh_in,/rstrt,/intep'
  return
endif

em=9.10953d-28
bk=1.38066d-16
hh=6.62618d-27
ee=1.60219d-12

tab=fltarr(5,ndep)
;grph=2.26653d-24
grph=2.38049d-24     
if(n_elements(ixx) eq 0) then ixx=indgen(nx) else ixx=[ixx]
if(n_elements(iyy) eq 0) then iyy=indgen(ny) else iyy=[iyy]

for jj=0,n_elements(iyy)-1 do begin
  iy=iyy[jj]
  for ii=0,n_elements(ixx)-1 do begin
    ix=ixx[ii]
    cnr=string3(ix)+'_'+string3(iy)
    atm_file='atmos.'+ext+'_'+cnr
    dsc_file='dscale.'+ext+'_'+cnr
    rstrt_file='rstrt.'+ext+'_'+cnr
    atmoid='multi3d '+ext+' (ix,iy)='+string(ix,iy,format="(2i4)")
    openw,luatm,atm_file,/get_lun
    printf,luatm,atmoid
    printf,luatm,'Height scale'
    printf,luatm,' 4.4400'
    printf,luatm,ndep
    if(keyword_set(intep) and n_elements(dscale) ne 0) then begin
      height1=height*1.e-5              ; original height in km
      height2=reform(dscale[ix,iy,*])   ; new height scale
      tmp=alog(reform(temp[ix,iy,*]))
      intep,height1,tmp,height2,yp
      temp2=exp(yp)
      tmp=alog(reform(nne[ix,iy,*]))
      intep,height1,tmp,height2,yp
      nne2=exp(yp)
      tmp=reform(vz[ix,iy,*])
      intep,height1,tmp,height2,yp
      vz2=yp
      tmp=alog(reform(rho[ix,iy,*]))
      intep,height1,tmp,height2,yp
      rho2=exp(yp)
    endif else begin  
      height2=height*1.e-5
      if(height[1] gt height[0]) then height2=-height2
      temp2=reform(temp[ix,iy,*])
      nne2=reform(nne[ix,iy,*])
      vz2=reform(vz[ix,iy,*])
      rho2=reform(rho[ix,iy,*])
    endelse
    tab[0,*]=height2
    tab[1,*]=temp2
    tab[2,*]=nne2
    tab[3,*]=vz2
    printf,luatm,tab
    if(n_elements(nh_in) ne 0) then begin
      printf,luatm,transpose(reform(nh_in[ix,iy,*,*])),format='(6e12.4)'
    endif else begin
      nh_in=fltarr(6,ndep)
      for k=0,NDEP-1 do begin
	phit=(2.d0*!pi*em*bk*temp2[k]/hh/hh)^1.5*exp(-13.6d0*ee/bk/temp2[k])
        ratio=phit/nne2[k]
        totnh=rho2[k]/grph
        nh_in[5,k]=ratio/(1.d0+ratio)*totnh
        ratio=nne2[k]/phit
        nh_in[0,k]=ratio/(1.d0+ratio)*totnh
        i=findgen(4)+2
        NH_IN[1:4,K]=i*i*exp(-13.6d0*(1.d0-1.d0/i/i)*ee/bk/temp2[k])*nh_in[0,k]
      endfor
      printf,luatm,nh_in,format='(6e12.4)'
    endelse
    free_lun,luatm
    if(keyword_set(rstrt)) then begin
      openw,lurstrt,rstrt_file,/get_lun
      printf,lurstrt,atmoid
      printf,lurstrt,'Height scale'
      printf,lurstrt,' 4.4400'
      printf,lurstrt,ndep
      tab[0,*]=height3d[ix,iy,*]
      if(height3d[ix,iy,1] gt height3d[ix,iy,0]) then tab[0,*]=-tab[0,*]
      xx=height*1.e-5
      tab[1,*]=interpol(temp[ix,iy,*],xx,height3d[ix,iy,*])
      tab[2,*]=interpol(nne[ix,iy,*],xx,height3d[ix,iy,*])
      tab[3,*]=interpol(vz[ix,iy,*],xx,height3d[ix,iy,*])
      printf,lurstrt,tab
      im=transpose(reform(n3d[ix,iy,*,*]))
      for k=0,ndep-1 do begin
        printf,lurstrt,im[*,k],format='(6e12.4)'
      endfor
      free_lun,lurstrt
    endif
    openw,ludsc,dsc_file,/get_lun
    printf,ludsc,atmoid
    printf,ludsc,'Height scale'
    if(n_elements(dscale) eq 0) then begin
      printf,ludsc,ndep,-12.0
      printf,ludsc,tab[0,*]
    endif else begin
      siz=size(dscale)
      ndep=siz[siz[0]]
      printf,ludsc,ndep,-12.0
      if(siz[0] eq 1) then begin
        for i=0,ndep-1 do begin
          printf,ludsc,dscale[i]
        endfor
      endif else begin
        for i=0,ndep-1 do begin
          printf,ludsc,dscale[ix,iy,i]
        endfor
      endelse
    endelse
    free_lun,ludsc
  endfor
endfor

end


