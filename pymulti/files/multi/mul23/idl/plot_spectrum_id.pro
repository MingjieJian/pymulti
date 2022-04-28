pro plot_spectrum_id ,kr=kr,ps=ps,gif=gif,xrange=xrange,yrange=yrange,color_sp=color_sp,limit=limit,image=image,$
 xlv=xlv,int0=int0,atmos=atmos,vmac=vmac,nne0=nne0,temp0=temp0,reverse=reverse,$
 ftsatlas=ftsatlas,nograv=nograv,abslin_file=abslin_file,even_dx=even_dx,dv0=dv0,$
 relabn=relabn,debug=debug,airxl=airxl,icont=icont,trady=trady,p_sym=p_sym,_extra=e,labaa=labaa
;+
; NAME:
;       PLOT_SPECTRUM_ID
;
; PROJECT:
;       La Palma, Multi
;
; PURPOSE:
;       Plots spectrum with line identifications
;
; CATEGORY:
;       Plotting
;
; SYNTAX:
;	plot_spectrum_id ,kr=kr,ps=ps,gif=gif,xrange=xrange,limit=limit,$
;       xl=xl,int=int,atmos=atmos,vmac=vmac,nne0=nne0,temp0=temp0,$
;       /ftsatlas,/nograv,abslin_file=abslin_file,/relabn,/airxl,/even_dx
;
; INPUTS:
;       only through keywords. Either KR or XL,INT or /FTSATLAS has to be given
;
; OUTPUTS:
;       Plots spectrum with line identifications
;
; OPTIONAL OUTPUTS:
;       None.
;
; KEYWORDS:
;       KR         - line number for Multi transition to plot
;       XLV        - wavelength in vaccum for spectrum to plot [nm]
;       INT0       - intensity to plot
;       AIRXL      - wavelength given in xl is air wavelength [nm]
;       ATMOS      - atmospheric ID string
;	FTSATLAS   - if set, solar atlas is plotted
;       NOGRAV     - do not apply gravitational redshift to ftsatlas, defaults to false
;       ICONT      - continuum value for ftsatlas, defaults to max(INT)
;       TRADY      - plot with y-scale in radiation temperature 
;       XRANGE     - plotting range in nm, defaults to range in Multi or XL if given
;       YRANGE     - plotting range in int, defaults to range in int
;       COLOR_SP   - color index for spectrum plot, defaults to 0
;       REVERSE    - plot decreasing wavelength if xrange is not given
;       ABSLIN_FILE- file with ABSLIN data, defaults to ~/multi/run/ABSLIN
;       VMAC       - macroturbulence to convolve spectrum with, in km/s, defaults to 0
;	NNE0       - Ne for Saha/Boltzmann calculation. Defaults to solar value at tau_500=0.1
;	TEMP0      - T  for Saha/Boltzmann calculation. Defaults to solar value at tau_500=0.1
;       LIMIT      - limit of f*N(i) compared with 1. for inclusion,
;                    defaults to 1.e-5
;       EVEN_DX    - evenly spaced line-ID labels over plot, default to false
;       DV0        - wavelength for line center from which to
;                    calculate velocities
;       RELABN     - use abundance relative to maximum value in ABSLIN
;       IMAGE      - image to overplot spectrum on
;       PS         - name of Postscript file to receive plot. /PS gives fig_spectrum_id.ps
;       GIF        - name of Gif file to receive plot. /GIF gives
;                    fig_spectrum_id.gif
;       LABAA      - labels in Aangstrom
;
; COMMON:
;       Common_multi - multi common blocks
;
; RESTRICTIONS:
;       None.
;
; SIDE EFFECTS:
;       None.
;
; CALLS
;	rabslin,cstrip,convl,iondata,saha,interpol,label
;
; HISTORY:
;       Version 1, April 28 2006, Mats Carlsson, UiO. Written
;
; CONTACT:
;       Mats Carlsson (mats.carlsson@astro.uio.no)
;-
@common_multi
common cfts_atlas,xl_fts,int_fts

if(!p.charsize eq 0.0) then pcharsize=1.0 else pcharsize=!p.charsize
cc=2.9979247e10
hh=6.62618e-27
bk=1.38066e-16
if(n_elements(nne0) eq 0) then nne0=7.19e12  ; from tau_500=0.1 in VAL3C
if(n_elements(temp0) eq 0) then temp0=5238.  ; from tau_500=0.1 in VAL3C
if(n_elements(abslin_file) eq 0) then abslin_file='ABSLIN'
if(n_elements(color_sp) eq 0) then color_sp=0
atlas_only=0
if(n_elements(int0) ne 0) then int=int0
if(n_elements(xlv) eq 0) then begin       ; no xl given, Multi or only ftsatlas
  if(n_elements(kr) ne 0) then begin     ; kr given, Multi spectrum to plot
    if(n_elements(atmos) eq 0) then begin
      dum=strsplit(def_ext,'_',/extract)
      atmos=dum[n_elements(dum)-1]
    endif
    nwide=0
    for i=0,nline-1 do begin
      if(iwide[i]) then nwide=nwide+1
    endfor
    if(kr ge nline) then begin
      xl_vac=cc*1.e8/frq[1:nq[kr],kr-nline+nwide]*0.1   ; xl in nm
    endif else begin
      double,kr,outint,xx,yy
      xl_vac=(alamb[kr]+xx)*0.1
      int=yy
    endelse
    xl_air=convl(xl_vac*10.)*0.1
    yscale=1.e-3
    ytitle='Intensity [J m!u-2!n s!u-1!n Hz!u-1!n sr!u-1!n]'
    if(keyword_set(trady)) then ytitle='Trad [K]'
  endif else begin
    atlas_only=1
    yscale=1.0
    ytitle='Intensity'
  endelse
endif else begin  ; xlv given
  if(n_elements(atmos) eq 0) then atmos=''
  xl_vac=xlv
  if(not keyword_set(airxl)) then xl_air=convl(xl_vac*10.)*0.1 else begin
    xl_air=xl_vac
    if(keyword_set(trady)) then begin  ; need vacuum wavelengths for trady option
      for i=0,n_elements(xl_vac)-1 do xl_vac[i]=inv_convl(xl_air[i]*10.)*0.1
    endif
  endelse
  yscale=1.0
  if(keyword_set(trady)) then ytitle='Trad [K]' else ytitle='Intensity'
endelse
thick=1
nspace=18
dnm=0.7
xsize=28
ysize=17.78
overlap=0.1
if(n_elements(ps) ne 0) then begin         ; set postscript options
  set_plot,'ps'
  file_ps=strtrim(string(ps),2)
  if(file_ps eq '1') then begin
    file_ps='fig_spectrum_id.ps'
  endif
  thick=3
  nspace=20
  device,file=file_ps,xsize=24,ysize=16,/landscape,/color,bits=8
endif
if(n_elements(gif) ne 0) then begin        ; set gif-file options
  if(n_elements(ps) ne 0) then begin
    print,'/ps and /gif cannot both be given'
    return
  endif
  set_plot,'z'
  device,set_resolution=[1000,600]
  file_gif=strtrim(string(gif),2)
  if(file_gif eq '1') then begin
    file_gif='fig_spectrum_id.gif'
  endif
endif

; get intensity and wavelength scales

if(not atlas_only) then begin
  xl0=min(xl_air)
  xl1=max(xl_air)
  if(n_elements(xrange) eq 0) then begin
    xrange=[xl0,xl1]
    if(keyword_set(reverse)) then xrange=reverse(xrange)
  endif
  if(n_elements(int) eq 0) then begin
    int=outint[1:nq[kr],nmu-1,kr]
  endif
  if(n_elements(icont) eq 0) then icont=max(int)
endif else begin
  if(n_elements(xrange) eq 0) then begin
    print,'with /ftsatlas xrange (in nm) has to be given'
    return
  endif
  xl0=min(xrange)
  xl1=max(xrange)
endelse
if(keyword_set(ftsatlas)) then begin
  if(n_elements(xl_fts) eq 0) then begin
    restore,'../idl/ftsatlas.idlsave'
  endif
  iw=where((xl_fts gt xl0) and (xl_fts lt xl1),count)
  if(count eq 0) then begin
    print,'FTS atlas only for range 329 to 1251 nm'
    stop
  endif
  xl_solar=xl_fts[iw]
  if(not keyword_set(nograv)) then begin
    xl_solar=xl_solar*(1.d0-633.d2/cc)  ; gravitational redshift subtracted
  endif
  int_solar=int_fts[iw]
endif

; macroturbulence convolution, assumes equidistant wavelengths!!

if(not atlas_only) then begin
  if(n_elements(vmac) eq 0) then vmac=0.0
  if(vmac gt 0.) then begin
    xl=xl_air
    xk = xl[0:30] - mean(xl[0:30]) ; kernel
    ld = vmac*1.e5/cc*mean(xl)           ;delta lambda vmac
    yk = exp(-(xk/ld)^2)            ;kernel gaussian
    yk = yk/total(yk)
    ymac = convol(int, yk, /center) ;convolution
    int=ymac
  endif
endif

; find lines from ABSLIN and abundances from Saha/Boltzmann
; Include lines with relative abundance gt limit

dum=file_search(abslin_file,count=count)
if(count gt 0) then begin
  if(n_elements(limit) eq 0) then limit=1.e-5
  rabslin,lines,file=abslin_file
  nlines=n_elements(lines)
  xll_air=convl(reform(lines[*].alamll))*0.1
  abn=fltarr(nlines)
  for i=0,nlines-1 do begin ; go through lines to calculate LTE population times f
    id=strmid(lines[i].abnill,0,2)
    iondata,id,awgt,abund,xi,u0,u1
    if(abund gt 0.0) then begin
      saha,temp0,nne0,xi,u0,u1,n1_ntot,n0_ntot
      if(lines[i].ionill eq 1) then begin
        saha_frac=n0_ntot
        bolz_frac=lines[i].gill/u0
      endif else begin
        saha_frac=n1_ntot
        bolz_frac=lines[i].gill/u1
      endelse
      bolz_frac=bolz_frac*exp(-lines[i].chill*cc*hh/bk/temp0)
      abn[i]=abund*saha_frac*bolz_frac*lines[i].fll
    endif else begin
      abn[i]=0.0
    endelse
  endfor
endif else begin
  nlines=0
  message,abslin_file+' not found, no line-ids plotted',/info
endelse

if(keyword_set(debug)) then begin
  openw,luw,'lines.tmp',/get_lun
  printf,luw,'  abund     saha_frac bolz_frac f         abn       abn/max(abn)'
  for i=0,nlines-1 do begin ; go through lines to calculate LTE population times f
    id=strmid(lines[i].abnill,0,2)
    iondata,id,awgt,abund,xi,u0,u1
    saha,temp0,nne0,xi,u0,u1,n1_ntot,n0_ntot
    if(lines[i].ionill eq 1) then begin
      saha_frac=n0_ntot
      bolz_frac=lines[i].gill/u0
    endif else begin
      saha_frac=n1_ntot
      bolz_frac=lines[i].gill/u1
    endelse
    bolz_frac=bolz_frac*exp(-lines[i].chill*cc*hh/bk/temp0)
    printf,luw,id,' ',lines[i].ionill,xll_air[i],format='(a2,a1,i1,f10.3)'
    if(keyword_set(relabn)) then rel_abn=abn[i]/max(abn) else rel_abn=abn[i]/4e-10
    printf,luw,format='(6e10.2)',abund,saha_frac,bolz_frac,lines[i].fll,abn[i],rel_abn
  endfor
  free_lun,luw
endif

if(nlines gt 0) then begin
  if(keyword_set(relabn)) then begin
    abn=abn/max(abn)
  endif else begin
    abn=abn/4e-10
  endelse
  abn=abn < 1.0

  iwl=where((abn gt limit) and (xll_air gt min(xrange)) and (xll_air lt max(xrange)),nlines)
endif

; load colour lookup table

tvlct,r0,g0,b0,/get
r=indgen(256)
g=r
b=r
r[1]=255
tvlct,r,g,b

if(not atlas_only) then begin

; plot multi or given spectrum

  if(keyword_set(trady)) then y=tradb(int,xl_vac*10.) else y=int*yscale
  if(n_elements(dv0) eq 0) then begin
    xx=xl_air
    xtitle='Wavelength [nm]'
  endif else begin
    xx=(xl_air-dv0)/dv0*300000.
    xrange=(xrange-dv0)/dv0*300000.
    xtitle='Velocity [km/s]'
  endelse
  if(n_elements(image) ne 0) then begin
    size_image=size(image)
    ny_image=size_image[2]
    bimage=bytscl(image,top=253)+2
    if(n_elements(yrange) eq 0) then yrange=[0,max(y)]
    yy=findgen(ny_image)/(ny_image-1)*(yrange[1]-yrange[0])+yrange[0]
    mplot_image,bimage,xx,yy,norm=0,/int,xtitle=xtitle,$
     ytitle=ytitle,$
     yrange=yrange,col=0,back=255,$
     xrange=xrange,xstyle=1,thick=thick,xmargin=[12,4],ymargin=[4,15],_extra=e
    oplot,xl_air,y,col=color_sp
  endif else begin
    plot,xx,y,xtitle=xtitle,$
     ytitle=ytitle,$
     yrange=yrange,col=0,back=255,$
     xrange=xrange,xstyle=1,thick=thick,xmargin=[12,4],ymargin=[4,15],_extra=e,/nodata
    oplot,xx,y,col=color_sp,psym=p_sym
  endelse

; overplot FTS spectrum

  if(keyword_set(ftsatlas)) then begin
    if(keyword_set(trady)) then begin
      xl_solar_vac=xl_solar
      for i=0,n_elements(xl_solar)-1 do xl_solar_vac[i]=inv_convl(xl_solar[i]*10.)*0.1
      y=tradb(int_solar*icont,xl_solar_vac*10.) 
    endif else y=int_solar*yscale*icont
    if(n_elements(dv0) eq 0) then begin
      xx=xl_solar
    endif else begin
      xx=(xl_solar-dv0)/dv0*300000.
    endelse
    oplot,xx,y,col=1
    label,[0,0],[atmos,'FTS solar atlas'],color=[0,1],x0=0.02,y0=0.10,dx=0.05
  endif else begin
    label,[0],[atmos],color=[0],x0=0.02,y0=0.10,dx=0.05
  endelse
endif else begin
  if(n_elements(dv0) eq 0) then begin
    xx=xl_solar
    xtitle='Wavelength [nm]'
  endif else begin
    xx=(xl_solar-dv0)/dv0*300000.
    xrange=(xrange-dv0)/dv0*300000.
    xtitle='Velocity [km/s]'
  endelse
  if(n_elements(image) ne 0) then begin
    size_image=size(image)
    ny_image=size_image[2]
    bimage=bytscl(image,top=253)+2
    yy=findgen(ny_image)/(ny_image-1)*1.05
    mplot_image,bimage,[min(xrange),max(xrange)],yy,norm=0,xtitle=xtitle,$
     ytitle=ytitle,$
     yrange=[0,1.05],col=0,back=255,$
     xrange=xrange,xstyle=1,ystyle=1,thick=thick,xmargin=[12,4],ymargin=[4,15],_extra=e
    oplot,xx,int_solar,col=1
    label,[0],['FTS solar atlas'],color=[1],x0=0.02,y0=0.10,dx=0.05
  endif else begin
    plot,xx,int_solar,xtitle=xtitle,$
     ytitle=ytitle,$
     yrange=[0,1.05],col=0,back=255,$
     xrange=xrange,xstyle=1,ystyle=1,thick=thick,xmargin=[12,4],ymargin=[4,15],/nodata
     oplot,xx,int_solar,col=1
    label,[0],['FTS solar atlas'],color=[1],x0=0.02,y0=0.10,dx=0.05
  endelse
endelse

; write line identifications

y0=!y.crange[1]
y1=y0+0.02*(!y.crange[1]-!y.crange[0])
y2=y0+0.10*(!y.crange[1]-!y.crange[0])
y3=y0+0.13*(!y.crange[1]-!y.crange[0])
mindx=float(!d.y_ch_size)/float(!d.x_size)*1.1
for i=0,nlines-1 do begin
  limit0=-alog10(limit)+0.1
  y=!y.crange[1]-0.5*(alog10(abn[iwl[i]])+limit0)/limit0*(!y.crange[1]-!y.crange[0])
  x=xll_air(iwl[i])
  if(n_elements(dv0) ne 0) then x=(x-dv0)/dv0*300000.
  if(keyword_set(even_dx)) then begin
    x2=min(xrange)+float(i)/(nlines-1)*abs(xrange[1]-xrange[0])
  endif else begin
    if(i eq 0) then x2old=-1.e8
    if(x-x2old gt mindx*abs(xrange[1]-xrange[0])) then begin
      x2=x
      x2old=x2
    endif else begin
      x2=abs(xrange[1]-xrange[0])*mindx+x2old
      x2old=x2
    endelse
  endelse
  plots,[x,x],[y0,y1],col=0
  plots,[x,x2],[y1,y2],col=0
  plots,[x2,x2],[y2,y3],col=0
  plots,[x,x],[y0,y],col=0
  if(x lt 1000.) then begin
    lspace=float(!d.x_ch_size)/float(!d.y_size)*(!y.crange[1]-!y.crange[0])*nspace*pcharsize*pcharsize
    xyouts,x2,y3+0.1*(y3-y0),lines[iwl[i]].abnill,$
     orient=90,col=0
    if(keyword_set(labaa)) then begin
      xyouts,x2,y3+0.1*(y3-y0)+lspace,string(10.d0*xll_air[iwl[i]],format='(f9.3)'),$
       orient=90,col=0,align=1.0
    endif else begin
      xyouts,x2,y3+0.1*(y3-y0)+lspace,string(xll_air[iwl[i]],format='(f8.3)'),$
       orient=90,col=0,align=1.0
    endelse
  endif else begin
    lspace=float(!d.x_ch_size)/float(!d.y_size)*(!y.crange[1]-!y.crange[0])*(nspace+1)*pcharsize*pcharsize
    xyouts,x2,y3+0.1*(y3-y0),lines[iwl[i]].abnill,$
     orient=90,col=0
    if(keyword_set(labaa)) then begin
      xyouts,x2,y3+0.1*(y3-y0)+lspace,string(10.d0*xll_air[iwl[i]],format='(f10.3)'),$
       orient=90,col=0,align=1.0
    endif else begin
      xyouts,x2,y3+0.1*(y3-y0)+lspace,string(xll_air[iwl[i]],format='(f9.3)'),$
       orient=90,col=0,align=1.0
    endelse
  endelse
endfor

tvlct,r0,g0,b0

if(n_elements(gif) ne 0) then begin
  write_gif,file_gif,tvrd(),r,g,b
  print,'figure file is:',file_gif
  set_plot,'x'
endif
if(n_elements(ps) ne 0) then begin
  print,'figure postscript file is: ',file_ps
  device,/close
  set_plot,'x'
endif

end

