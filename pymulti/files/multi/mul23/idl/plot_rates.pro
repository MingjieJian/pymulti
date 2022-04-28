pro plot_rates,level,xrange=xrange,yrange=yrange,min_rate=min_rate,text=text,$
               labn=labn,norm=norm,xscale=xscale,tauplot=tauplot,nocalc=nocalc,$
               include=include,nolab_header=nolab_header,$
               x0=x0,y0=y0,dy=dy,csize=csize,_extra=e
;+
;   plot_rates,level,xrange=xrange
;
;            plots statistical equilibrium for a given level
;            level    level to consider
;            xrange   range in x
;            yrange   range in y
;            min_rate minimum rate to be plotted
;            text     gives labels, default is nl
;            /labn    gives labels=indgen(nk)
;            norm     gives normalization of rate by n
;            xscale   give x-axis scale : 'tau','cmass',or 'height'
;            tauplot  plot optical depth scale
;            nocalc   uses the already calculated net rates pr and pc
;            x0,y0,dy,csize parameters for label placement
;            include  levels to show
;            no_labheader  gives no label header
;-
@common_multi

if(n_params(0) lt 1) then begin
  print,'plot_rates,level,xrange=xrange,yrange=yrange,min_rate=min_rate,$'
  print,'           text=text,/labn,/norm,xscale=xscale,/tauplot,/nocalc,$'
  print,'           x0=x0,y0=y0,dy=dy,include=include,/nolab_header'
  return
endif

; defaults
;
;  30-Nov-1992 PGJ modifications
;
xtitle='lg !4s!3!d500!n'
zz=taulg
if(keyword_set(xscale)) then begin
   case strupcase(xscale) of
     'HEIGHT': begin
        zz=height
        xtitle='Height [km]'
      end
     'CMASS': begin
        zz=alog10(cmass)
        xtitle='lg m [g cm!U-2!N]'
      end
     'TAU': begin
        xtitle='lg !4s!3!d500!n'
      end
   else: begin
        xtitle='lg !4s!3!d500!n'
      end
   endcase
endif
;
;
if(n_elements(xrange) eq 0) then xrange=[min(zz),max(zz)]
minxr=min(xrange)
maxxr=max(xrange)
iw=where((zz ge minxr) and (zz le maxxr), count)
if(count lt 1) then begin
  print,'No data in range'
  return
endif
;
;  30-Nov-1992 PGJ modifications end
;
if(n_elements(min_rate) eq 0) then min_rate=0.05
if(n_elements(labn) ne 0) then begin
  text=strtrim(string(indgen(nk)),2)
endif
if(n_elements(text) eq 0) then begin
  c_n,label,nn
  c_l,label,ll
  text=strarr(nk)
  for i=0,nk-2 do text(i)=strtrim(nn(i),2)+strmid('spdfghikl',ll(i),1)
  text(nk-1)='cont'
endif
if(n_elements(norm) eq 0) then norm=0
if(n_elements(nocalc) eq 0) then nocalc=0
if(n_elements(include) eq 0) then include=indgen(nk)

; check input

if(level lt 0) or (level gt nk-1) then begin
  print,'level out of range, level=',level
  return
endif


; calculate rates per particle

if(nocalc eq 0) or (n_elements(pr) eq 0) or (n_elements(pc) eq 0) then begin
  net_rates

; multiply with number of particles if norm=0

  if(norm eq 0) then begin
    for k=0,ndep-1 do begin
      pc(k,*,*)=pc(k,*,*)*totn(k)
      pr(k,*,*)=pr(k,*,*)*totn(k)
    endfor
  endif
endif

; find maximum value of any net rate

max_abs_rate=max(abs(pr(iw,*,level))) > max(abs(pc(iw,*,level)))
ymin=min(pr(iw,*,level)) < min(pc(iw,*,level))
ymax=max(pr(iw,*,level)) > max(pc(iw,*,level))

; find total sum of all net rates

total_sum=fltarr(ndep)
for i=0,nk-1 do begin
  total_sum=total_sum+pr(*,i,level)+pc(*,i,level)
endfor

ymin=ymin < min(total_sum(iw))
ymax=ymax > max(total_sum(iw))

if(n_elements(yrange) eq 0) then $
  if(max(abs(!y.range)) eq 0) then yrange=[ymin,ymax] else yrange=!y.range
if(norm eq 0) then ytitle='dn/dt' else ytitle='1/n dn/dt'
if(!y.title eq ' ') then ytitle=''
plot,[xrange(0),xrange(1)],[ymin,ymax],/nodata,xrange=xrange,$
 xtitle=xtitle,ytitle=ytitle,yrange=yrange,$
 title=label(level) + '  ' + strtrim(strmid(atmoid,0,25),2),$
 xstyle=1,ystyle=3,_extra=e
;
;  30-Nov-1992 PGJ modifications begin
;
if(keyword_set(tauplot)) then begin
  mint=fix(min(alog10(tauq)))
  xlabel0=-5+[indgen(15)]
  y00=0.15
  scale2,zz,alog10(tauq),y00,1,0.7,xlabel=xlabel0
endif
;
;  30-Nov-1992 PGJ modifications end
;
linestyle=0
sum_rate=fltarr(ndep)   ; sum of all plotted rates
if(n_elements(x0) eq 0) then x0=0.15  ; relative position for first label
if(n_elements(y0) eq 0) then y0=0.94
if(n_elements(dy) eq 0) then dy=0.06
if(n_elements(csize) eq 0) then csize=1.0
if(n_elements(nolab_header) eq 0) then begin
  label_bob,x0,y0,0,'Collisions, no markers: Radiation',psym=4,csize=csize
  y0=y0-dy
endif
linestyle=[0,1,2,3,4,5]
iplot=0
for i=0,nk-1 do begin
  plot_i=0
  dum=where(i eq include,count)
  if(i ne level) and (count gt 0) then begin
    if(max(abs(pr(iw,i,level))) gt min_rate*max_abs_rate) then begin
      oplot,zz(iw),pr(iw,i,level),linestyle=linestyle(iplot)
      sum_rate=sum_rate+pr(*,i,level)
      plot_i=1
    endif
    if(max(abs(pc(iw,i,level))) gt min_rate*max_abs_rate) then begin
      oplot,zz(iw),pc(iw,i,level),linestyle=linestyle(iplot)
      oplot,zz(iw),pc(iw,i,level),psym=4
      sum_rate=sum_rate+pc(*,i,level)
      plot_i=1
    endif
    if(plot_i eq 1) then begin
      print,label(i),iplot
      textlabel=text(i)+' -> '+text(level)
      label_bob,x0,y0,linestyle(iplot),textlabel,csize=csize
      y0=y0-dy
      iplot=iplot+1
      if(iplot gt 5) then iplot=0
    endif
  endif
endfor

oplot,zz(iw),sum_rate(iw),psym=1
label_bob,x0,y0,0,'Sum',psym=1,csize=csize

max_dev=max(abs(sum_rate-total_sum))
if(min_rate ne 0.0) then print,'max deviation (sum - total_net_rate) =',max_dev
return
end

