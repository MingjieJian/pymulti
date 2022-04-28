pro plot_response_sotbfi,id ,atmosid=atmosid,xrange=xrange,ps=ps,gif=gif

if(n_params() lt 1) then begin
  print,'plot_response_sotbfi,id ,atmosid=atmosid,xrange=xrange,ps=ps,gif=gif'
  return
endif

if(n_elements(atmosid) eq 0) then atmosid=id
if(keyword_set(gif)) then begin
  set_plot,'z'
  device,set_res=[960,600]
endif else if(keyword_set(ps)) then begin
  set_plot,'ps'
  file_ps='response_'+id+'.ps'
  device,file=file_ps,/color,bits=8,/portrait
endif

if(n_elements(xrange) eq 0) then xrange=[-100,1300]

restore,'resp_'+id+'.idlsave'

colf=[0,0,254,200,110,1]
dx=([0,0,10,-20,0,20]+10)*abs(xrange[1]-xrange[0])/1400.
spec_ct,r,g,b
!p.charsize=1.5
sotbfir=sotbfir > 0.
plot,xx,sotbfir[*,3],back=255,col=0,xrange=xrange,xst=1,/nodata,$
 xtitle='Height [km]',ytitle='Response function to !4D!3I/I given !4D!3T/T'
for i=1,5 do oplot,xx,sotbfir[*,i],col=colf[i]
l_sty=[0,0,0,0,0]
text=['Ca H','G-band','Blue','Green','Red']
label,l_sty,text,x0=0.7,y0=0.8,size=1.5,color=colf[1:*]
xyouts,0.25,0.87,'Background atmosphere:'+strupcase(atmosid),/norm,col=0,size=!p.charsize
for i=1,5 do begin
  yy=sotbfir[*,i]
  hresp=trapez(xx,xx*yy)/trapez(xx,yy)
  oplot,[1,1]*hresp,[0,0.005],col=colf[i]
  xyouts,hresp+dx[i],0.006,string(fix(hresp),format='(i4)'),orient=90.,col=colf[i],size=!p.charsize
endfor

text=''
if(keyword_set(gif)) then begin
  file_gif='response_sotbfi_'+id+'.gif'
  write_gif,file_gif,tvrd(),r,g,b
  print,'gif file: '+file_gif
endif else if(keyword_set(ps)) then begin
  device,/close
  print,'ps file:',file_ps
  file_ps='response_'+id+'_ad.ps'
  device,file=file_ps
endif else read,'<cr> to continue',text

n=n_elements(xx)
scale=[1,sqrt(rho[n-2]/rho)]
plot,xx,sotbfir[*,3]*scale,back=255,col=0,xrange=xrange,xst=1,/nodata,$
 xtitle='Height [km]',ytitle='Response function to !4D!3I/I given !4D!3T/T scaled with 1/sqrt(rho)'
for i=1,5 do oplot,xx,sotbfir[*,i]*scale,col=colf[i]
l_sty=[0,0,0,0,0]
text=['Ca H','G-band','Blue','Green','Red']
label,l_sty,text,x0=0.7,y0=0.8,size=1.5,color=colf[1:*]
xyouts,0.25,0.87,'Background atmosphere:'+strupcase(atmosid),/norm,col=0,size=!p.charsize
for i=1,5 do begin
  yy=sotbfir[*,i]*scale
  hresp=trapez(xx,xx*yy)/trapez(xx,yy)
  oplot,[1,1]*hresp,[0,0.005],col=colf[i]
  xyouts,hresp+dx[i],0.006,string(fix(hresp),format='(i4)'),orient=90.,col=colf[i],size=!p.charsize
endfor

if(keyword_set(gif)) then begin
  file_gif='response_sotbfi_'+id+'_ad.gif'
  write_gif,file_gif,tvrd(),r,g,b
  print,'gif file: '+file_gif
endif else if(keyword_set(ps)) then begin
  device,/close
  print,'ps file:',file_ps
endif

set_plot,'x'
end
