pro label_bob,x,y,line,text,psym=psym,csize=csize

if(x lt -1.) or (y lt -1.) then return   ; no label if x or y lt -1

; Label a plot at rel coordinates x,y with line type and text
; Have to do this way because oplot ignores /normal
      if n_elements(psym) eq 0 then psym=0

	x0=!x.crange(0)+(x+0.02)*(!x.crange(1)-!x.crange(0))
	x1=!x.crange(0)+x*(!x.crange(1)-!x.crange(0))
	x2=!x.crange(0)+(x-0.1)*(!x.crange(1)-!x.crange(0))
	if (!x.type eq 1) then begin
	  x0=exp(2.3026*x0)
	  x1=exp(2.3026*x1)
	  x2=exp(2.3026*x2)
	endif
	y1=!y.crange(0)+y*(!y.crange(1)-!y.crange(0))
	y2=!y.crange(0)+(y+0.008)*(!y.crange(1)-!y.crange(0))
	if (!y.type eq 1) then begin
	  y1=exp(2.3026*y1)
	  y2=exp(2.3026*y2)
	endif
	xyouts,x0,y1,text,size=csize
	oplot,[x1,.5*(x1+x2),x2],[y2,y2,y2],line=line,psym=psym
end
