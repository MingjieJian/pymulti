pro xmovie_ev, event

COMMON cxdata,x1,y1,x2,y2,id_frame
COMMON cxmovie,go,speed,n_frames,frame_nr,frame_increment,$
 movie_base,bg_base,speedslider,frameslider,incrementslider,$
 xminslider,xmaxslider,ixmin,ixmax,xmin,xmax,$
 yminslider,ymaxslider,iymin,iymax,ymin,ymax
COMMON cfixdx,dx_fix

WIDGET_CONTROL, event.id, GET_UVALUE = eventval

CASE eventval OF
  "XMOVIE_DONE" : begin
                    WIDGET_CONTROL, event.top, /DESTROY
                    WIDGET_CONTROL, bg_base, /DESTROY
                  endcase
  "XMOVIE_START": xmovie_start
  "XMOVIE_STOP" : xmovie_stop
  "XMOVIE_PSYM" : xmovie_psym
  "XMOVIE_DX_FIX" : xmovie_dx_fix
  "XMOVIE_HELP" : xdisplayfile,"xmovie.hlp", $
                      TITLE = "Xmovie Help", $
                      GROUP = event.top, $
                      WIDTH = 70, $
                      HEIGHT= 20, $
                      TEXT  = [                                               $
  "Xmovie is used for displaying a series of plots as a time-series.",        $
  "The user can select the speed or specific frames in the animation.",       $
  " ",                                                                        $
  "Xmovie can be called with the following arguments:",                       $
  "  xmovie,y1                 x is taken as the index",                      $
  "  xmovie,x1,y1              x1 is of form (nx) or (nx,n_frames)",          $
  "                            y1 is of form (nx,n_frames) with any",         $
  "                               number of dimensions of size 1",            $
  " ",                                                                        $
  "  xmovie,x1,y1,x2,y2        x2,y2 is overplotted over x1,y1",              $
  "                            with dashed lines",                            $
  " ",                                                                        $
  "keyword ID can be given with an array with frame-id of length n_frames",   $
  "/NORM  normalizes y1 and y2 to range [0,1]",                               $
  " ",                                                                        $
  "The start button will start the animation.",                               $
  "The stop button will pause the animation making possible the view",        $
  "of individual frames.",                                                    $
  " ",                                                                        $
  "The top slider is used to control the speed of",                           $
  "the animation.  Moving it to the far right is one hundred",                $
  "percent, as fast as the animation can go.",                                $
  "When the animation is stopped, the speed slider is inactive.",             $
  "It is then possible to select individual frames with the",                 $
  "Frame number slider. The frame is shown when the mouse-button is ",        $
  "released.",                                                                $
  "The [+] button will move Frame increment frames forward,",                 $
  "the [-] button moves backwards. The increment can be changed",             $
  "with the Frame increment slider."]

  "SPEED_SLIDER": speed=event.value
  "FRAME_SLIDER": begin
                  frame_nr=event.value
                  xmovie_show
                  endcase
  "XMOVIE_FWD"  : xmovie_show,frame_nr+frame_increment
  "XMOVIE_BKW"  : xmovie_show,frame_nr-frame_increment
  "INCR_SLIDER" : frame_increment=event.value
  "XMIN_SLIDER" : begin
                  ixmin=event.value < 999
                  if(dx_fix ne 0) then begin
                    if(ixmin gt 1000-dx_fix) then begin
                      ixmin=ixmin < (1000-dx_fix)
                      WIDGET_CONTROL, xminslider, SET_VALUE=ixmin
                    endif
                    ixmax=(ixmin+dx_fix) < 1000
                    WIDGET_CONTROL, xmaxslider, SET_VALUE=ixmax
                  endif
                  if(ixmin ge ixmax) then begin
                    ixmax=ixmin+1
                    WIDGET_CONTROL, xmaxslider, SET_VALUE=ixmax
                  endif
                  if(go eq 0) then xmovie_show
                  endcase
  "XMAX_SLIDER" : begin
                  ixmax=event.value > 1
                  if(ixmin ge ixmax) then begin
                    ixmin=ixmax-1
                    WIDGET_CONTROL, xminslider, SET_VALUE=ixmin
                  endif
                  if(dx_fix ne 0) then dx_fix=ixmax-ixmin
                  if(go eq 0) then xmovie_show
                  endcase
  "YMIN_SLIDER" : begin
                  iymin=event.value < 999
                  if(iymin ge iymax) then begin
                    iymax=iymin+1
                    WIDGET_CONTROL, ymaxslider, SET_VALUE=iymax
                  endif
                  if(go eq 0) then xmovie_show
                  endcase
  "YMAX_SLIDER" : begin
                  iymax=event.value > 1
                  if(iymin ge iymax) then begin
                    iymin=iymax-1
                    WIDGET_CONTROL, yminslider, SET_VALUE=iymin
                  endif
                  if(go eq 0) then xmovie_show
                  endcase
  ELSE          : MESSAGE, "Event User Value Not Found"
ENDCASE

return
end

pro xmovie_show,frame
;+
;   xmovie_show,frame
;
;            show individual frame
;
;-

COMMON cxdata,x1,y1,x2,y2,id_frame
COMMON cxmovie,go,speed,n_frames,frame_nr,frame_increment,$
 movie_base,bg_base,speedslider,frameslider,incrementslider,$
 xminslider,xmaxslider,ixmin,ixmax,xmin,xmax,$
 yminslider,ymaxslider,iymin,iymax,ymin,ymax

if(go eq 1) then return                             ; must be paused
if(n_elements(frame) ne 0) then frame_nr=frame

if(frame_nr ge n_frames) then frame_nr=frame_nr-n_frames
if(frame_nr lt 0) then frame_nr=frame_nr+n_frames
WIDGET_CONTROL, frameslider, SET_VALUE=frame_nr    ; display frame_nr

; plot frame

xrange=[xmin+ixmin/1000.*(xmax-xmin),xmin+ixmax/1000.*(xmax-xmin)]
yrange=[ymin+iymin/1000.*(ymax-ymin),ymin+iymax/1000.*(ymax-ymin)]

dum=size(x1)
if(dum(0) eq 1) then $
  plot,x1(*),y1(*,frame_nr),$
    xrange=xrange,yrange=yrange,xstyle=1,ystyle=1 $
else $
  plot,x1(*,frame_nr),y1(*,frame_nr),$
    xrange=xrange,yrange=yrange,xstyle=1,ystyle=1
if(n_elements(y2) gt 1) then begin
  dum=size(x2)
  if(dum(0) eq 1) then $
    oplot,x2(*),y2(*,frame_nr),line=2 $ 
  else $
    oplot,x2(*,frame_nr),y2(*,frame_nr),line=2
endif

xyouts,20,20,strtrim(string(id_frame(frame_nr)),2),/device

return
end

pro xmovie_start
;+
;   xmovie_start
;
;            start movie
;
;-

COMMON cxdata,x1,y1,x2,y2,id_frame
COMMON cxmovie,go,speed,n_frames,frame_nr,frame_increment,$
 movie_base,bg_base,speedslider,frameslider,incrementslider,$
 xminslider,xmaxslider,ixmin,ixmax,xmin,xmax,$
 yminslider,ymaxslider,iymin,iymax,ymin,ymax

WIDGET_CONTROL, speedslider, /SENSITIVE            ; enable speed slider
;WIDGET_CONTROL, frameslider, SENSITIVE=0           ; disable frame slider
WIDGET_CONTROL, xminslider, SENSITIVE=0            ; disable xrange slider
WIDGET_CONTROL, xmaxslider, SENSITIVE=0            ; disable xrange slider
WIDGET_CONTROL, yminslider, SENSITIVE=0            ; disable yrange slider
WIDGET_CONTROL, ymaxslider, SENSITIVE=0            ; disable yrange slider

go=1

WIDGET_CONTROL,bg_base,timer=0.0  ; start xmovie_bck

return
end

pro xmovie_bck, dummy
;+
;   xmovie_bck
;
;            background process for movie plotting 
;
;-

COMMON cxdata,x1,y1,x2,y2,id_frame
COMMON cxmovie,go,speed,n_frames,frame_nr,frame_increment,$
 movie_base,bg_base,speedslider,frameslider,incrementslider,$
 xminslider,xmaxslider,ixmin,ixmax,xmin,xmax,$
 yminslider,ymaxslider,iymin,iymax,ymin,ymax
COMMON cdirection,direction


if(go eq 0) then return                      ; paused ?
WIDGET_CONTROL, frameslider, SET_VALUE=frame_nr    ; display frame_nr
dum=size(x1)
if(dum(0) eq 1) then $                       ; erase previous graph
  oplot,x1(*),y1(*,frame_nr),color=0 $
else $
  oplot,x1(*,frame_nr),y1(*,frame_nr),color=0
if(n_elements(y2) gt 1) then begin
  dum2=size(x2)
  if(dum2(0) eq 1) then $
    oplot,x2(*),y2(*,frame_nr),line=2,color=0 $
  else $
    oplot,x2(*,frame_nr),y2(*,frame_nr),line=2,color=0
endif
xyouts,20,20,strtrim(string(id_frame(frame_nr)),2),color=0,/device
frame_nr = frame_nr + direction              ; increment frame_nr
if(frame_nr eq n_frames) then frame_nr=0     ; cyclic movie
if(frame_nr eq -1) then frame_nr=n_frames-1
if(dum(0) eq 1) then $                       ; erase previous graph
  oplot,x1(*),y1(*,frame_nr) $
else $
  oplot,x1(*,frame_nr),y1(*,frame_nr)
if(n_elements(y2) gt 1) then begin
  dum2=size(x2)
  if(dum2(0) eq 1) then $
    oplot,x2(*),y2(*,frame_nr),line=2 $
  else $
    oplot,x2(*,frame_nr),y2(*,frame_nr),line=2
endif
xyouts,20,20,strtrim(string(id_frame(frame_nr)),2),/device

WIDGET_CONTROL, bg_base, timer=1./speed

return
end

pro xmovie_stop

COMMON cxdata,x1,y1,x2,y2,id_frame
COMMON cxmovie,go,speed,n_frames,frame_nr,frame_increment,$
 movie_base,bg_base,speedslider,frameslider,incrementslider,$
 xminslider,xmaxslider,ixmin,ixmax,xmin,xmax,$
 yminslider,ymaxslider,iymin,iymax,ymin,ymax

IF(go eq 1) then begin                             ; pause xmovie_bck
  go=0                                               
endif
WIDGET_CONTROL, speedslider, SENSITIVE=0           ; disable speed slider
WIDGET_CONTROL, frameslider, SET_VALUE=frame_nr    ; display frame_nr
WIDGET_CONTROL, frameslider, /SENSITIVE            ; enable frame slider
WIDGET_CONTROL, xminslider, /SENSITIVE             ; enable xrange slider
WIDGET_CONTROL, xmaxslider, /SENSITIVE             ; enable xrange slider
WIDGET_CONTROL, yminslider, /SENSITIVE             ; enable yrange slider
WIDGET_CONTROL, ymaxslider, /SENSITIVE             ; enable yrange slider


return
end

pro xmovie_psym

COMMON cxmovie,go,speed,n_frames,frame_nr,frame_increment,$
 movie_base,bg_base,speedslider,frameslider,incrementslider,$
 xminslider,xmaxslider,ixmin,ixmax,xmin,xmax,$
 yminslider,ymaxslider,iymin,iymax,ymin,ymax

!psym=-1-!psym         ; toggle !psym
xmovie_show,frame_nr

return
end

pro xmovie_dx_fix

COMMON cfixdx,dx_fix

dx_fix=1000-dx_fix         ; toggle dx_fix

return
end

pro xmovie,xt1,yt1 ,xt2,yt2, id=id,norm=norm,reverse=reverse,psym=psym,$
 xrange=xrange,yrange=yrange,fixdx=fixdx,GROUP=GROUP
;+
;   xmovie,x1,y1 [,x2,y2,id=id,norm=norm,/reverse,/psym]
;   xmovie,y1 [,id=id,norm=norm,/reverse,/psym]
;
;            X-version of movie routine
;
;-
COMMON cxdata,x1,y1,x2,y2,id_frame
COMMON cxmovie,go,speed,n_frames,frame_nr,frame_increment,$
 movie_base,bg_base,speedslider,frameslider,incrementslider,$
 xminslider,xmaxslider,ixmin,ixmax,xmin,xmax,$
 yminslider,ymaxslider,iymin,iymax,ymin,ymax
COMMON cdirection,direction
COMMON cfixdx,dx_fix

if(n_params(0) eq 0) then begin
  print,'xmovie,x1,y1 ,x2,y2,id=id,/norm,/reverse,/fixdx,/psym'
  return
endif

y2=0
IF (n_params(0) EQ 1) THEN BEGIN   ; one parameter = y1
  y1=reform(xt1)                   ; remove dimensions of size 1
  dum=size(y1)
  maxt=dum(dum(0))                 ; time index is last  index
  mdep=dum(1)                      ; x-index    is first index
  x1=indgen(mdep)                  ; set x1=index
ENDIF ELSE BEGIN
  x1=reform(xt1)                   ; remove dimensions of size 1
  y1=reform(yt1)                   ; remove dimensions of size 1
  if(n_elements(xt2) ne 0) then begin
    x2=reform(xt2)
    y2=reform(yt2)
  endif
ENDELSE
if(n_elements(norm) eq 0) then norm=0
if(n_elements(fixdx) eq 0) then dx_fix=0 else dx_fix=1000

if(norm eq 1) then begin
  y1=y1-min(y1)
  y1=y1/max(y1)
  if(n_elements(y2) gt 1) then begin
    y2=y2-min(y2)
    y2=y2/max(y2)
  endif
endif

; set n_frames

dum=size(y1)
n_frames=dum(dum(0))
if(n_elements(y2) gt 1) then begin
  dum2=size(y2)
  n_frames=n_frames < dum2(dum2(0))
endif
n_x=dum(1)
if(n_elements(y2) gt 1) then begin
  n_x=n_x < dum2(1)
endif
if(n_elements(id) eq 0) then id_frame=indgen(n_frames) else id_frame=id

; check validity of inputs

if(n_elements(id_frame) ne n_frames) then begin
  print,'n_elements(id_frame) ne n_frames'
  return
endif
if(dum(0) lt 2) then begin
   n_frames=1
;  print,'y1 must have at least two dimensions'
;  return
endif
;if(n_elements(y2) gt 1) then begin
;  if(max(abs(size(y2)-size(y1))) ne 0) then begin
;    print,'y1 and y2 have different dimensions'
;    return
;  endif
;endif

; set defaults

speed=5
if(n_elements(reverse) ne 0) then begin
  direction=-1 
  frame_nr=n_frames-1
endif else begin
  direction=+1
  frame_nr=0
endelse
frame_increment=1
go=0
if(n_elements(xrange) eq 0) then begin
  if(!x.range(0) eq !x.range(1)) then begin
    if(n_elements(y2) le 1) then begin
      minx=min(x1)
      maxx=max(x1)
    endif else begin
      minx=min(x1) < min(x2)
      maxx=max(x1) > max(x2)
    endelse
    dx=maxx-minx
    xrange=[minx-0.05*dx,maxx+0.05*dx] 
  endif else xrange=!x.range
endif
xmin=xrange(0)
xmax=xrange(1)
ixmin=0
ixmax=1000
if(n_elements(yrange) eq 0) then begin
  if(!y.range(0) eq !y.range(1)) then begin
    if(n_elements(y2) le 1) then begin
      miny=min(y1)
      maxy=max(y1)
    endif else begin
      miny=min(y1) < min(y2)
      maxy=max(y1) > max(y2)
    endelse
    dy=maxy-miny
    yrange=[miny-0.01*dy,maxy+0.01*dy] 
  endif else yrange=!y.range
endif
ymin=yrange(0)
ymax=yrange(1)
iymin=0
iymax=1000
psymold=!psym
if(n_elements(psym) ne 0) then !psym=-1

movie_base = WIDGET_BASE(TITLE = "Xmovie")  ; create the main base
r_base = WIDGET_BASE(movie_base,/row)
c_base = WIDGET_BASE(r_base,/column)
XPdMenu, [ '" Done "		XMOVIE_DONE',$
           '" Start"		XMOVIE_START',$
           '" Stop "		XMOVIE_STOP',$
           '" Psym "		XMOVIE_PSYM',$
           '" Help "		XMOVIE_HELP'], c_base  ; top buttons

XPdMenu, [ '" Dxfix "		XMOVIE_DX_FIX'], c_base  ; next buttons

speedslider = WIDGET_SLIDER(c_base, $
      XSIZE = 256, $
      MINIMUM = 1, $
      MAXIMUM = 100, $
      VALUE = speed, $
      /DRAG, $
      TITLE = 'Animation speed', $
      UVALUE = "SPEED_SLIDER")

frameslider = WIDGET_SLIDER(c_base, $
      XSIZE = 256, $
      MINIMUM = 0, $
      MAXIMUM = (n_frames-1) > 1, $
      VALUE = frame_nr, $
      /DRAG, $
      TITLE = 'Frame number', $
      UVALUE = "FRAME_SLIDER")

XPdMenu, [ '"   -   "		XMOVIE_BKW',$
           '"   +   "		XMOVIE_FWD'], c_base

incrementslider = WIDGET_SLIDER(c_base, $
      XSIZE = 256, $
      MINIMUM = 1, $
      MAXIMUM = (n_frames/2) > 2, $
      VALUE = frame_increment, $
      TITLE = 'Frame increment', $
      UVALUE = "INCR_SLIDER")

draw = WIDGET_DRAW(r_base, $
      XSIZE=700, $
      YSIZE=550)

xminslider = WIDGET_SLIDER(c_base, $
      XSIZE = 256, $
      MINIMUM = 0, $
      MAXIMUM = 1000, $
      VALUE = 0, $
      /DRAG,$
      TITLE = 'x-min', $
      UVALUE = "XMIN_SLIDER")

xmaxslider = WIDGET_SLIDER(c_base, $
      XSIZE = 256, $
      MINIMUM = 0, $
      MAXIMUM = 1000, $
      VALUE = 1000, $
      /DRAG,$
      TITLE = 'x-max', $
      UVALUE = "XMAX_SLIDER")

yminslider = WIDGET_SLIDER(c_base, $
      XSIZE = 256, $
      MINIMUM = 0, $
      MAXIMUM = 1000, $
      VALUE = 0, $
      /DRAG,$
      TITLE = 'y-min', $
      UVALUE = "YMIN_SLIDER")

ymaxslider = WIDGET_SLIDER(c_base, $
      XSIZE = 256, $
      MINIMUM = 0, $
      MAXIMUM = 1000, $
      VALUE = 1000, $
      /DRAG,$
      TITLE = 'y-max', $
      UVALUE = "YMAX_SLIDER")

bg_base = WIDGET_BASE(event_pro="xmovie_bck") ; set up background widget
WIDGET_CONTROL, bg_base, /REALIZE

WIDGET_CONTROL, movie_base, /REALIZE          ; create the widgets

xmovie_show,frame_nr

XManager, "xmovie", movie_base, $
      EVENT_HANDLER = "xmovie_ev", $
      GROUP_LEADER = GROUP

!psym=psymold

return
end


