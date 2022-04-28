PRO SPEC_ct,RED,GREEN,BLUE
;+
;   spec_ct,red,green,blue
;
;            produce and load a spectrum colour lookup table
;
;-
RED=INTARR(256)
GREEN=INTARR(256)
BLUE=INTARR(256)
INT100=FINDGEN(100)*255./99.
INT60=FINDGEN(60)*255./59.
INT30=FINDGEN(30)*255./29.
RED(1:60)=255.
GREEN(1)=INT60
;RED(61:71)=255.
RED(61)=255.-INT60
GREEN(61:150)=255.
GREEN(151)=255.-INT60
BLUE(101)=INT60
BLUE(161:250)=255.
RED(191)=INT60*0.7
RED(251:254)=255*0.7
GREEN(251:254)=0
BLUE(251:254)=255
RED(255)=255
GREEN(255)=255
BLUE(255)=255

n_colors=!d.n_colors < 256
red=interpol(red,n_colors)
green=interpol(green,n_colors)
blue=interpol(blue,n_colors)

TVLCT,RED,GREEN,BLUE

RETURN
END




