pro read_vald,filein,fileout ,verbose=verbose
;+
;   read_vald,filein,fileout ,/verbose
;
;            reads the VALD lines into an abslin-configuration
;-
if(n_params() lt 2) then begin
  print,'read_vald,filein,fileout ,/verbose'
  return
endif

c = 2.9979247e+8                ;velocity of light [m/s]
h = 6.6261758e-34               ;planck constant [Js]
ev = 1.602117e-19               ;electron volt [J]

elements=['H ','HE','C ','N ','O ','NE','NA','MG','AL','SI','S ','K ','CA','SC','TI','V ','CR','MN','FE','CO','NI']

;first find how many lines, i:
openr, lu, filein, /get_lun 
header = ''
while(strpos(header,'Lande') lt 0) do begin 
  readf, lu, header
  if(keyword_set(verbose)) then print, header
endwhile
readf, lu, header

i = 0
on_ioerror, jumpto1
while 1 do begin 
    line_1 = strarr(1)
    line_2 = strarr(1)
    readf, lu, line_1
    readf, lu, line_2
    element = strmid(line_1, 1, 2)
    element = strupcase(element)
    ion = strmid(line_1, 4, 1)
    ion = fix(ion)
    if ion gt 2 then begin
        if(keyword_set(verbose)) then print, 'ION > 2'
    endif else if(where(element[0] eq elements) eq -1) then begin
       if(keyword_set(verbose)) then print, 'Unknown element '+element[0]
    endif else begin
        i = i+1      
    endelse
endwhile

jumpto1:
free_lun, lu


;read the lines into the abslin-configuration file, 'multi_lines.dat':
openr, lu, filein, /get_lun
openw, lur, fileout, /get_lun
openw, lur_2, 'lines_lambda.dat', /get_lun
printf, lur, i
printf, lur_2, i
header = ''
while(strpos(header,'Lande') lt 0) do begin 
  readf, lu, header
endwhile
readf, lu, header

on_ioerror, jumpto
while 1 do begin 
    line1 = strarr(1)
    readf, lu, line1
    dum=str_sep(line1,',')
    if(keyword_set(verbose)) then print, line1
    
    element = strmid(dum[0], 1, 2)
    element = strupcase(element)
    if(keyword_set(verbose)) then print, element
    ion = strmid(dum[0], 4, 1)
    ion = fix(ion)
    if(keyword_set(verbose)) then print, ion
    wl = float(dum[1])          ;lambda in air (lambda > 200nm : air)
    wl = inv_convl(wl)          ;lambda in vac
    if(keyword_set(verbose)) then print, wl
    lg_gf = dum[2]
    lg_gf = float(lg_gf)
    if(keyword_set(verbose)) then print, lg_gf
    exc_lo = dum[3]
    exc_lo = double(exc_lo)
    if(keyword_set(verbose)) then print, exc_lo
    lambda_lo = ((h/ev)*c)/exc_lo
    lambda_lo = 1./(lambda_lo*100)
    j_lo = dum[4]
    j_lo = float(j_lo)
    g_lo = 2*j_lo +1          
    if(keyword_set(verbose)) then print, j_lo
    exc_up = dum[5]
    exc_up = double(exc_up)
    if(keyword_set(verbose)) then print, exc_up
    lambda_up = ((h/ev)*c)/exc_up
    lambda_up = 1./(lambda_up*100)
    j_up = dum[6]
    j_up = float(j_up)
    g_up = 2*j_up + 1
    if(keyword_set(verbose)) then print, j_up
    lande_lo = dum[7]
    if(keyword_set(verbose)) then print, lande_lo
    lande_up = dum[8]
    if(keyword_set(verbose)) then print, lande_up
    lande_mean = dum[9]
    if(keyword_set(verbose)) then print, lande_mean
    rad = dum[10]
    rad = float(rad)
    stark = dum[11]
    stark = float(stark)
    waals = dum[12]
    waals = float(waals)
    el1 = strmid(dum[0], 0, 3)
    id = ''
    id = element[0]
if(keyword_set(verbose)) then print, rad, stark, waals
    if waals eq 0.00 then begin
        vdw = 0.00
    endif else begin
;; EXC_LO AND EXC_UP IN EV!!!
      if(ion le 2) then begin
        vdw = damp_vdw(waals, ion, exc_lo, exc_up, id)  
      endif else begin
        vdw = 1.0
      endelse
    endelse
        
    ;vdw = damp_vdw(waals, ion, exc_lo, exc_up, id)
    if(keyword_set(verbose)) then print, vdw
    if(keyword_set(verbose)) then print, rad
    if(keyword_set(verbose)) then print, stark
    if(keyword_set(verbose)) then print, waals
    
    if stark eq 0.00 then begin
        stark = 0.00
    endif else begin
        stark = 10^(stark)
    endelse
    
    if rad eq 0.00 then begin
        rad = 0.00
    endif else begin
        rad = 10^(rad)
    endelse
   
    if(keyword_set(verbose)) then print, rad
    if(keyword_set(verbose)) then print, stark
    if(keyword_set(verbose)) then print, waals
    
    line2 = strarr(1)
    readf, lu, line2
    if(keyword_set(verbose)) then print, line2
    
;    ref = strmid(line2, 43, 1)
;    ref = fix(ref)
;    if(keyword_set(verbose)) then print, ref
    
    if ion gt 2 then begin
        if(keyword_set(verbose)) then print, 'ION > 2'
    endif else if(where(element[0] eq elements) eq -1) then begin
       if(keyword_set(verbose)) then print, 'Unknown element '+element[0]
    endif else begin
        printf, lur, ''       
        if ion eq 1 then begin
            ion_str = 'I'
            printf, lur, element, ion_str, format = '(a2, 1x, a1)'
            printf, lur, ion, lambda_lo, g_lo, format = '(3x, I1, F18.3, 1x, F4.1)'
            printf, lur, element, ion_str, format = '(a2, 1x, a1)'
            printf, lur, ion, lambda_up, g_up, format = '(3x, I1, F18.3, 1x, F4.1)'
            printf, lur_2, element, ion_str, wl, format = '(a2, 2x, a1, 2x, F8.3)'
        endif
        if ion eq 2 then begin
            ion_str = 'II'
            printf, lur, element, ion_str, format = '(a2, 1x, a2)'
            printf, lur, ion, lambda_lo, g_lo, format = '(3x, I1, F18.3, 1x, F4.1)'
            printf, lur, element, ion_str, format = '(a2, 1x, a2)' 
            printf, lur, ion, lambda_up, g_up, format = '(3x, I1, F18.3, 1x, F4.1)'
            printf, lur_2, element, ion_str, wl, format = '(a2, 1x, a2, 2x, F8.3)'
        endif
        printf, lur, wl, wl-2, wl+2, format = '(3(1x, F9.3))'
;        printf, lur, (10^(lg_gf))/g_lo, rad, '1.00', stark, format = '(3x, E10.4, 1x, E9.3, 2x, F4.2, 2x, E9.3)'
        if(vdW gt 0.0) and (vdW lt 100.) then vdW_out=vdW else vdW_out=1.0
        printf, lur, (10^(lg_gf))/g_lo, rad, vdW_out, stark, format = '(3x, E10.4, 1x, E9.3, 2x, F4.2, 2x, E9.3)'
                                ; endelse
    endelse
endwhile
jumpto:
if(keyword_set(verbose)) then print, i
free_lun, lu
free_lun, lur
free_lun, lur_2
end
