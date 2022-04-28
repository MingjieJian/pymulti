pro iondata,id,awgt,abnd,xi,u0,u1 ,verbose=verbose
;+
;   iondata,id,awgt,abnd,xi,u0,u1 ,/verbose
;
;            returns ionic data for element id
;            awgt    atomic weight in amu
;            abnd    solar abundance according to Asplund, Grevesse, Sauval
;            xi     ionization potential for ground state (eV)
;            u0      partition function for ground state at 5040 K
;            u1      partition function for ground state at 5040 K
;-
if(n_params() lt 6) then begin
  print,'iondata,id,awgt,abnd,xi,u0,u1 ,/verbose'
  return
endif

; element id array
;  0    1    2    3    4    5    6    7    8    9
id_el=['H ','He','Li','Be','B ','C ','N ','O ','F ',$  ; 0x
  'Ne','Na','Mg','Al','Si','P ','S ','Cl','Ar','K ',$  ; 1x
  'Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu',$  ; 2x
  'Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ',$  ; 3x
  'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In',$  ; 4x
  'Sn','Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr',$  ; 5x
  'Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm',$  ; 6x
  'Yb','Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au',$  ; 7x
  'Hg','Tl','Pb','Bi','Po','At','Rn','Fa','Ra','Ac',$  ; 8x
  'Th','Pa','U ']

; abundance
;   0     1     2     3     4     5     6     7     8     9
abnd_el=[ 12., 10.93, 1.05, 1.38, 2.70, 8.39, 7.78, 8.66, 4.56,$ ;0x
    7.84, 6.17, 7.53, 6.37, 7.51, 5.36, 7.14, 5.50, 6.18, 5.08,$ ;1x
    6.31, 3.05, 4.90, 4.00, 5.64, 5.39, 7.45, 4.92, 6.23, 4.21,$ ;2x
    4.60, 2.88, 3.58, 2.29, 3.33, 2.56, 3.28, 2.60, 2.92, 2.21,$ ;3x
    2.59, 1.42, 1.92,-9.99, 1.84, 1.12, 1.69, 0.94, 1.77, 1.60,$ ;4x
    2.00, 1.00, 2.19, 1.51, 2.27, 1.07, 2.17, 1.13, 1.58, 0.71,$ ;5x
    1.45,-9.99, 1.01, 0.52, 1.12, 0.28, 1.14, 0.51, 0.93, 0.00,$ ;6x
    1.08, 0.06, 0.88,-0.17, 1.11, 0.23, 1.45, 1.38, 1.64, 1.01,$ ;7x
    1.13, 0.90, 2.00, 0.65,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99,$ ;8x
    0.06,-9.99,-0.52]                                            ;9x
iw=where(abnd_el ne -9.99)
abnd_el[iw]=10.^(abnd_el[iw]-12.)
iw=where(abnd_el eq -9.99)
abnd_el[iw]=0.00

; atomic weights
;   0     1     2     3     4     5     6     7     8     9
awgt_el=[  1.0,  4.0, 6.94, 9.01, 10.8, 12.0, 14.0, 16.0, 19.0,$ ; 0x
    20.2, 23.0, 24.3, 27.0, 28.1, 31.0, 32.1, 35.5, 39.9, 39.1,$ ; 1x
    40.1, 45.0, 47.9, 50.9, 52.0, 54.9, 55.8, 58.9, 58.7, 63.5,$ ; 2x
    65.4, 69.7, 72.6, 74.9, 79.0, 79.9, 83.8, 85.5, 87.6, 88.9,$ ; 3x
    91.2, 92.9, 95.9, 99.0,101.1,102.9,106.4,107.9,112.4,114.8,$ ; 4x
   118.7,121.8,127.6,126.9,131.3,132.9,137.3,138.9,140.1,140.9,$ ; 5x
   144.2,145.0,150.4,152.0,157.3,158.9,162.5,164.9,167.3,168.9,$ ; 6x
   173.0,175.0,178.5,180.9,183.9,186.2,190.2,192.2,195.1,197.0,$ ; 7x
   200.6,204.4,207.2,209.0,210.0,210.0,222.0,223.0,226.0,227.0,$ ; 8x
   232.0,231.0,238.0]                                            ; 9x

; ionization potential in eV for neutral stage
;   0     1     2     3     4     5     6     7     8     9
xi_el= [13.60,24.58, 5.39, 9.32, 8.30,11.26,14.53,13.61,17.42,$ ; 0x
   21.56, 5.14, 7.64, 5.98, 8.15,10.48,10.36,13.01,15.75, 4.34,$ ; 1x
    6.11, 6.54, 6.82, 6.74, 6.76, 7.43, 7.87, 7.86, 7.63, 7.72,$ ; 2x
    9.39, 6.00, 7.88, 9.81, 9.75,11.84,14.00, 4.18, 5.69, 6.38,$ ; 3x
    6.84, 6.88, 7.10, 7.28, 7.36, 7.46, 8.33, 7.57, 8.99, 5.78,$ ; 4x
    7.34, 8.64, 9.01,10.45,12.13, 3.89, 5.21, 5.61, 6.90, 5.80,$ ; 5x
    6.30,99.99, 5.60, 5.67, 6.16, 6.70, 6.80, 6.00, 6.00, 6.00,$ ; 6x
    6.20, 6.10, 7.00, 7.88, 7.98, 7.87, 8.70, 9.00, 9.00, 9.22,$ ; 7x
   10.43, 6.11, 7.42, 7.29, 8.43, 9.30,10.75, 4.00, 5.28, 6.90,$ ; 8x
    6.00,99.99, 6.00]                                            ; 9x

; partition function for neutral element
; taken from Allen: Astrophysical Quantities
;   0     1     2     3     4     5     6     7     8     9
u0_el=[   0.30, 0.00, 0.32, 0.01, 0.78, 0.97, 0.61, 0.94, 0.75,$ ; 0x
    0.00, 0.31, 0.01, 0.77, 0.98, 0.65, 0.91, 0.72, 0.00, 0.34,$ ; 1x
    0.07, 1.08, 1.48, 1.62, 1.02, 0.81, 1.43, 1.52, 1.47, 0.36,$ ; 2x
    0.00, 0.73, 0.91,-9.99, 0.83,-9.99, 0.00, 0.36, 0.10, 1.08,$ ; 3x
    1.53,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99, 0.00,-9.99,$ ; 4x
    0.73,-9.99,-9.99,-9.99,-9.99,-9.99, 0.36, 1.41,-9.99,-9.99,$ ; 5x
   -9.99,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99,$ ; 6x
    0.02,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99,$ ; 7x
   -9.99,-9.99, 0.26,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99,$ ; 8x
   -9.99,-9.99,-9.99]                                            ; 9x
iw=where(u0_el ne -9.99)
u0_el[iw]=10.^u0_el[iw]
iw=where(u0_el eq -9.99)
u0_el[iw]=0.00

; partition function for singly ionized element
; taken from Allen: Astrophysical Quantities
;   0     1     2     3     4     5     6     7     8     9
u1_el=[   0.00, 0.30, 0.00, 0.30, 0.00, 0.78, 0.95, 0.60, 0.92,$ ; 0x
    0.73, 0.00, 0.31, 0.00, 0.76, 0.91, 0.62, 0.89, 0.69, 0.00,$ ; 1x
    0.34, 1.36, 1.75, 1.64, 0.86, 0.89, 1.63, 1.46, 1.02, 0.01,$ ; 2x
    0.30, 0.00, 0.64,-9.99, 0.60,-9.99, 0.62, 0.00, 0.34, 1.18,$ ; 3x
    1.66,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99, 0.30,-9.99,$ ; 4x
    0.52,-9.99,-9.99,-9.99,-9.99,-9.99, 0.62, 1.47,-9.99,-9.99,$ ; 5x
   -9.99,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99,$ ; 6x
    0.30,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99,$ ; 7x
   -9.99,-9.99, 0.32,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99,$ ; 8x
   -9.99,-9.99,-9.99]                                            ; 9x
iw=where(u1_el ne -9.99)
u1_el[iw]=10.^u1_el[iw]
iw=where(u1_el eq -9.99)
u1_el[iw]=0.00


; find element index

iel=where(strtrim(strupcase(id),2) eq strtrim(strupcase(id_el),2),count)
if(count ne 1) then begin
  if(keyword_set(verbose)) then print,id,' not found among ',id_el
  abnd=0.0
  return
endif

awgt=awgt_el[iel]
abnd=abnd_el[iel]
xi=xi_el[iel]
u0=u0_el[iel]
u1=u1_el[iel]
xi=xi[0]
u0=u0[0]
u1=u1[0]
return
end