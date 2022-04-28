;+
; NAME:
;	DEFAULT
;
; PURPOSE:
;	This procedure sets default extension for input files
;	idl1
;	idlcnt
;	idlny
;	idlopc
;	dumc
;	jny
;
;	extension='none' sets original upper case file names.
;
; CATEGORY:
;	Multi
;
; CALLING SEQUENCE:
;
;	DEFAULT, Extension
;
; INPUTS:
;	Extension:	File extension, 'none' to set original upper case
;			file names
;
; OUTPUTS:
;	in common block:
;	def_ext		file extension
;       openfile=0	to signal reopening of files4
;
; COMMON BLOCKS:
;	common_multi
;
; MODIFICATION HISTORY:
; 	Written by:	Mats Carlsson
;       03-03-14        added common block cnyrd to force opening of
;                       new files in nyrd when changing atmosphere
;-
pro default,extension

@common_multi
common cnyrd,openfile,lu2,ljny

if(n_params(0) lt 1) then begin
  print,'default,extension'
  return
endif

def_ext=strtrim(string(extension),2)
openfile=0

return
end
