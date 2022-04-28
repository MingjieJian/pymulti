;+
; NAME:
;	NYCLOSE
;
; PURPOSE:
;	This procedure closes the idlny and jny files.
;
; CATEGORY:
;	Multi
;
; CALLING SEQUENCE:
;
;	NYCLOSE
;
; INPUTS:
;	None
;
; OUTPUTS:
;	In common:
;	openfile=0	signals closed files
;	lu2=0		signals closed files
;	ljny=0		signals closed files
;
; COMMON BLOCKS:
;	common_multi
;
; MODIFICATION HISTORY:
; 	Written by:	Mats Carlsson
;-
pro nyclose,dummy

@common_multi

if(n_elements(lu2) ne 0) then $
 if(lu2 ne 0) then free_lun,lu2
if(n_elements(ljny) ne 0) then $
 if(ljny ne 0) then free_lun,ljny
openfile=0
lu2=0
ljny=0

return
end
