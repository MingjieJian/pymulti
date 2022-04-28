pro rout3d,file ,verbose=verbose,swap_endian=swap_endian
;+
;   rout3d,file ,/verbose
;
;            reads OUT3D file from mul23guc_3d
;-
@common_multi3d
common cout3d,id3d,taulg3d,cmasslg3d,xnorm3d,icont,x3d,height3d
common civ3d,iv3d
common cn3d,n3d

if(n_params() lt 1) then begin
  print,'rout3d,file ,/verbose'
  return
endif

; find endianness of file
openr,lur,file,/get_lun
rec=assoc(lur,lonarr(1))
test=rec[0]
free_lun,lur
if(test eq 16) then swap_endian=0 else swap_endian=1

isize=0L
itype=0L
cname='12345678'
openr,lur,file,/get_lun,/f77,swap_endian=swap_endian
while(not eof(lur)) do begin
  readu,lur,itype,isize,cname
  if(keyword_set(verbose)) then print,itype,isize,' ',cname
  case(strtrim(cname,2)) of
    'id': begin
          id3d=string(' ',format='(a80)')
          readu,lur,id3d
          id3d=strtrim(id3d,2)
	end
    'dim': begin
	  nx=0L & ny=0L & ndep=0L & mq=0L & nrad=0L
	  if(isize eq 5) then version=1 else version=2
	  if(version eq 1) then begin
	    readu,lur,nx,ny,ndep,mq,nrad
	  endif else if(version eq 2) then begin
	    nq=lonarr(isize-5)
	    readu,lur,nx,ny,ndep,mq,nrad,nq
            nqtot=total(nq)+nrad
          endif
        end
    'q': begin
          q=fltarr(mq,nrad)
          readu,lur,q
        end
    'xl': begin
          xl=dblarr(nqtot)
          readu,lur,xl
        end
    'taulg3d': begin
          taulg3d=fltarr(nx,ny,ndep)
          readu,lur,taulg3d
        end
    'cmass3d': begin
          cmasslg3d=fltarr(nx,ny,ndep)
          readu,lur,cmasslg3d
        end
    'dscal2': begin
          dpnew3d=fltarr(nx,ny,ndep)
          readu,lur,dpnew3d
        end
    'icont': begin
	  if(version eq 1) then begin
	    icont=fltarr(mq,nx,ny,nrad)
	  endif else begin
	    icont=fltarr(total(nq)+nrad,nx,ny)
	  endelse
	  readu,lur,icont
        end
    'Iv': begin
          iv3d=fltarr(total(nq)+nrad,nx,ny)
	  readu,lur,iv3d
        end
    'xnorm3d': begin
          xnorm3d=fltarr(nx,ny,ndep)
          readu,lur,xnorm3d
        end
    'x3d': begin
	  x3d=fltarr(nx,ny,ndep)
	  readu,lur,x3d
        end
    'height3d': begin
	  height3d=fltarr(nx,ny,ndep)
	  readu,lur,height3d
        end
    'n3d': begin
          nk=isize/(nx*ny*ndep)
	  n3d=fltarr(nx,ny,ndep,nk)
	  readu,lur,n3d
        end
    'nk': begin
          nk=0L
	  readu,lur,nk
	  n3d=fltarr(nx,ny,ndep,nk)
          dum=fltarr(nx,ny,ndep)
          ilev=0
        end
    else: begin
      if(strmid(cname,0,3) eq 'n3d') then begin
        readu,lur,dum
        n3d[*,*,*,ilev]=dum
        ilev=ilev+1
      endif else begin
          message,'case not found:'+cname,/info
          readu,lur,cname
      endelse
    end
  endcase
endwhile
free_lun,lur
end

@common_multi3d
common cout3d,id3d,taulg3d,cmasslg3d,xnorm3d,icont,x3d,height3d
common civ3d,iv3d
common cn3d,n3d
end
