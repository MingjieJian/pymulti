      PROGRAM MULTI
C
C  MAIN CALLING PROGRAM FOR SOLVING MULTI-LEVEL NON-LTE PROBLEMS
C  FOR DOCUMENTATION, COMMENTS ON EXPERIENCE WITH THE CODE AND
C  A DESCRIPTION OF THE METHOD SEE
C  CARLSSON 1986: 'UPPSALA ASTRONOMICAL OBSERVATORY, REPORT NO. 33'.
C  A MORE DETAILED DESCRIPTION OF THE METHOD IS IN
C  SCHARMER, CARLSSON 1985: 'JOURNAL OF COMPUTATIONAL PHYSICS',59,56
C  FOR LIST OF VARIABLES AND FILES, SEE  ROUTINE  START.
C
C  PARAMETER VARIABLES:
C  MDEP1  NUMBER OF DEPTH POINTS. MUST BE THE NUMBER ACTUALLY USED.
C  MK1    NUMBER OF LEVELS INCLUDING CONTINUUM LEVELS. -''-.
C  MDEP   MAXIMUM NUMBER OF DEPTH POINTS.
C  MK     MAXIMUM NUMBER OF LEVELS INCLUDING CONTINUUM LEVELS
C  MLINE  MAXIMUM NUMBER OF BOUND-BOUND LINES
C  MWIDE  MAXIMUM NUMBER OF BROAD LINES + CONTINUA IN DETAIL
C  MRAD   MLINE+MWIDE
C  MRFIX  MAXIMUM NUMBER OF FIXED TRANSITIONS
C  MQ     MAXIMUM NUMBER OF FREQUENCY POINTS IN ONE TRANSITION
C  MMU    MAXIMUM NUMBER OF ANGLE POINTS
C
C  ONLY MDEP1 AND MK1 MUST BE EQUAL TO THE VALUES ACTUALLY USED, NDEP1
C  AND NK1. THE OTHERS MAY BE LARGER THAN VALUES ACTUALLY USED.
C
C: MULTI  90-07-31  MODIFICATIONS: (MATS CARLSSON)
C:        CHANGED TO CONTAIN DECLARATION OF W AND E
C:        THESE PASSED IN ARGUMENT LIST INSTEAD OF IN COMMON TO
C:        AVOID NK=MK1 AND NDEP=MDEP1 REQUIREMENT
C:
C:        09-03-04  MODIFICATIONS: (MATS CARLSSON)
C:        THIS IS A WRAPPER TO ALLOW CALCULATIONS COLUMN BY COLUMN
C:        IN A 3D ATMOSPHERE
C:
      use var_3d
      INCLUDE 'PREC'
      INCLUDE 'PARAM'
      INCLUDE 'PARAMW'
      INCLUDE 'CATOM'
      INCLUDE 'CATMOS'
      INCLUDE 'CATMO2'
      INCLUDE 'CTRAN'
      INCLUDE 'CSLINE'
      INCLUDE 'CGAUSI'
      INCLUDE 'CCONST'
      INCLUDE 'CINPUT'
      INCLUDE 'CLGMX'
      INCLUDE 'CLU'
      INCLUDE 'COPCL'
      INCLUDE 'C3D'
      integer :: irun,i,ix,iy,ixx,iyy,k0,k1,k2,nu1,kr,icrsw0,iconv0
      integer :: lent,iredo_switch,id,nk1,ndep1,ldum,lopnmx,count
      integer :: count_conv,count_div,count_convl,count_divl
      DIMENSION W(MK1,MDEP1,MK1,MDEP1),E(MK1,MDEP1)
      dimension dpnew(mdep)
c
      debug=.false.
c
c  init MPI. set mpnproc and mpmyid for serial runs before initializing
c
      mpnproc=1
      mpmyid=0
      call mpi_init(mperr)
      call mpi_comm_size(MPI_COMM_WORLD,mpnproc,mperr)
      call mpi_comm_rank(MPI_COMM_WORLD,mpmyid,mperr)
      allocate(icxy(0:mpnproc),stat=ierr)
      if(ierr.ne.0) write(*,*) 'allocate of icxy failed'
      master=mpmyid.eq.masterid
      slave=.not.master
      mpmyid2=mpmyid
c
c  read INPUT_3D file
c
      call rinput_3d
c
c  calculate nx, ny and icdo
c
      call initiate_3d
c
      do 995 irun=1,nruns     ! loop over atmos3d files
c
c  open global files always needed for multi compatability
c
        CALL MOPEN(LOUT,'OUT',1,'MYNEW')
        CALL MOPEN(LJOBLO,'JOBLOG',1,'MYNEW')
        CALL MOPEN(LDUMS,'DUMS',1,'MYNEW')
        do lphi=1,maxlu
          if(.not.lu(lphi)) goto 5
        enddo
        call stop('main: no free unit number for lphi')
    5   continue
        lu(lphi)=.true.
C
C  read input parameters
C
        run_ext=bext(irun)
        CALL RINPUT
c
c  set printout parameters to zero to avoid excessive output
c
        call disable_printout
c        
        icrsw0=icrsw         ! save icrsw for reset
        iconv0=iconv         ! save iconv for reset
        istart2=istart       ! for transfer in var_3d
        mq2=mq               ! for transfer in var_3d
c
c  read out3d0 file to set icdo array if icont_opt.eq.-1
c
        count_convl=0
        count_divl=0
        if(icont_opt.eq.-1) then
          call rout3d0(count)
          if(count.eq.0) goto 990   ! no redo, next run
        endif
c
c  distribute columns over processors
c
        call distrib_columns
c
c  read atmos3d file and save variables in processor local variables
c
        call ratmos3d
c
c  go through all columns, calculate and save output to local variables
c
        do i=1,icxy(mpmyid)
          ix=ixc(i,mpmyid)
          iy=iyc(i,mpmyid)
          ixx=ixxc(i,mpmyid)
          iyy=iyyc(i,mpmyid)
c
c  interpolate to new depthscale according to dscopt option
c
          write(*,'(a,5i5)') 'id,ix,iy,ixx,iyy ',
     &     mpmyid,ix,iy,ixx,iyy
          write(ljoblo,'(a,5i5)') 'id,ix,iy,ixx,iyy ',
     &     mpmyid,ix,iy,ixx,iyy
          if(dscopt(1:lent(dscopt)).eq.'dtrho') then
            call ipol_dscale(ndep2,height_in,temp3d(:,i),
     &       ne3d(:,i),vz3d(:,i),rho3d(:,i))
          else if(dscopt(1:6).eq.'dscal2') then
            call ipol_dscal2(ndep2,height_in,temp3d(:,i),
     &       ne3d(:,i),vz3d(:,i),rho3d(:,i),
     &       dpnew3d(:,i))
          else
            temp=temp3d(:,i)
            ne=ne3d(:,i)
            vel=vz3d(:,i)
            rho=rho3d(:,i)
            height=height_in
          endif
          if(istart2.eq.-1) then
            call ipol_n(ndep2,nk2,height_in,n3d(:,:,i))
          else if(istart2.eq.-2) then
            do k=1,ndep2
              do j=1,nk2
                n(j,k)=n3d(k,j,i)
              enddo
            enddo
          endif
          ndep=ndep2
c
c  convert height to km and make height increase upwards
c
          if(height(2).gt.height(1)) then
            height=-height(1:ndep2)*1.d-5
          else
            height=height(1:ndep2)*1.d-5
          endif
          call rewind(lopc)
          call rewind(lphi)
c
          iredo_switch=1
          iconv=iconv0
          CALL START
          IF(NK.GT.MK1) CALL STOP('MULTI: NK.GT.MK1')
          IF(NDEP.GT.MDEP1) CALL STOP('MULTI: NDEP.GT.MDEP1')
          NK1=NK
          NDEP1=NDEP
          nk2=nk           ! for transfer through var_3d
          ndep2=ndep
 97       continue
          CALL ITER(W,E,NK1,NDEP1)
c
c  check for NaN
c
          do k=1,ndep
            do j=1,nk
              if(n(j,k).ne.n(j,k)) iconv=0
            enddo
          enddo
c
c  if not converged and cr switching:
c    increase switching steps per decade by factor of two
c    this is tried iredo_max times
c  if converged, reset icrsw to INPUT values
c
          if(iconv.ne.1 .and. icrsw.gt.0 .and. 
     &     iredo_switch.le.iredo_max) then
            icrsw=icrsw*2
            iredo_switch=iredo_switch+1
            write(ljoblo,'(a,3i5)') 'retrying: id,ix,iy',
     &       mpmyid,ix,iy
            call initia
            goto 97
          endif
          icrsw=icrsw0
          if(iconv.ne.1 .and. itmax.gt.0) then    ! no convergence
            tau(:)=1.d-10                         ! fill output variables
            cmass(:)=1.d-10                       ! with values to avoid
            outint(:,:,:)=0.d0                    ! printout of NaN
            xnorm(:)=0.d0
            xsave(:)=0.d0
            dpnew(:)=height(:)
            n(:,:)=0.d0
            write(*,'(a,3i5)') 'diverged: id,ix,iy',
     &       mpmyid,ix,iy
            write(ljoblo,'(a,3i5)') 'diverged: id,ix,iy',
     &       mpmyid,ix,iy
            count_divl=count_divl+1
          else
            write(*,'(a,3i5)') 'converged: id,ix,iy',mpmyid,ix,iy
            write(ljoblo,'(a,3i5)') 'converged: id,ix,iy',mpmyid,ix,iy
            CALL FORMAL
            call dscal2(dpnew)
            count_convl=count_convl+1
          endif
c 
c  calculate total number of frequencies for packing of Iv
c  done here because need to be after call to start
c  only done once per processor per atmosphere
c  if(master) fill xl with wavelengths in AA in vacuum
c  allocate ouput local variables if first column
c
          if(i.eq.1) then
            if(nrad_write.eq.-1) nrad_write=nrad
            nqtot=0
            do kr=1,nrad
              if(lkr_write(kr)) nqtot=nqtot+nq(kr)+1
            enddo
            call calc_xl
            call allocate_out3dvar
          endif
c
c  save results in local output variables
c
          taulg3d(:,i)=log10(tau(1:ndep))
          cmasslg3d(:,i)=log10(cmass(1:ndep))
          nu1=1
          do kr=1,nrad
            if(lkr_write(kr)) then
              Iv(nu1:nu1+nq(kr),i)=outint(0:nq(kr),nmu,kr)
              nu1=nu1+nq(kr)+1
            endif
          enddo
          xnorm3d(:,i)=xnorm(1:ndep)
          x3d(:,i)=xsave(1:ndep)*xnorm(1:ndep)
          height3d(:,i)=height(1:ndep)
          dpnew3d(:,i)=dpnew(1:ndep)
          do j=1,nk
            n3d(:,j,i)=n(j,1:ndep)
            b3d(:,j,i)=n(j,1:ndep)/nstar(j,1:ndep)
          enddo
C
          write(ljoblo,'(a)') 'NORMAL END'
          CALL REWIND(LOPC)
          CALL REWIND(LPHI)
        enddo   ! end of column loop
c
c  calculations ready, synchronize
c
        call mpi_barrier(MPI_COMM_WORLD,mperr)
c
c  write output to files if any column converged
c
        if(mpnproc.gt.1) then 
          call mpi_allreduce(count_convl,count_conv,1,MPI_INTEGER,
     &     MPI_SUM,MPI_COMM_WORLD,mperr)
          call mpi_allreduce(count_divl,count_div,1,MPI_INTEGER,
     &     MPI_SUM,MPI_COMM_WORLD,mperr)
        else
          count_conv=count_convl
          count_div=count_divl
        endif
        if(master) write(*,*) ' main: '//run_ext(1:lent(run_ext))//
     &   ' columns conv:',count_conv,' div:',count_div
        if(count_conv.gt.0) then
          call wout3d
          call wtaulg3d
          call wcmass3d
          call wxnorm3d
          call wx3d
          call wheight3d
          call wIv3d
          call wdscal2
          call wn3d
          call wb3d
        endif
c
c  close files and delete scratch files
c
  990   continue
        lopnmx=1
        call mopen(ldumi,'DUMI',1,'MYOLD')
        do i=1,maxlu
          if(lu(i)) then
            lopnmx=lopnmx+1
            if(i.eq.ldums .or. i.eq.lout .OR.
     *       i.eq.latmos .or. i.eq.ljoblo .OR.
     *       i.eq.ldumi) then
              close(i,status='delete')
              lu(i)=.false.
            else
              close(I)
              lu(i)=.false.
            endif
          endif
        enddo
        if(icont_opt.eq.-1 ) then
          if(count.gt.0) call deallocate_3dvar
        else
          call deallocate_3dvar
        endif
        call mpi_barrier(MPI_COMM_WORLD,mperr)
  995 continue   ! end of irun loop
c
c  deallocate variables spanning runs
c  variables allocated here and in rinput_3d and initiate_3d
c
      deallocate(icxy,stat=ierr)
      if(ierr.ne.0) write(*,*) 'deallocation of icxy failed'
      deallocate(kr_write,nq_write,lkr_write,stat=ierr)
      if(ierr.ne.0) write(*,*) 'deallocation of kr_write... failed'
      deallocate(idxy,icdo,stat=ierr)
      if(ierr.ne.0) write(*,*) 'deallocation of idxy,icdo failed'
      deallocate(bext)
      if(debug) then
        write(*,'(a,i4,a,i5)') 'id=',mpmyid,' max open files:',lopnmx
      endif
      call mpi_finalize(mperr)
      stop
      end
c
c************************************************************************
c
      subroutine rinput_3d
c
c  read parameters from INPUT_3D file
c  if krsave is negative, the number of transitions to write is
c  read and their indices. logical array lkr_write keeps track
c
      use var_3d
      include 'PREC'
      include 'PARAM'
      include 'C3D'
      integer :: linput_3d,kr,i
c
      allocate(kr_write(mrad),nq_write(mrad),lkr_write(mrad),stat=ierr)
      if(ierr.ne.0) write(*,*) 'allocation of kr_write... failed'
      call mopen(linput_3d,'INPUT_3D',1,'OLD')
      read(linput_3d,*) icont_opt,iredo_max
      read(linput_3d,*) dscopt
      read(linput_3d,*) ix0,ix1,ixstep
      read(linput_3d,*) iy0,iy1,iystep
      read(linput_3d,*) krsave,nusave
      if(krsave.lt.0) then
        krsave=abs(krsave)
        read(linput_3d,*) nrad_write
        read(linput_3d,*) (kr_write(i),i=1,nrad_write)
        do kr=1,mrad
          lkr_write(kr)=.false.
        enddo
        do i=1,nrad_write
          lkr_write(kr_write(i))=.true.
        enddo
      else
        nrad_write=-1  ! nrad not known yet, set nrad_write before output
        do kr=1,mrad
          lkr_write(kr)=.true.
        enddo
      endif
      read(linput_3d,*) nruns
      allocate(bext(nruns),stat=ierr)
      if(ierr.ne.0) write(*,*) 'allocation of bext failed'
      do i=1,nruns
        read(linput_3d,*) bext(i)
      enddo
c
      call close(linput_3d)
c
      end subroutine rinput_3d
c
c**********************************************************************
c
      subroutine initiate_3d
c
      use var_3d
      implicit none
      integer :: i,ix,iy,ierr,lout3d,itype,isize,lent
      character :: cdummy
      character(len=8) :: cname
      real(4) :: dummy
c
c  initializations that can be done before reading INPUT and ATMOS
c
      nx=0
      ny=0
      do ix=ix0,ix1,ixstep
        nx=nx+1
      enddo
      do iy=iy0,iy1,iystep
        ny=ny+1
      enddo
      allocate(idxy(nx,ny),stat=ierr)
      if(ierr.ne.0) write(*,*) 'allocation of idxy failed'
      allocate(icdo(nx,ny),stat=ierr)
      if(ierr.ne.0) write(*,*) 'allocation of icdo failed'
      icdo(:,:)=1  ! default is to do all columns

      end subroutine initiate_3d
c
c**********************************************************************
c
      subroutine rout3d0(count)
c
      use var_3d
      implicit none
      integer :: i,ix,iy,ierr,lout3d,itype,isize,lent,count
      character :: cdummy
      character(len=8) :: cname
      real(4) :: dummy
c
c  if icont_opt=-1 out3d and Iv3d files are read and only columns where 
c  Iv.eq.0 are calculated
c  master reads the file sets the
c  array icdo to 1 (redo column) or 0
c
      if(master) then
        call mopen(lout3d,'out3d.'//run_ext(1:lent(run_ext)),0,
     &   'old')
        read(lout3d) itype,isize,cname    ! id
        write(*,*) itype,isize,cname
        read(lout3d) cdummy
        read(lout3d) itype,isize,cname    ! dim
        write(*,*) itype,isize,cname
        read(lout3d) nx,ny,ndep2
        write(*,*) itype,isize,cname
        read(lout3d) itype,isize,cname    ! xl
        write(*,*) itype,isize,cname
        call close(lout3d)
        nqtot=isize
        write(*,*) 'rout3d: nx,ny,ndep,nqtot=',nx,ny,ndep2,nqtot
        allocate(sp3d(nqtot,nx,ny),stat=ierr)
        if(ierr.ne.0) write(*,*) 'allocation of sp3d(Iv) failed'
        call mopen(lout3d,'Iv3d.'//run_ext(1:lent(run_ext)),0,
     &    'old')
        read(lout3d) itype,isize,cname
        read(lout3d) sp3d
        call close(lout3d)
        count=0
        do iy=1,ny
          do ix=1,nx
            if(sp3d(1,ix,iy).eq.0) then
              icdo(ix,iy)=1
              count=count+1
            else
              icdo(ix,iy)=0
            endif
          enddo
        enddo
        write(*,*) 'rout3d: '//run_ext(1:lent(run_ext))//
     &   ' number of redos:',count
        deallocate(sp3d,stat=ierr)
        if(ierr.ne.0) write(*,*) 'deallocation of sp3d(Iv) failed'
      endif
c
      call mpi_bcast(icdo,nx*ny,MPI_INTEGER,masterid,
     &   MPI_COMM_WORLD,mperr)
      call mpi_bcast(count,1,MPI_INTEGER,masterid,
     &   MPI_COMM_WORLD,mperr)
      call mpi_barrier(MPI_COMM_WORLD,mperr)
c
      end subroutine rout3d0
c
c**********************************************************************
c
      subroutine distrib_columns
c
c  distribute columns over processors
c  ixx,iyy is point in original atmosphere
c  ix,iy is point in calculated atmosphere
c  idxy(ix,iy)  id of processor to handle column ix,iy
c  ixc(ic,j)      ix of column ic of processor j
c  iyc(ic,j)      iy of column ic of processor j
c  ixxc(ic,j)     ixx of column ic of processor j
c  iyyc(ic,j)     iyy of column ic of processor j
c  icxy(j)        number of columns handled by processor j
c  ncxy           maximum of icxy over processors
c
      use var_3d
      implicit none
      integer :: id,ix,iy,ixx,iyy,ierr,ic
c
c  calculate maximum number of columns of any processor
c
      id=0
      ix=0
      icxy(:)=0
      do ixx=ix0,ix1,ixstep
        ix=ix+1
        iy=0
        do iyy=iy0,iy1,iystep
          iy=iy+1
          if(icont_opt.eq.-1) then
            if(icdo(ix,iy).eq.0.) cycle
          endif
          idxy(ix,iy)=id
          icxy(id)=icxy(id)+1
          id=id+1
          if(id.gt.mpnproc-1) id=0
        enddo
      enddo
      ncxy=maxval(icxy)
c
c  allocate variables and set indices
c
      allocate(ixc(ncxy,0:mpnproc-1),iyc(ncxy,0:mpnproc-1),
     & ixxc(ncxy,0:mpnproc-1),iyyc(ncxy,0:mpnproc-1),stat=ierr)
      if(ierr.ne.0) write(*,*) 'allocation of ixc... failed'
      id=0
      ix=0
      ic=1
      do ixx=ix0,ix1,ixstep
        ix=ix+1
        iy=0
        do iyy=iy0,iy1,iystep
          iy=iy+1
          if(icont_opt.eq.-1) then
            if(icdo(ix,iy).eq.0.) cycle
          endif
          ixc(ic,id)=ix
          iyc(ic,id)=iy
          ixxc(ic,id)=ixx
          iyyc(ic,id)=iyy
          id=id+1
          if(id.gt.mpnproc-1) then
            id=0
            ic=ic+1
          endif
        enddo
      enddo
c
      end subroutine distrib_columns
c
c**********************************************************************
c
      subroutine ratmos3d
c
c  read atmos3d and possibly dscal2 and n3d files
c  only one full cube (sp3d) allocated to save memory
c  interpolation to given scale done later
c
      use var_3d
      implicit none
      integer :: latom
      integer :: i,latmos3d,ldsc3d,ln3d,itype,isize,nx2,ny2,ierr,lent
      character(len=8) cname
      character(len=80) text
c
      if(master) then
        call mopen(latmos3d,'atmos3d.'//run_ext(1:lent(run_ext)),0,
     &   'old')
        read(latmos3d) itype,isize,cname
        write(*,*) itype,isize,cname
        read(latmos3d) nx2,ny2,ndep2
        write(*,*) nx2,ny2,ndep2
        allocate(sp3d(nx2,ny2,ndep2),stat=ierr)
        if(ierr.ne.0) write(*,*) 'allocation of sp3d failed'
        allocate(sp2d_tmp(ndep2,ncxy),stat=ierr)
        if(ierr.ne.0) write(*,*) 'allocation of sp2d_tmp failed'
        sp2d_tmp=0.
        allocate(widthx(nx2),stat=ierr)
        if(ierr.ne.0) write(*,*) 'allocation of widthx failed'
        allocate(widthy(ny2),stat=ierr)
        if(ierr.ne.0) write(*,*) 'allocation of widthy failed'
      endif
      call mpi_bcast(ndep2,1,MPI_INTEGER,masterid,MPI_COMM_WORLD,mperr)
c
c  allocate local input variables
c
      if(icxy(mpmyid).gt.0) then
        allocate(sp2d(ndep2,icxy(mpmyid)),stat=ierr)
        if(ierr.ne.0) write(*,*) 'allocation of sp2d failed'
        allocate(temp3d(ndep2,icxy(mpmyid)),stat=ierr)
        if(ierr.ne.0) write(*,*) 'allocation of temp3d failed'
        allocate(ne3d(ndep2,icxy(mpmyid)),stat=ierr)
        if(ierr.ne.0) write(*,*) 'allocation of ne3d failed'
        allocate(vz3d(ndep2,icxy(mpmyid)),stat=ierr)
        if(ierr.ne.0) write(*,*) 'allocation of vz3d failed'
        allocate(rho3d(ndep2,icxy(mpmyid)),stat=ierr)
        if(ierr.ne.0) write(*,*) 'allocation of rho3d failed'
        if(dscopt(1:6).eq.'dscal2') then
          allocate(dpnew3d(ndep2,icxy(mpmyid)),stat=ierr)
          if(ierr.ne.0) write(*,*) 'allocation of dpnew3d failed'
        endif
      endif
      allocate(height_in(ndep2),stat=ierr)
      if(ierr.ne.0) write(*,*) 'allocation of height_in failed'
c
c  read and send
c
      if(master) then
        read(latmos3d) itype,isize,cname
        write(*,*) itype,isize,cname
        read(latmos3d) widthx
        read(latmos3d) itype,isize,cname
        write(*,*) itype,isize,cname
        read(latmos3d) widthy
        read(latmos3d) itype,isize,cname
        write(*,*) itype,isize,cname
        read(latmos3d) height_in
      endif
      call mpi_bcast(height_in,ndep2,MPI_REAL,masterid,
     & MPI_COMM_WORLD,mperr)
      if(master) then
        read(latmos3d) itype,isize,cname
        write(*,*) itype,isize,cname
        read(latmos3d) sp3d      ! ne
      endif
      call scatter_3d_atmos
      if(icxy(mpmyid).gt.0) ne3d=sp2d
      if(master) then
        read(latmos3d) itype,isize,cname
        write(*,*) itype,isize,cname
        read(latmos3d) sp3d      ! temp
      endif
      call scatter_3d_atmos
      if(icxy(mpmyid).gt.0) temp3d=sp2d
      if(master) then
        read(latmos3d) itype,isize,cname
        write(*,*) itype,isize,cname
        read(latmos3d) sp3d      ! vx
        read(latmos3d) itype,isize,cname
        write(*,*) itype,isize,cname
        read(latmos3d) sp3d      ! vy
        read(latmos3d) itype,isize,cname
        write(*,*) itype,isize,cname
        read(latmos3d) sp3d      ! vz
      endif
      call scatter_3d_atmos
      if(icxy(mpmyid).gt.0) vz3d=sp2d
      if(master) then
        read(latmos3d) itype,isize,cname
        write(*,*) itype,isize,cname
        read(latmos3d) sp3d      ! rho
      endif
      call scatter_3d_atmos
      if(icxy(mpmyid).gt.0) rho3d=sp2d
      if(master) then
        deallocate(sp3d,stat=ierr)
        if(ierr.ne.0) write(*,*) 'deallocation of sp3d failed'
        deallocate(sp2d_tmp,stat=ierr)
        if(ierr.ne.0) write(*,*) 'deallocation of sp2d_tmp failed'
        call close(latmos3d)
      endif
c
c  dscal2 file
c
      if(dscopt(1:6).eq.'dscal2') then
        if(master) then
          allocate(sp3d(nx,ny,ndep2),stat=ierr)
          if(ierr.ne.0) write(*,*) 'allocation of sp3d failed'
          allocate(sp2d_tmp(ndep2,ncxy),stat=ierr)
          if(ierr.ne.0) write(*,*) 'allocation of sp2d_tmp failed'
          call mopen(ldsc3d,'dscal2.'//run_ext(1:lent(run_ext)),0,'old')
          read(ldsc3d) itype,isize,cname
          read(ldsc3d) sp3d
          call close(ldsc3d)
        endif
        call scatter_3d
        if(icxy(mpmyid).gt.0) dpnew3d=sp2d
        if(master) then
          deallocate(sp3d,stat=ierr)
          if(ierr.ne.0) write(*,*) 'deallocation of sp3d failed'
          deallocate(sp2d_tmp,stat=ierr)
          if(ierr.ne.0) write(*,*) 'deallocation of sp2d_tmp failed'
        endif
      endif
c
c  n3d file if restart
c
      if(istart2.le.-1) then
        if(master) then
          allocate(sp3d(nx,ny,ndep2),stat=ierr)
          if(ierr.ne.0) write(*,*) 'allocation of sp3d failed'
          allocate(sp2d_tmp(ndep2,ncxy),stat=ierr)
          if(ierr.ne.0) write(*,*) 'allocation of sp2d_tmp failed'
          call mopen(ln3d,'n3d.'//run_ext(1:lent(run_ext)),
     &      -nx*ny*ndep2,'UNKNOWN')
          call mopen(latom,'ATOM',1,'OLD')  ! get nk from ATOM file, 3rd line
          do i=1,3
 100        continue
              read(latom,'(a)') text
            if(text(1:1).eq.'*') goto 100
          enddo
          read(text,*) nk2
          call close(latom)
        endif
        call mpi_bcast(nk2,1,MPI_INTEGER,masterid,MPI_COMM_WORLD,mperr)
        if(icxy(mpmyid).gt.0) then
          allocate(n3d(ndep2,nk2,icxy(mpmyid)),stat=ierr)
          if(ierr.ne.0) write(*,*) 'allocation of n3d failed'
        endif
        do i=1,nk2
          if(master) then
            read(ln3d,rec=i) sp3d
          endif
          call scatter_3d
          if(icxy(mpmyid).gt.0) n3d(:,i,:)=sp2d(:,:)
        enddo
        if(master) then
          call close(ln3d)
          deallocate(sp3d,stat=ierr)
          if(ierr.ne.0) write(*,*) 'deallocation of sp3d failed'
          deallocate(sp2d_tmp,stat=ierr)
          if(ierr.ne.0) write(*,*) 'deallocation of sp2d_tmp failed'
        endif
      endif
c
      end subroutine ratmos3d
c
c**********************************************************************
c
      subroutine scatter_3d_atmos
c
c  splits a 3d variable from global form sp3d on master to local sp2d
c  sp3d is full atmos3d cube so indices are given in ixxc, iyyc
c  sp3d, sp2d_tmp and sp2d must be allocated outside routine
c
      use var_3d
      implicit none
      integer :: id,i,itag,ierr,status(MPI_STATUS_SIZE)
c
      if(master) then             ! sending part
        do id=0,mpnproc-1
          itag=id
          if(id.eq.masterid) then ! master to master, pure copy
            do i=1,icxy(id)
              sp2d(:,i)=sp3d(ixxc(i,id),iyyc(i,id),:)
            enddo
          else if(icxy(id).gt.0) then
            do i=1,icxy(id)
              sp2d_tmp(:,i)=sp3d(ixxc(i,id),iyyc(i,id),:)
            enddo
            call mpi_send(sp2d_tmp,ndep2*icxy(id),
     &       MPI_REAL,id,itag,MPI_COMM_WORLD,mperr)
            if(debug) write(*,*) 'sending (from, to, tag)',
     &       masterid,id,itag
          endif
        enddo
      else if(icxy(mpmyid).gt.0) then                  ! receiving end
        itag=mpmyid
        call mpi_recv(sp2d,ndep2*icxy(mpmyid),MPI_REAL,
     &   masterid,itag,MPI_COMM_WORLD,status,mperr)
        if(debug) write(*,*) 'receiving (from, to, tag)',
     &   masterid,mpmyid,itag
      endif
c
      call mpi_barrier(MPI_COMM_WORLD,mperr)  ! wait for all send/recv to complete
      end subroutine scatter_3d_atmos
c
c**********************************************************************
c
      subroutine scatter_3d
c
c  splits a 3d variable from global form sp3d on master to local sp2d
c  sp3d is calculation cube so indices are given in ixc, iyc
c  sp3d, sp2d_tmp and sp2d must be allocated outside routine
c
      use var_3d
      implicit none
      integer :: id,i,itag,ierr,status(MPI_STATUS_SIZE)
c
      if(master) then             ! sending part
        do id=0,mpnproc-1
          itag=id
          if(id.eq.masterid) then ! master to master, pure copy
            do i=1,icxy(id)
              sp2d(:,i)=sp3d(ixc(i,id),iyc(i,id),:)
            enddo
          else if(icxy(id).gt.0) then
            do i=1,icxy(id)
              sp2d_tmp(:,i)=sp3d(ixc(i,id),iyc(i,id),:)
            enddo
            call mpi_send(sp2d_tmp,ndep2*icxy(id),
     &       MPI_REAL,id,itag,MPI_COMM_WORLD,mperr)
            if(debug) write(*,*) 'sending (from, to, tag)',
     &       masterid,id,itag
          endif
        enddo
      else if(icxy(mpmyid).gt.0) then                  ! receiving end
        itag=mpmyid
        call mpi_recv(sp2d,ndep2*icxy(mpmyid),MPI_REAL,
     &   masterid,itag,MPI_COMM_WORLD,status,mperr)
        if(debug) write(*,*) 'receiving (from, to, tag)',
     &   masterid,mpmyid,itag
      endif
c
      call mpi_barrier(MPI_COMM_WORLD,mperr)  ! wait for all send/recv to complete
      end subroutine scatter_3d
c
c**********************************************************************
c
      subroutine gather_3d
c
c  collects piecec of a 3d variable from local copies into global sp3d
c
      use var_3d
      implicit none
      integer :: id,i,itag,ierr,status(MPI_STATUS_SIZE)
c
      if(master) then             ! receiving part
        do id=0,mpnproc-1
          itag=id
          if(id.eq.masterid) then ! master to master, pure copy
            do i=1,icxy(id)
              sp3d(ixc(i,id),iyc(i,id),:)=sp2d(:,i)
            enddo
          else if(icxy(id).gt.0) then
            call mpi_recv(sp2d_tmp,ndep2*icxy(id),MPI_REAL,
     &       id,itag,MPI_COMM_WORLD,status,mperr)
            do i=1,icxy(id)
              sp3d(ixc(i,id),iyc(i,id),:)=sp2d_tmp(:,i)
            enddo
            if(debug) write(*,*) 'receiving (from, to, tag)',
     &       id,masterid,itag
          endif
        enddo
      else if(icxy(mpmyid).gt.0) then                  ! sending end
        itag=mpmyid
        call mpi_send(sp2d,ndep2*icxy(mpmyid),MPI_REAL,
     &   masterid,itag,MPI_COMM_WORLD,mperr)
        if(debug) write(*,*) 'sending (from, to, tag)',
     &   mpmyid,masterid,itag
      endif
c
      call mpi_barrier(MPI_COMM_WORLD,mperr)  ! wait for all send/recv to complete
      end subroutine gather_3d
c
c**********************************************************************
c
      subroutine gather_Iv
c
c  collects pieces of Iv from local copies into global Iv3d
c
      use var_3d
      implicit none
      integer :: id,i,itag,ierr,status(MPI_STATUS_SIZE)
c
      if(master) then             ! receiving part
        do id=0,mpnproc-1
          itag=id
          if(id.eq.masterid) then ! master to master, pure copy
            do i=1,icxy(id)
              Iv3d(:,ixc(i,id),iyc(i,id))=Iv(:,i)
            enddo
          else if(icxy(id).gt.0) then
            call mpi_recv(Iv2d_tmp,nqtot*icxy(id),MPI_REAL,
     &       id,itag,MPI_COMM_WORLD,status,mperr)
            do i=1,icxy(id)
              Iv3d(:,ixc(i,id),iyc(i,id))=Iv2d_tmp(:,i)
            enddo
            if(debug) write(*,*) 'Iv receiving (from, to, tag)',
     &       id,masterid,itag
          endif
        enddo
      else if(icxy(mpmyid).gt.0) then                  ! sending end
        itag=mpmyid
        call mpi_send(Iv,nqtot*icxy(mpmyid),MPI_REAL,
     &   masterid,itag,MPI_COMM_WORLD,mperr)
        if(debug) write(*,*) 'Iv sending (from, to, tag)',
     &   mpmyid,masterid,itag
      endif
c
      call mpi_barrier(MPI_COMM_WORLD,mperr)  ! wait for all send/recv to complete
      end subroutine gather_Iv
c
c**********************************************************************
c
      subroutine calc_xl
c
c  if(master) fill xl with wavelengths in AA in vacuum
c
      use var_3d
      include 'PREC'
      include 'PARAM'
      include 'CATOM'
      include 'CCONST'
      include 'CSLINE'
      integer :: nu0,kt,ierr,kr,nu,ii
c
      if(master) then
        qn=qnorm*1.d5/cc
        allocate(xl(nqtot),stat=ierr)
        xl(:)=0.0
        if(ierr.ne.0) write(*,*) 'allocation of xl failed'
        nu0=1
        kt=1
        do kr=1,nline
          if(lkr_write(kr)) then
            do nu=1,nq(kr)
              xl(nu0+nu)=alamb(kr)-alamb(kr)*q(nu,kr)*qn/
     &           (1.0d0+q(nu,kr)*qn)
            enddo
            nu0=nu0+nq(kr)+1
            if(iwide(kr)) kt=kt+1
          endif
        enddo
        do kr=nline+1,nrad
          if(lkr_write(kr)) then
            do nu=1,nq(kr)
              xl(nu0+nu)=cc*1.e8/frq(nu,kt)
            enddo
            nu0=nu0+nq(kr)+1
            kt=kt+1
          endif
        enddo
        ii=1
        do kr=1,nrad
          if(lkr_write(kr)) then
            nq_write(ii)=nq(kr)
            ii=ii+1
          endif
        enddo
      endif
c
      end subroutine calc_xl
c
c**********************************************************************
c
      subroutine allocate_out3dvar
c
c  allocate local variables for output to out3d file
c
      use var_3d
      implicit none
      integer :: ierr
c
      if(icxy(mpmyid).gt.0) then
        allocate(Iv(nqtot,icxy(mpmyid)),stat=ierr)
        if(ierr.ne.0) write(*,*) 'allocation of Iv failed'
        Iv(:,:)=0.
        allocate(taulg3d(ndep2,icxy(mpmyid)),stat=ierr)
        if(ierr.ne.0) write(*,*) 'allocation of taulg3d failed'
        allocate(cmasslg3d(ndep2,icxy(mpmyid)),stat=ierr)
        if(ierr.ne.0) write(*,*) 'allocation of cmasslg3d failed'
        allocate(xnorm3d(ndep2,icxy(mpmyid)),stat=ierr)
        if(ierr.ne.0) write(*,*) 'allocation of xnorm3d failed'
        allocate(x3d(ndep2,icxy(mpmyid)),stat=ierr)
        if(ierr.ne.0) write(*,*) 'allocation of x3d failed'
        allocate(height3d(ndep2,icxy(mpmyid)),stat=ierr)
        if(ierr.ne.0) write(*,*) 'allocation of height3d failed'
        if(dscopt(1:6).ne.'dscal2') then
          allocate(dpnew3d(ndep2,icxy(mpmyid)),stat=ierr)
          if(ierr.ne.0) write(*,*) 'allocation of dpnew3d failed'
        endif
        if(istart2.ne.-1) then
          allocate(n3d(ndep2,nk2,icxy(mpmyid)),stat=ierr)
          if(ierr.ne.0) write(*,*) 'allocation of n3d failed'
        endif
        allocate(b3d(ndep2,nk2,icxy(mpmyid)),stat=ierr)
        if(ierr.ne.0) write(*,*) 'allocation of b3d failed'
      endif

      end subroutine allocate_out3dvar
c
c**********************************************************************
c
      subroutine deallocate_3dvar
c
c  dellocate all variables (input and output) that are allocated for
c  each atmos3d run
c
      use var_3d
      implicit none
      integer :: ierr
c
c  variables allocated in distrib_columns
c
      deallocate(ixc,iyc,ixxc,iyyc,stat=ierr)
      if(ierr.ne.0) write(*,*) 'deallocation of ixc... failed'
c
c  variables allocated in ratmos3d
c
      deallocate(height_in)
      if(icxy(mpmyid).gt.0) then
        deallocate(temp3d,ne3d,vz3d,rho3d,dpnew3d,
     &   sp2d,stat=ierr)
        if(ierr.ne.0) write(*,*) 'deallocation of temp3d... failed'
      endif
      if(master) then
        deallocate(widthx,widthy,stat=ierr)
        if(ierr.ne.0) write(*,*) 'deallocation of widthx... failed'
      endif
c
c  variables allocated in calc_xl
c
      if(master) then
        deallocate(xl,stat=ierr)
        if(ierr.ne.0) write(*,*) 'deallocation of xl failed'
      endif
c
c  variables allocated in allocate_out3dvar
c
      if(icxy(mpmyid).gt.0) then
        deallocate(Iv,stat=ierr)
        if(ierr.ne.0) write(*,*) 'deallocation of Iv failed'
        deallocate(taulg3d,stat=ierr)
        if(ierr.ne.0) write(*,*) 'deallocation of taulg3d failed'
        deallocate(cmasslg3d,stat=ierr)
        if(ierr.ne.0) write(*,*) 'deallocation of cmasslg3d failed'
        deallocate(xnorm3d,stat=ierr)
        if(ierr.ne.0) write(*,*) 'deallocation of xnorm3d failed'
        deallocate(x3d,stat=ierr)
        if(ierr.ne.0) write(*,*) 'deallocation of x3d failed'
        deallocate(height3d,stat=ierr)
        if(ierr.ne.0) write(*,*) 'deallocation of height3d failed'
        deallocate(n3d,stat=ierr)
        if(ierr.ne.0) write(*,*) 'deallocation of n3d failed'
        deallocate(b3d,stat=ierr)
        if(ierr.ne.0) write(*,*) 'deallocation of b3d failed'
      endif
c
      end subroutine deallocate_3dvar
c
c**********************************************************************
c
      subroutine wout3d
c
c  write output to file
c
      use var_3d
      integer :: itype,isize,lout,lent
      character(len=8) :: cname
c
      if(master) then
        call mopen(lout,'out3d.'//run_ext(1:lent(run_ext)),0,
     &   'MASTERNEW')
        itype=7
        isize=80
        cname='id      '
        write(lout) itype,isize,cname
        write(*,*) itype,isize,cname
        write(lout) run_ext
        itype=3
        isize=5+nrad_write
        cname='dim     '
        write(lout) itype,isize,cname
        write(*,*) itype,isize,cname
        write(lout) nx,ny,ndep2,mq2,nrad_write,nq_write(1:nrad_write)
        itype=5
        isize=nqtot
        cname='xl      '
        write(lout) itype,isize,cname
        write(*,*) itype,isize,cname
        write(lout) xl(1:nqtot)
        call close(lout)
      endif
c
      end subroutine wout3d
c
c**********************************************************************
c
      subroutine wtaulg3d
c
c  write output to file
c  master collects the data from the slaves and does the writing
c  sp3d is used to collect the data from the slaves local sp2d
c  sp2d_tmp is used for receive
c
      use var_3d
      integer :: itype,isize,lout,ierr,lent
      character(len=8) :: cname
c
c  output to taulg3d file
c
      if(master) then
        allocate(sp3d(nx,ny,ndep2),stat=ierr)
        if(ierr.ne.0) write(*,*) 'allocation of sp3d failed'
        allocate(sp2d_tmp(ndep2,ncxy),stat=ierr)
        if(ierr.ne.0) write(*,*) 'allocation of sp2d_tmp failed'
        sp2d_tmp=0.
        if(icont_opt.eq.-1) then
          call mopen(lout,'taulg3d.'//run_ext(1:lent(run_ext)),0,
     &    'old')
          read(lout) itype,isize,cname
          read(lout) sp3d
          call close(lout)
        endif
      endif
      if(icxy(mpmyid).gt.0) sp2d=taulg3d
      call gather_3d
      if(master) then
        call mopen(lout,'taulg3d.'//run_ext(1:lent(run_ext)),0,
     &   'MASTERNEW')
        itype=4
        isize=nx*ny*ndep2
        cname='taulg3d '
        write(lout) itype,isize,cname
        write(*,*) itype,isize,cname
        write(lout) sp3d
        call close(lout)
        deallocate(sp3d,stat=ierr)
        if(ierr.ne.0) write(*,*) 'deallocation of sp3d failed'
        deallocate(sp2d_tmp,stat=ierr)
        if(ierr.ne.0) write(*,*) 'deallocation of sp2d_tmp failed'
      endif
c
      end subroutine wtaulg3d
c
c**********************************************************************
c
      subroutine wcmass3d
c
c  write output to file
c  master collects the data from the slaves and does the writing
c  sp3d is used to collect the data from the slaves local sp2d
c  sp2d_tmp is used for receive
c
      use var_3d
      integer :: itype,isize,lout,ierr,lent
      character(len=8) :: cname
c
c  output to cmass3d file
c
      if(master) then
        allocate(sp3d(nx,ny,ndep2),stat=ierr)
        if(ierr.ne.0) write(*,*) 'allocation of sp3d failed'
        allocate(sp2d_tmp(ndep2,ncxy),stat=ierr)
        if(ierr.ne.0) write(*,*) 'allocation of sp2d_tmp failed'
        sp2d_tmp=0.
        if(icont_opt.eq.-1) then
          call mopen(lout,'cmass3d.'//run_ext(1:lent(run_ext)),0,
     &    'old')
          read(lout) itype,isize,cname
          read(lout) sp3d
          call close(lout)
        endif
      endif
      if(icxy(mpmyid).gt.0) sp2d=cmasslg3d
      call gather_3d
      if(master) then
        call mopen(lout,'cmass3d.'//run_ext(1:lent(run_ext)),0,
     &   'MASTERNEW')
        itype=4
        isize=nx*ny*ndep2
        cname='cmass3d  '
        write(lout) itype,isize,cname
        write(*,*) itype,isize,cname
        write(lout) sp3d
        call close(lout)
        deallocate(sp3d,stat=ierr)
        if(ierr.ne.0) write(*,*) 'deallocation of sp3d failed'
        deallocate(sp2d_tmp,stat=ierr)
        if(ierr.ne.0) write(*,*) 'deallocation of sp2d_tmp failed'
      endif
c
      end subroutine wcmass3d
c
c**********************************************************************
c
      subroutine wxnorm3d
c
c  write output to file
c  master collects the data from the slaves and does the writing
c  sp3d is used to collect the data from the slaves local sp2d
c  sp2d_tmp is used for receive
c
      use var_3d
      integer :: itype,isize,lout,ierr,lent
      character(len=8) :: cname
c
c  output to xnorm3d file
c
      if(master) then
        allocate(sp3d(nx,ny,ndep2),stat=ierr)
        if(ierr.ne.0) write(*,*) 'allocation of sp3d failed'
        allocate(sp2d_tmp(ndep2,ncxy),stat=ierr)
        if(ierr.ne.0) write(*,*) 'allocation of sp2d_tmp failed'
        sp2d_tmp=0.
        if(icont_opt.eq.-1) then
          call mopen(lout,'xnorm3d.'//run_ext(1:lent(run_ext)),0,
     &    'old')
          read(lout) itype,isize,cname
          read(lout) sp3d
          call close(lout)
        endif
      endif
      if(icxy(mpmyid).gt.0) sp2d=xnorm3d
      call gather_3d
      if(master) then
        call mopen(lout,'xnorm3d.'//run_ext(1:lent(run_ext)),0,
     &   'MASTERNEW')
        itype=4
        isize=nx*ny*ndep2
        cname='xnorm3d '
        write(lout) itype,isize,cname
        write(*,*) itype,isize,cname
        write(lout) sp3d
        call close(lout)
        deallocate(sp3d,stat=ierr)
        if(ierr.ne.0) write(*,*) 'deallocation of sp3d failed'
        deallocate(sp2d_tmp,stat=ierr)
        if(ierr.ne.0) write(*,*) 'deallocation of sp2d_tmp failed'
      endif
c
      end subroutine wxnorm3d
c
c**********************************************************************
c
      subroutine wx3d
c
c  write output to file
c  master collects the data from the slaves and does the writing
c  sp3d is used to collect the data from the slaves local sp2d
c  sp2d_tmp is used for receive
c
      use var_3d
      integer :: itype,isize,lout,ierr,lent
      character(len=8) :: cname
c
c  output to x3d file
c
      if(master) then
        allocate(sp3d(nx,ny,ndep2),stat=ierr)
        if(ierr.ne.0) write(*,*) 'allocation of sp3d failed'
        allocate(sp2d_tmp(ndep2,ncxy),stat=ierr)
        if(ierr.ne.0) write(*,*) 'allocation of sp2d_tmp failed'
        sp2d_tmp=0.
        if(icont_opt.eq.-1) then
          call mopen(lout,'x3d.'//run_ext(1:lent(run_ext)),0,
     &    'old')
          read(lout) itype,isize,cname
          read(lout) sp3d
          call close(lout)
        endif
      endif
      if(icxy(mpmyid).gt.0) sp2d=x3d
      call gather_3d
      if(master) then
        call mopen(lout,'x3d.'//run_ext(1:lent(run_ext)),0,
     &   'MASTERNEW')
        itype=4
        isize=nx*ny*ndep2
        cname='x3d     '
        write(lout) itype,isize,cname
        write(*,*) itype,isize,cname
        write(lout) sp3d
        call close(lout)
        deallocate(sp3d,stat=ierr)
        if(ierr.ne.0) write(*,*) 'deallocation of sp3d failed'
        deallocate(sp2d_tmp,stat=ierr)
        if(ierr.ne.0) write(*,*) 'deallocation of sp2d_tmp failed'
      endif
c
      end subroutine wx3d
c
c**********************************************************************
c
      subroutine wheight3d
c
c  write output to file
c  master collects the data from the slaves and does the writing
c  sp3d is used to collect the data from the slaves local sp2d
c  sp2d_tmp is used for receive
c
      use var_3d
      integer :: itype,isize,lout,ierr,lent
      character(len=8) :: cname
c
c  output to height3d file
c
      if(master) then
        allocate(sp3d(nx,ny,ndep2),stat=ierr)
        if(ierr.ne.0) write(*,*) 'allocation of sp3d failed'
        allocate(sp2d_tmp(ndep2,ncxy),stat=ierr)
        if(ierr.ne.0) write(*,*) 'allocation of sp2d_tmp failed'
        sp2d_tmp=0.
        if(icont_opt.eq.-1) then
          call mopen(lout,'height3d.'//run_ext(1:lent(run_ext)),0,
     &    'old')
          read(lout) itype,isize,cname
          read(lout) sp3d
          call close(lout)
        endif
      endif
      if(icxy(mpmyid).gt.0) sp2d=height3d
      call gather_3d
      if(master) then
        call mopen(lout,'height3d.'//run_ext(1:lent(run_ext)),0,
     &   'MASTERNEW')
        itype=4
        isize=nx*ny*ndep2
        cname='height3d '
        write(lout) itype,isize,cname
        write(*,*) itype,isize,cname
        write(lout) sp3d
        call close(lout)
        deallocate(sp3d,stat=ierr)
        if(ierr.ne.0) write(*,*) 'deallocation of sp3d failed'
        deallocate(sp2d_tmp,stat=ierr)
        if(ierr.ne.0) write(*,*) 'deallocation of sp2d_tmp failed'
      endif
c
      end subroutine wheight3d
c
c**********************************************************************
c
      subroutine wIv3d
c
c  write output to file
c  master collects the data from the slaves and does the writing
c  sp3d is used to collect the data from the slaves local sp2d
c  sp2d_tmp is used for receive
c
      use var_3d
      integer :: itype,isize,lout,ierr,lent
      character(len=8) :: cname
c
      if(master) then
        allocate(Iv3d(nqtot,nx,ny),stat=ierr)
        if(ierr.ne.0) write(*,*) 'allocation of Iv3d failed'
        allocate(Iv2d_tmp(nqtot,ncxy),stat=ierr)   ! local receive variable
        if(ierr.ne.0) write(*,*) 'allocation of Iv2d_tmp failed'
        Iv2d_tmp=0.
        if(icont_opt.eq.-1) then
          call mopen(lout,'Iv3d.'//run_ext(1:lent(run_ext)),0,
     &    'old')
          read(lout) itype,isize,cname
          read(lout) Iv3d
          call close(lout)
        endif
      endif
      call gather_Iv
      if(master) then
        call mopen(lout,'Iv3d.'//run_ext(1:lent(run_ext)),0,
     &   'MASTERNEW')
        itype=4
        isize=nx*ny*nqtot
        cname='Iv      '
        write(lout) itype,isize,cname
        write(*,*) itype,isize,cname
        write(lout) Iv3d
        deallocate(Iv3d,stat=ierr)              ! Iv written
        if(ierr.ne.0) write(*,*) 'deallocation of Iv3d failed'
        deallocate(Iv2d_tmp,stat=ierr)
        if(ierr.ne.0) write(*,*) 'deallocation of Iv2d_tmp failed'
        call close(lout)
      endif
c
      end subroutine wIv3d
c
c**********************************************************************
c
      subroutine wdscal2
c
c  write output to dscal2 file
c  master collects the data from the slaves and does the writing
c  sp3d is used to collect the data from the slaves local sp2d
c  sp2d_tmp is used for receive
c
      use var_3d
      integer :: itype,isize,ierr,lent,lout
      character(len=8) :: cname
c
c  output to dscale file
c
      if(master) then
        allocate(sp3d(nx,ny,ndep2),stat=ierr)
        if(ierr.ne.0) write(*,*) 'allocation of sp3d failed'
        allocate(sp2d_tmp(ndep2,ncxy),stat=ierr)
        if(ierr.ne.0) write(*,*) 'allocation of sp2d_tmp failed'
        sp2d_tmp=0.
        if(icont_opt.eq.-1) then
          call mopen(lout,'dscal2.'//run_ext(1:lent(run_ext)),0,
     &    'old')
          read(lout) itype,isize,cname
          read(lout) sp3d
          call close(lout)
        endif
      endif
      if(icxy(mpmyid).gt.0) sp2d=dpnew3d
      call gather_3d
      if(master) then
        call mopen(lout,'dscal2.'//run_ext(1:lent(run_ext)),0,
     &   'MASTERNEW')
        itype=4
        isize=nx*ny*ndep2
        cname='dscal2  '
        write(lout) itype,isize,cname
        write(*,*) itype,isize,cname
        write(lout) sp3d
        call close(lout)
        deallocate(sp3d,stat=ierr)              ! dpnew3d written
        if(ierr.ne.0) write(*,*) 'deallocation of sp3d failed'
        deallocate(sp2d_tmp,stat=ierr)
        if(ierr.ne.0) write(*,*) 'deallocation of sp2d_tmp failed'
      endif
c
      end subroutine wdscal2
c
c**********************************************************************
c
      subroutine wn3d
c
c  write output to n3d file
c  master collects the data from the slaves and does the writing
c  sp3d is used to collect the data from the slaves local sp2d
c  sp2d_tmp is used for receive
c
      use var_3d
      integer :: itype,isize,ierr,lent,i
      integer :: ln3d,ln3d0,nk0
      character(len=8) :: cname
c
c  output to n3d file
c
      if(master) then
        call mopen(ln3d,'n3d.'//run_ext(1:lent(run_ext)),-nx*ny*ndep2,
     &   'MASTERNEW')
        allocate(sp3d(nx,ny,ndep2),stat=ierr)
        if(ierr.ne.0) write(*,*) 'allocation of sp3d failed'
        allocate(sp2d_tmp(ndep2,ncxy),stat=ierr)
        if(ierr.ne.0) write(*,*) 'allocation of sp2d_tmp failed'
      endif
      call mpi_bcast(nk2,1,MPI_INTEGER,masterid,
     &   MPI_COMM_WORLD,mperr)
      do i=1,nk2
        if(master .and. icont_opt.eq.-1) then
          read(ln3d,rec=i) sp3d
        endif
        if(icxy(mpmyid).gt.0) sp2d(:,:)=n3d(:,i,:)
        call gather_3d
        if(master) then
          write(*,*) 4,nx*ny*ndep2,'n3d, i=',i
          write(ln3d,rec=i) sp3d
        endif
      enddo
      if(master) then
        call close(ln3d)
        deallocate(sp3d,stat=ierr)              ! n3d written
        if(ierr.ne.0) write(*,*) 'deallocation of sp3d failed'
        deallocate(sp2d_tmp,stat=ierr)
        if(ierr.ne.0) write(*,*) 'deallocation of sp2d_tmp failed'
      endif
c
      end subroutine wn3d
c
c**********************************************************************
c
      subroutine wb3d
c
c  write output to b3d file
c  master collects the data from the slaves and does the writing
c  sp3d is used to collect the data from the slaves local sp2d
c  sp2d_tmp is used for receive
c
      use var_3d
      integer :: itype,isize,ierr,lent,i
      integer :: lb3d,lb3d0,nk0
      character(len=8) :: cname
c
c  output to b3d file
c
      if(master) then
        call mopen(lb3d,'b3d.'//run_ext(1:lent(run_ext)),-nx*ny*ndep2,
     &   'MASTERNEW')
        allocate(sp3d(nx,ny,ndep2),stat=ierr)
        if(ierr.ne.0) write(*,*) 'allocation of sp3d failed'
        allocate(sp2d_tmp(ndep2,ncxy),stat=ierr)
        if(ierr.ne.0) write(*,*) 'allocation of sp2d_tmp failed'
      endif
      call mpi_bcast(nk2,1,MPI_INTEGER,masterid,
     &   MPI_COMM_WORLD,mperr)
      do i=1,nk2
        if(master .and. icont_opt.eq.-1) then
          read(lb3d,rec=i) sp3d
        endif
        if(icxy(mpmyid).gt.0) sp2d(:,:)=b3d(:,i,:)
        call gather_3d
        if(master) then
          write(*,*) 4,nx*ny*ndep2,'b3d, i=',i
          write(lb3d,rec=i) sp3d
        endif
      enddo
      if(master) then
        call close(lb3d)
        deallocate(sp3d,stat=ierr)              ! b3d written
        if(ierr.ne.0) write(*,*) 'deallocation of sp3d failed'
        deallocate(sp2d_tmp,stat=ierr)
        if(ierr.ne.0) write(*,*) 'deallocation of sp2d_tmp failed'
      endif
c
      end subroutine wb3d
c
c**********************************************************************
c
      character*3 function string3(i)
c
      if(i.lt.10) then 
        write(string3,'(A2,i1)') '00',i
      else if(i.lt.100) then
        write(string3,'(A1,i2)') '0',i
      else 
        write(string3,'(i3)') i
      endif

      return
      end
c
c**********************************************************************
c
      function lent(text)
c
c  finds the index of the last non-blank character of a text-string
c
      character*(*) text
c
      l=len(text)
      do 10 i=l,1,-1
        if(text(i:i).ne.' ') goto 100
   10 continue
  100 continue
      lent=i
c
      end
c
c**********************************************************************
c
      subroutine ipol_dscale(ndep2,height_in,temp_in,ne_in,vz_in,rho_in)
c
c  finds a new more optimal depth-scale and interpolates input atmosphere to it
c  use multi variables for output
c
      include 'PREC'
      include 'PARAM'
      include 'CATMOS'
      include 'CATMO2'
      include 'CCONST'
c
      integer ndep2
      real(4) height_in(ndep2),temp_in(ndep2),ne_in(ndep2),vz_in(ndep2)
      real(4) rho_in(ndep2)
      dimension aind2(mdep),f(mdep),x(mdep),thmbf(mdep),xhmbf(mdep)
      logical debug
c
      debug=.false.
      ttop=5.d4
      tbottom=1.2d4 ! not used now
      taumax=100
      grph=2.26e-24
c
c  calculate tau scale from H-bf opacity alone
c  assume all hydrogen to be neutral
c
      thmbf(1)=1.e-6
      crhmbf=2.9256e-17
      do k=1,ndep2
        xhmbf(k)=1.03526d-16*ne_in(k)*crhmbf/temp_in(k)**1.5*
     &  exp(0.754d0*ee/bk/temp_in(k))*rho_in(k)/grph
      enddo
      do k=2,ndep2
        thmbf(k)=thmbf(k-1)+0.5*(xhmbf(k)+xhmbf(k-1))*
     &   abs(height_in(k-1)-height_in(k))
      enddo
c
c  set atmospheric range from temperature of ttop and tau(hmbf)=taumax
c
      k0=1
      k1=1
      aind2(:)=0.0d0
      do k=2,ndep2
        if(temp_in(k).gt.ttop .and. k0.eq.k-1) k0=k
c        if(temp_in(k).lt.tbottom) k1=k  ! old criterion
        if(thmbf(k).lt.taumax) k1=k
      enddo
c
      do k=k0+1,k1
        tdiv=abs(log10(temp_in(k))-log10(temp_in(k-1)))/log10(1.1)
        rdiv=abs(log10(rho_in(k))-log10(rho_in(k-1)))/log10(1.1)
        taudiv=(log10(thmbf(k))-log10(thmbf(k-1)))/0.1
c        write(*,'(i4,3f10.4)') k,tdiv,rdiv,taudiv
        aind2(k)=aind2(k-1)+max(max(tdiv,rdiv),taudiv)
      enddo
c
      aind2=aind2/aind2(k1)*(ndep2-1)
c
c  make new height scale and interpolate atmosphere
c
      nn=k1-k0+1
      f(1:nn)=height_in(k0:k1)
      do k=1,ndep2
        xp=(k-1)
        call intep(k,aind2(k0:k1),f,nn,iup,xp,height(k))
c        if(k.gt.1) then
c          if(height(k).gt.height(k-1)) debug=.true.
c        endif
      enddo
      f(1:nn)=alog(temp_in(k0:k1))
      x(1:nn)=height_in(k0:k1)
      do 400 k=1,ndep2
        xp=height(k)
        call intep(k,x,f,nn,iup,xp,temp(k))
        temp(k)=exp(temp(k))
        if(temp(k).lt.1.d3) then
          do kk=2,nn
            if((xp-x(kk-1))*
     &         (xp-x(kk)) .le. 0.d0) goto 350
          enddo
  350     continue
          temp(k)=f(kk-1)+(xp-x(kk-1))/
     &     (x(kk)-x(kk-1))*(f(kk)-f(kk-1))
          temp(k)=exp(temp(k))
          write(*,*) 'linear interpolation, temp(k)=',temp(k)
        endif  
  400 continue
      f(1:nn)=alog(ne_in(k0:k1))
      do 500 k=1,ndep2
        xp=height(k)
        call intep(k,x,f,nn,iup,xp,ne(k))
        ne(k)=exp(ne(k))
  500 continue
      f(1:nn)=vz_in(k0:k1)
      do 600 k=1,ndep2
        xp=height(k)
        call intep(k,x,f,nn,iup,xp,vel(k))
  600 continue
      f(1:nn)=alog(rho_in(k0:k1))
      do 700 k=1,ndep2
        xp=height(k)
        call intep(k,x,f,nn,iup,xp,rho(k))
        rho(k)=exp(rho(k))
  700 continue
c
      if(debug) then
        write(*,'(a)') 'input variables'
        write(*,'(a)') '  k       height        temp          ne'//
     &   '          vz         rho'
        do 750 k=1,ndep2
          write(*,'(i3,1p,5e12.4)') k,height_in(k),temp_in(k),ne_in(k),
     &     vz_in(k),rho_in(k)
  750   continue
        write(*,*) ' '
        write(*,'(a)') 'interpolation array'
        write(*,'(a)') '  k       aind2'
        do 760 k=k0,k1
          write(*,'(i3,1p,e12.4)') k,aind2(k)
  760   continue
        write(*,*) ' '
        write(*,'(a)') 'interpolated variables'
        write(*,'(a)') '  k       height        temp          ne'//
     &   '          vz         rho'
        do 800 k=1,ndep2
          write(*,'(i3,1p,5e12.4)') k,height(k),temp(k),ne(k),
     &     vel(k),rho(k)
  800   continue
      endif
c
      end
c
c*****************************************************************
c
      subroutine ipol_dscal2(ndep2,height_in,temp_in,ne_in,vz_in,rho_in,
     & height2)
c
c  interpolate to height2 scale
c  use multi variables for output
c
      include 'PREC'
      include 'PARAM'
      include 'CATMOS'
      include 'CATMO2'
c
      integer ndep2
      real(4) height_in(ndep2),temp_in(ndep2),ne_in(ndep2),vz_in(ndep2)
      real(4) rho_in(ndep2),height2(ndep2)
      dimension f(mdep),x(mdep)
      logical debug
c
      debug=.false.
c
c  make new height scale and interpolate atmosphere
c  height_in and height2 in cm, converted to km when written to ATMOS_ files
c
      height(1:ndep2)=height2(1:ndep2)*1.e5
      f(1:ndep2)=alog(temp_in(1:ndep2))
      x(1:ndep2)=height_in(1:ndep2)
      do 400 k=1,ndep2
        xp=height(k)
        call intep(k,x,f,ndep2,iup,xp,temp(k))
        temp(k)=exp(temp(k))
  400 continue
      f(1:ndep2)=alog(ne_in(1:ndep2))
      do 500 k=1,ndep2
        xp=height(k)
        call intep(k,x,f,ndep2,iup,xp,ne(k))
        ne(k)=exp(ne(k))
  500 continue
      f(1:ndep2)=vz_in(1:ndep2)
      do 600 k=1,ndep2
        xp=height(k)
        call intep(k,x,f,ndep2,iup,xp,vel(k))
  600 continue
      f(1:ndep2)=alog(rho_in(1:ndep2))
      do 700 k=1,ndep2
        xp=height(k)
        call intep(k,x,f,ndep2,iup,xp,rho(k))
        rho(k)=exp(rho(k))
  700 continue
c
      if(debug) then
        write(*,'(a)') 'input variables'
        write(*,'(a)') '  k       height        temp          ne'//
     &   '          vz         rho'
        do 750 k=1,ndep2
          write(*,'(i3,1p,5e12.4)') k,height_in(k),temp_in(k),ne_in(k),
     &     vz_in(k),rho_in(k)
  750   continue
        write(*,*) ' '
        write(*,'(a)') 'interpolated variables'
        write(*,'(a)') '  k       height        temp          ne'//
     &   '          vz         rho'
        do 800 k=1,ndep2
          write(*,'(i3,1p,5e12.4)') k,height(k),temp(k),ne(k),
     &     vel(k),rho(k)
  800   continue
      endif
c
      end
c
c*****************************************************************
c
      subroutine ipol_n(ndep2,nk2,height_in,n_in)
c
c  interpolate n from height_in to height scale
c  use multi variables for output
c
      include 'PREC'
      include 'PARAM'
      include 'CATOM'
      include 'CATMOS'
      include 'CATMO2'
c
      integer ndep2,nk2
      real(4) height_in(ndep2),n_in(ndep2,nk2)
      dimension y(mdep),x(mdep)
      logical debug
c
      debug=.false.
c
c  make new height scale and interpolate atmosphere
c
      x(1:ndep2)=height_in(1:ndep2)
      do 500 i=1,nk2
        y(1:ndep2)=alog(n_in(1:ndep2,i))
        do 400 k=1,ndep2
          xp=height(k)
          call intep(k,x,y,ndep2,iup,xp,n(i,k))
          n(i,k)=exp(n(i,k))
  400   continue
 500  continue
c
      if(debug) then
        write(*,'(a)') 'input variables'
        write(*,'(a)') '  k       height           n'
        do 750 k=1,ndep2
          write(*,'(i3,1p,7e10.2)') k,height_in(k),n_in(k,:)
  750   continue
        write(*,*) ' '
        write(*,'(a)') 'interpolated variables'
        write(*,'(a)') '  k       height           n'
        do 800 k=1,ndep2
          write(*,'(i3,1p,7e10.2)') k,height(k),n(1:nk2,k)
  800   continue
      endif
c
      end
c
c*****************************************************************
c
      subroutine intep(it,x,f,n,iup,xp,p)
c
c  interpolation routine based on hermite polynomials
c  ref: publications of the dominion astrophysical observatory, xvi,6,67
c       graham hill: intep, an effective interpolation subroutine
c
      include 'PREC'
      integer it,n,iup
      real(8) lp1,lp2,l1,l2
      dimension f(n),x(n)
      save
c
      if(it.eq.1) then
        ier=1
        io=1
        iup=0
        if(x(2).lt.x(1)) iup=1
        n1=n-1
        if((xp.ge.x(n).and.iup.eq.0).or.(xp.le.x(n).and.iup.eq.1)) then
          p=f(n)
          ier=2
          return
        else if((xp.le.x(1).and.iup.eq.0).or.
     &    (xp.ge.x(1).and.iup.eq.1)) then
          p=f(1)
          ier=2
          return
        endif
      endif
    8 do 1 i=io,n
        if(xp.lt.x(i).and.iup.eq.0) goto 2
        if(xp.gt.x(i).and.iup.eq.1) goto 2
    1 continue
      p=f(n)
      ier=2
      return
    2 i=i-1
      if(i.eq.io-1) goto 4
      io=i+1
      lp1=1.d0/(x(i)-x(i+1))
      lp2=1.d0/(x(i+1)-x(i))
      if(i.eq.1) fp1=(f(2)-f(1))/(x(2)-x(1))
      if(i.eq.1) goto 3
      fp1=(f(i+1)-f(i-1))/(x(i+1)-x(i-1))
    3 if(i.ge.n1) fp2=(f(n)-f(n-1))/(x(n)-x(n-1))
      if(i.ge.n1) goto 4
      fp2=(f(i+2)-f(i))/(x(i+2)-x(i))
    4 xpi1=xp-x(i+1)
      xpi=xp-x(i)
      l1=xpi1*lp1
      l2=xpi*lp2
      p=f(i)*(1.d0-2.d0*lp1*xpi)*l1*l1+f(i+1)*(1.d0-2.d0*lp2*xpi1)*
     & l2*l2+fp2*xpi1*l2*l2+fp1*xpi*l1*l1
c
      return
      end   
c
c*********************************************************************
c
      subroutine disable_printout
c
c  disable printout to avoid excessive printout
c
      INCLUDE 'CINPUT'

      iwabnd=0
      iwatms=0
      iwatom=0
      iwchan=0
      iwdamp=0
      iwemax=0
      iweqw=0
      iwevec=0
      iwhead=0
      iwhse=0
      iwlgmx=0
      iwline=0
      iwlte=0
      iwn=0
      iwniit=0
      iwopac=0
      iwrad=0
      iwrate=0
      iwstrt=0
      iwtauq=0
      iwtest=0
      iwwmat=0
      idl1=0
      idlny=0
      idlcnt=0
      idlopc=0

      end
