      MODULE var_3d
c
c  module containing wrapper 3d varibles
c
c  mpi variables
c
      include 'mpif.h'
      integer masterid,mpmyid,mpnproc,mperr,nx,ny
      parameter (masterid=0)
      logical master, slave
      character*80 dscopt,ext
c
c  3d wrapper variables
c
      real(4), dimension(:,:,:), allocatable :: sp3d     ! full 3d cube for read/write
      real(4), dimension(:,:),   allocatable :: sp2d     ! local cube for send/receive
      real(4), dimension(:,:),   allocatable :: sp2d_tmp ! temporary cube for send/receive
      real(4), dimension(:,:,:), allocatable :: Iv3d     ! full 3d cube for read/write
      real(4), dimension(:,:),   allocatable :: Iv2d     ! local cube for send/receive
      real(4), dimension(:,:),   allocatable :: Iv2d_tmp ! temporary cube for send/receive
      real(4), dimension(:,:,:), allocatable :: n3d      ! restart/output
      real(4), dimension(:,:,:), allocatable :: b3d      ! restart/output
      real(4), dimension(:,:),   allocatable :: n2d      ! local cube
      real(4), dimension(:,:),   allocatable :: temp3d,ne3d,vz3d,rho3d   ! atmos3d variables
      real(4), dimension(:),     allocatable :: widthx,widthy,height_in  ! atmos3d
      real(4), dimension(:,:), allocatable :: taulg3d,cmasslg3d,xnorm3d  ! output
      real(4), dimension(:,:), allocatable :: height3d,dpnew3d,x3d,Iv ! output
      real(8), dimension(:),     allocatable :: xl
      integer, dimension(:,:), allocatable :: icdo,idxy         ! column admin variables
      integer, dimension(:,:), allocatable :: ixc,iyc,ixxc,iyyc
      integer, dimension(:), allocatable :: icxy    ! number of columns per processor
      integer :: ncxy                               ! max number of columns per processor
      integer :: nk2,ndep2,istart2,mq2              ! transfer of variables from multi
      external string3
      character*3 string3
      logical :: debug

c INPUT_3D variables (xsave,krsave,nusave,l3d in C3D)
      integer :: icont_opt,iredo_max,ix0,ix1,ixstep,iy0,iy1,iystep,nruns
      integer :: nrad_write,nqtot
      integer, dimension(:), allocatable :: kr_write,nq_write
      logical, dimension(:), allocatable :: lkr_write
      character*80, dimension(:), allocatable :: bext
c
      character*80 :: run_ext

      end module var_3d
