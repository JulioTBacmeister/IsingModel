program ising_NxN_drv
  implicit none
  real    :: x1,x2,x3,y1,y2,y3,t
  integer :: niter,i
  !integer, allocatable :: CCC(:,:,:)
  !real, allocatable :: ZZZc(:,:,:)
  integer :: CCC(4,4,2**16),NX,NC
  real :: ZZZc(2**16)

     
  !write(*,*) " x1,x2,x3 , NX"
  !read(*,*) x1,x2,x3,NX
  !write(*,*) x1,x2,x3 ,NX

  NX = 5
  x2=0.
  x3=0.

  do i=1,15
  !do i=8,8
     T=0.5*i
     x1=1./T
     !call rng_sub( x1,x2,x3,y1,y2,y3 )
     !call ising_4x4_partfun( x1,x2,x3,ZZZc,CCC)
     call ising_NxN_partfun( NX,x1,x2,x3)  !,ZZZc,CCC)
  end do
     
   end program ising_NxN_drv
