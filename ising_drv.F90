program ising_drv
  implicit none
  real    :: x1,x2,x3,y1,y2,y3
  integer :: niter,i
  !integer, allocatable :: CCC(:,:,:)
  !real, allocatable :: ZZZc(:,:,:)
  integer :: CCC(4,4,2**16),NX,NC
  real :: ZZZc(2**16)

     
  write(*,*) " x1,x2,x3 , NC , NX"
  read(*,*) x1,x2,x3,NC,NX
  write(*,*) x1,x2,x3 ,NC,NX


 
     !call rng_sub( x1,x2,x3,y1,y2,y3 )
     !call ising_4x4_partfun( x1,x2,x3,ZZZc,CCC)
     call ising_MC_partfun( NC,NX,x1,x2,x3)  !,ZZZc,CCC)

end program ising_drv
