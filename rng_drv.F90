program rng_drv
  implicit none
  real    :: x1,x2,x3,y1,y2,y3
  integer :: niter,i

     
  write(*,*) " x1,x2,x3,niter "
  read(*,*) x1,x2,x3,niter
  write(*,*) x1,x2,x3 


  
  do i=1,niter
     !call rng_sub( x1,x2,x3,y1,y2,y3 )
     call rng_sub2( x1,x2,x3,y1,y2,y3 )
     x1=y1
     x2=y2
     x3=y3
  end do

end program rng_drv
