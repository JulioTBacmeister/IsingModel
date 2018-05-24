program linsys
  implicit none
  integer, parameter :: N=4
  integer :: OK,i
  real :: AA(N,N),BB(N),KK0(N),KK1(N)
  real :: AAx(N,N)
  integer :: pivot(N)

  AAx(:,:) = 0.
  do i=1,N
     AAx( i : N , i ) = 1.
     BB(i) = i*1.
  end do
  
#if 0
  AA = AAx
  AAx(:,N) = AA(:,1)
  AAx(:,1)   = AA(:,N)
  BB(N)=(1)*1.
  BB(1)=(n)*1.
#endif
  
     !!!!!!!!!!!!!!!!!****************
     !!!!!!!!!!!!!!!!!!!!!************
     !  NOTE that LAPACK subr's seem
     !  to need *TRANSPOSE* of the way A
     !  is naturally written above
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     AA = transpose(AAx)  !!!!!!!!!!!!
     !AA = AAx ! transpose(AAx)  !!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          call DGESV(N, 1, AA, N, pivot, BB, N, ok)
          write(*,*) BB
          write(*,*) pivot
          

end program linsys
