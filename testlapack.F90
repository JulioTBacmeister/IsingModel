program testlapack
real :: a(3,3),b(3),bb(3)
integer :: pivot(3), OK,j


a(1,1)=3.
a(1,2)=1.
a(1,3)=0.

a(2,1)=0.
a(2,2)=2.
a(2,3)=0.

a(3,1)=0.
a(3,2)=1.
a(3,3)=3.

b(:)=(/ 0.,1.,2. /)
bb=b


write(61) A,B

do j=1,3
   write(6,'(3f6.2)' ) A(:,J)
end do
!! call DGETRS( 'N' , N, NRHS, A, LDA, IPIV, B, LDB, INFO	)
!! find the solution using the LAPACK routine SGESV
	call DGESV(3, 1, A, 3, pivot, b, 3, ok)
 write(61) B
 

end program testlapack
