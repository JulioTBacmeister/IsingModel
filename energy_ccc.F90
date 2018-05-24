function energy_ccc(ccc,x1,x2,x3,NX) result (energy)
  implicit none
  integer, intent(IN) :: NX
  integer, intent(IN) :: ccc(0:NX+1,0:NX+1)
  real, intent(IN) :: x1,x2,x3

  integer :: a,b,c,d,i,j,n
  real :: NN,NNN,NNNN,Energy
  !======================================================
  ! Input coupling coefficients are assumed to be in the
  ! form K = J/(kb*T) where J's come from Hamiltonian
  !
  !       H = SUM_(ij=NN) -J*S_i*S_j
  !
  ! Critical value: K=0.4407... (analytical, Onsager)
  !=====================================================
  !====================================
  ! Energy of initial spin configuration
  !=====================================   
  ! 2x2 block
  !---------
  ! a b 
  ! c d 
  !---------
     NN   = 0.
     NNN  = 0.
     NNNN = 0.
     do j=1,NX-1
     do i=1,NX-1
        a = ccc(i,j)
        b = ccc(i+1,j)
        c = ccc(i,j+1)
        d = ccc(i+1,j+1)

        NN   = NN    + a*b + a*c
        NNN  = NNN   + a*d + b*c
        NNNN = NNNN  + a*b*c+d
     end do
     end do
     ! Periodic Boundaries
     do j=1,NX-1
     i=NX
        a = ccc(i,j)
        b = ccc(1,j)
        c = ccc(i,j+1)
        NN   = NN    + a*c + a*b
     end do
     j=NX
     do i=1,NX-1
        a = ccc(i,j)
        b = ccc(i+1,j)
        c = ccc(i,1)
        NN   = NN    + a*b + a*c
     end do
     j=NX
     i=NX
        a = ccc(i,j)
        b = ccc(1,j)
        c = ccc(i,1)
        NN   = NN    + a*b + a*c

  Energy = x1*NN + x2*NNN + x3*NNNN

  
end function energy_ccc
