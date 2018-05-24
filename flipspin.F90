function flipspin(cc,x1,x2,x3,rnd) result (fspin)
  implicit none
  integer, intent(IN) :: cc(-1:1,-1:1)
  real, intent(IN) :: x1,x2,x3,rnd
  integer :: fspin

  real :: NNi,NNNi,NNNNi,EEi
  real :: NNo,NNNo,NNNNo,EEo,DEE
  
  integer :: spinsum
  !======================================================
  ! Input coupling coefficients are assumed to be in the
  ! form K = J/(kb*T) where J's come from Hamiltonian
  !
  !       H = SUM_(ij=NN) -J*S_i*S_j
  !
  ! Critical value: K=0.4407... (analytical, Onsager)
  !=====================================================

  NNi = x1*( cc(0,0)*cc(1,0) + cc(0,0)*cc(-1,0) &
         + cc(0,0)*cc(0,1) + cc(0,0)*cc(0,-1) )

  EEi = -NNi
  EEo = NNi

  DEE = EEo - EEi

  if (DEE < 0) then
     fspin = -1
  else
     if ( exp(-DEE) > rnd ) then
        fspin=-1
     else
        fspin=1
     end if
  end if
  
end function flipspin
