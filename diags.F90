module diags

  PUBLIC energy_ccc
  PUBLIC iblock
  PUBLIC flipspin
  PUBLIC corrftn
  
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function corrftn(ccc, NX,counts) result (corrs)
  implicit none
  integer, intent(IN) :: NX
  integer, intent(IN) :: ccc(0:NX+1,0:NX+1)
  integer, intent(OUT),optional :: counts(0:NX/2)
  real :: corrs(0:NX/2)
  integer :: i,i2,j,j2,ir
  real :: ss, r, smean, svar


   !  mean spin
  smean = sum( ccc(1:NX,1:NX ) ) /dble(NX*NX)


  
  corrs(:)=0.
  counts(:)=0
  do j2=1,NX
  do i2=1,NX
     do j=1,NX
     do i=1,NX
        r=SQRT( dble(i-i2)**2 + dble(j-j2)**2 )
        ir=INT( r )
        if (ir <= NX/2 ) then
           ss         = dble( (ccc(i,j)-smean)*(ccc(i2,j2)-smean) )
           counts(ir) = counts(ir) + 1
           corrs(ir)  = corrs(ir) + dble(ss)
        end if
     end do
     end do
   end do
   end do

   where( counts > 0)
      corrs = corrs/counts
   elsewhere
      corrs = 0.
   end where

   ! divide by variance
   svar  = corrs(0)
   corrs = corrs/svar

   
 end function corrftn
!==================================================
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function iblock(a,b,c,d) result (blockspin)
  implicit none
  integer, intent(IN) :: a,b,c,d
  integer :: blockspin

  integer :: spinsum

  spinsum  = a + b + c + d

  if (spinsum == 0) then
     blockspin = a / abs(a)
  else
     blockspin = spinsum / abs( spinsum )
  endif
  
end function iblock



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

end module diags
