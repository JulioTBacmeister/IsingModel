module renorm_utils

  PUBLIC zzNxNp
  PUBLIC blockrg
  
contains

!==================================================

  subroutine zzNxNp( NX, x0, x1,x2,x3, bltz,ccc,lclosedbc)
  implicit none

  real, intent(inout)  :: x1,x2,x3,x0
  integer, intent(in)  :: NX
  integer, intent(out) :: CCC(0:Nx+1,0:Nx+1,2**(Nx*Nx) )
  real, intent(out)    :: bltz( 2**(Nx*Nx) )
  logical, intent(in),optional :: lclosedbc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  logical :: llclosedbc
  integer :: a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,s
  integer :: i1,i2,i3,i4,Nspin,Nconf, ii,jj,kk
  integer :: NN,NNN,NNNN
  real    :: Econfig, boltzmann_wgt,w1,w2,w3,ZZZ,ZZp,mm1,mm2
  integer :: CC( Nx*Nx,2**(Nx*Nx) )
  
  ! Input coupling coefficients are assumed to be in the
  ! form K = J/(kb*T) where J's come from Hamiltonian
  !
  !       H = SUM_(ij=NN) J*S_i*S_j
  !
  ! Critical value: K=0.4407... (analytical, Onsager)
  !

  if (present(lclosedbc)) then
     llclosedbc = lclosedbc
  else
     llclosedbc = .false.
  end if
  Nspin=NX*NX
  Nconf=2**Nspin

  CCC(:, :, :)=0.
  do i=1,Nspin
  do c=1,Nconf
     L=2**(i-1)
     s=(-1)**INT( (c-1)/L )
     cc( i, c ) = s
  end do
  end do

  CCC(1:NX,1:NX,1:Nconf) = RESHAPE( CC, (/ NX,NX,Nconf /) )

  ZZZ = 0.
  if(llclosedbc) then
     do kk=1,Nconf
       ! == Closed BCs ===
       ccc( 1:NX , NX+1, kk)    = 0
       ccc( 1:NX , 0   , kk)    = 0
       ccc( NX+1 , 1:NX, kk)    = 0 
       ccc( 0    , 1:NX, kk)    = 0 
       ccc( NX+1 , Nx+1, kk)    = 0 
       ccc(  0   , Nx+1, kk)    = 0
       ccc(  Nx+1,  0  , kk)    = 0
       ccc(   0  ,  0  , kk)    = 0 
    end do
  else
     do kk=1,Nconf
       ! == Periodic BCs ===
       ccc( 1:NX , NX+1, kk)    = ccc( 1:NX  ,  1 , kk)
       ccc( 1:NX , 0   , kk)    = ccc( 1:NX  ,  NX, kk)
       ccc( NX+1 , 1:NX, kk)    = ccc(  1   , 1:NX, kk)
       ccc( 0    , 1:NX, kk)    = ccc(  NX   , 1:NX, kk)
       ccc( NX+1 , Nx+1, kk)    = ccc(  1   ,   1, kk)
       ccc(  0   , Nx+1, kk)    = ccc(  NX  ,   1, kk)
       ccc(  Nx+1,  0  , kk)    = ccc(  1   ,  NX, kk)
       ccc(   0  ,  0  , kk)    = ccc( NX   ,  NX, kk)
    end do
 end if

  do kk=1,Nconf
     NN   = 0.
     NNN  = 0.
     NNNN = 0.
     do j=1,NX
     do i=1,NX
        a = ccc(i,j,kk)
        b = ccc(i+1,j,kk)
        c = ccc(i,j+1,kk)
        d = ccc(i+1,j+1,kk)

        NN   = NN    + a*b + a*c
        NNN  = NNN   + a*d + b*c
        NNNN = NNNN  + a*b*c*d
     end do
     end do
     Econfig = x1*NN + x2*NNN + x3*NNNN + x0
     bltz(kk) =  exp( Econfig )
  end do

end subroutine zzNxNp
!=============================================

function blockrg(ccc,nx,b,lclosedbc) result (bccc)
  implicit none
  integer, intent(IN) :: ccc(0:nx+1,0:nx+1),nx,b
  logical, intent(in),optional :: lclosedbc
  integer :: bccc(0:nx/b+1,0:nx/b+1)

  integer :: nxb,i,j,spinsum,ii,jj
  logical :: llclosedbc

  if (present(lclosedbc)) then
     llclosedbc = lclosedbc
  else
     llclosedbc = .false.
  end if

  bccc=0.
  nxb=nx/b
  do jj=1,nxb
  do ii=1,nxb
     spinsum=0
     do j=(jj-1)*b+1,jj*b
        do i=(ii-1)*b+1,ii*b
           spinsum = spinsum + ccc(i,j)
        end do
     end do
     if ( spinsum /= 0 ) then
        bccc(ii,jj) = spinsum / abs( spinsum )
     else
        bccc(ii,jj) = ccc( (ii-1)*b+1 , b*jj )
     end if
  end do
  end do

  if ( .not.llclosedbc) then
    ! == Periodic BCs ===
       bccc( 1:NXB , NXB+1 )    = bccc( 1:NXB  ,  1  )
       bccc( 1:NXB , 0    )    = bccc( 1:NXB  ,  NXB )
       bccc( NXB+1 , 1:NXB )    = bccc(  1   , 1:NXB )
       bccc( 0    , 1:NXB )    = bccc(  NXB   , 1:NXB )
       bccc( NXB+1 , Nxb+1 )    = bccc(  1   ,   1 )
       bccc(  0   , Nxb+1 )    = bccc(  NXB  ,   1 )
       bccc(  Nxb+1,  0   )    = bccc(  1   ,  NXB )
       bccc(   0  ,  0   )    = bccc( NXB   ,  NXB )
    else
    ! == Closed BCs ===
       bccc( 1:NXB , NXB+1 )    = 0 
       bccc( 1:NXB , 0    )    = 0 
       bccc( NXB+1 , 1:NXB )    = 0 
       bccc( 0    , 1:NXB )    = 0
       bccc( NXB+1 , Nxb+1 )    = 0
       bccc(  0   , Nxb+1 )    = 0 
       bccc(  Nxb+1,  0   )    = 0 
       bccc(   0  ,  0   )    = 0 
    endif
  
end function blockrg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!==================================================
  function s4x3ftn(ccc, NX ) result (snm)
  implicit none
  integer, intent(IN) :: NX
  integer, intent(IN) :: ccc(0:NX+1,0:NX+1)
  real :: snm(4,3)
  integer :: i,i2,j,j2,ir,NN,NNN,NNNN,a,b,c,d,N11,N12,N13,N22,N23,N33
  real :: ss, r, smean, svar
  
  ! Go through spin configuration CCC and add up
  ! interaction terms:
  !       nearest neighbors - NN
  !       next nearest ..   - NNN
  !       square/quad       - NNNN
  NN   = 0.
  NNN  = 0.
  NNNN = 0.
  do j=1,NX
  do i=1,NX
        a = ccc(i,j)
        b = ccc(i+1,j)
        c = ccc(i,j+1)
        d = ccc(i+1,j+1)

        NN   = NN    + a*b + a*c
        NNN  = NNN   + a*d + b*c
        NNNN = NNNN  + a*b*c+d
   end do
   end do
   
   ! construct cross-correlations <S_n S_m>
   ! for CCC
       snm(1,1) = dble(NN * NN) 
       snm(2,1) = dble(NNN * NN)
       snm(3,1) = dble(NNNN * NN)

       snm(1,2) = dble(NN * NNN)
       snm(2,2) = dble(NNN * NNN)
       snm(3,2) = dble(NNNN * NNN)

       snm(1,3) = dble(NN * NNNN)
       snm(2,3) = dble(NNN * NNNN)
       snm(3,3) = dble(NNNN * NNNN)

   ! Append expectation values for correlations
       snm(4,1) = dble( NN )    
       snm(4,2) = dble( NNN )    
       snm(4,3) = dble( NNNN )

    ! Seems intuitive to rescale by NX**2 (but should we???)
       snm = snm /( NX**2 )
       
 end function s4x3ftn
!==================================================

end module renorm_utils
