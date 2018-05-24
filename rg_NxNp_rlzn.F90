subroutine rg_NxNp_rlzn( NX, BLK, x1,x2,x3, SS,SSB )
  use renorm_utils, only : blockrg 
  implicit none

  real, intent(inout)  :: x1,x2,x3
  integer, intent(in) :: NX,BLK
  real, intent(out) :: SS(4,3),SSB(4,3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real    :: z(-1:1,-1:1,-1:1,-1:1)
  integer :: a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,s
  integer :: i1,i2,i3,i4,Nspin,Nconf, ii,jj,kk
  integer :: NN,NNN,NNNN
  real    :: Econfig, boltzmann_wgt,w1,w2,w3,ZZZ,ZZp,mm1,mm2
  integer :: iblock
  integer :: CCC(0:Nx+1,0:Nx+1,2**(Nx*Nx) ),CC(Nx*Nx,2**(Nx*Nx) )
  integer :: bccc( 0:NX/2+1, 0:Nx/2+1, 2**(NX*NX) )
  real    :: ZZZc( 2**(Nx*Nx) ) , EEc( 2**(Nx*Nx) )
  real    :: corr(NX,NX,NX,NX),tps
  real    :: S1,S2,S3, S11,S12,S13,S21,S22,S23,S31,S32,S33
  
  w1=exp(x1)
  w2=exp(x2)
  w3=exp(x3)

  bccc(:,:,:)=0.

  S1=0.
  S2=0.
  S3=0.
  S11=0.
  S12=0.
  S13=0.
  S21=0.
  S22=0.
  S23=0.
  S31=0.
  S32=0.
  S33=0.

  ! Input coupling coefficients are assumed to be in the
  ! form K = J/(kb*T) where J's come from Hamiltonian
  !
  !       H = SUM_(ij=NN) J*S_i*S_j
  !
  ! Critical value: K=0.4407... (analytical, Onsager)
  !

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
  ! 2x2
  !--------- 
  ! c d
  ! a b 
  !---------
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

     Econfig = x1*NN + x2*NNN + x3*NNNN
     boltzmann_wgt =  exp( Econfig )


     ! Now calculate even <spin-spin> terms for config kk
     S1 = S1 + boltzmann_wgt*dble( NN )
     S2 = S2 + boltzmann_wgt*dble( NNN )
     S3 = S3 + boltzmann_wgt*dble( NNNN )

     S11 = S11 + boltzmann_wgt*dble( NN * NN )
     S21 = S21 + boltzmann_wgt*dble( NNN * NN )
     S31 = S31 + boltzmann_wgt*dble( NNNN * NN )

     S12 = S12 + boltzmann_wgt*dble( NN * NNN )
     S22 = S22 + boltzmann_wgt*dble( NNN * NNN )
     S32 = S32 + boltzmann_wgt*dble( NNNN * NNN )

     S13 = S13 + boltzmann_wgt*dble( NN * NNNN )
     S23 = S23 + boltzmann_wgt*dble( NNN * NNNN )
     S33 = S33 + boltzmann_wgt*dble( NNNN * NNNN )

     ZZZ = ZZZ + boltzmann_wgt
     ZZZc(kk)=boltzmann_wgt
     EEc(kk) = Econfig
  end do

  S1=S1/ZZZ
  S2=S2/ZZZ
  S3=S3/ZZZ
  SS(1,1) = S11 / ZZZ - S1*S1
  SS(1,2) = S12 / ZZZ - S1*S2
  SS(1,3) = S13 / ZZZ - S1*S3
  SS(2,1) = S21 / ZZZ - S2*S1
  SS(2,2) = S22 / ZZZ - S2*S2
  SS(2,3) = S23 / ZZZ - S2*S3
  SS(3,1) = S31 / ZZZ - S3*S1
  SS(3,2) = S32 / ZZZ - S3*S2
  SS(3,3) = S33 / ZZZ - S3*S3

  SS(4,1) = S1
  SS(4,2) = S2
  SS(4,3) = S3


if (BLK>=2) then
  S1=0.
  S2=0.
  S3=0.
  S11=0.
  S12=0.
  S13=0.
  S21=0.
  S22=0.
  S23=0.
  S31=0.
  S32=0.
  S33=0.
  do kk=1,Nconf
     bccc(:,:,kk) = blockrg( ccc(:,:,kk) , NX , BLK )
     NN   = 0.
     NNN  = 0.
     NNNN = 0.
     do j=1,NX/BLK
     do i=1,NX/BLK
        a = bccc(i,j,kk)
        b = bccc(i+1,j,kk)
        c = bccc(i,j+1,kk)
        d = bccc(i+1,j+1,kk)

        NN   = NN    + a*b + a*c
        NNN  = NNN   + a*d + b*c
        NNNN = NNNN  + a*b*c*d
     end do
     end do

     !Econfig = x1*NN + x2*NNN + x3*NNNN
     !boltzmann_wgt =  exp( Econfig )
     boltzmann_wgt=ZZZc(kk)


     ! Now calculate even <spin-spin> terms for config kk
     S1 = S1 + boltzmann_wgt*dble( NN )
     S2 = S2 + boltzmann_wgt*dble( NNN )
     S3 = S3 + boltzmann_wgt*dble( NNNN )

     S11 = S11 + boltzmann_wgt*dble( NN * NN )
     S21 = S21 + boltzmann_wgt*dble( NNN * NN )
     S31 = S31 + boltzmann_wgt*dble( NNNN * NN )

     S12 = S12 + boltzmann_wgt*dble( NN * NNN )
     S22 = S22 + boltzmann_wgt*dble( NNN * NNN )
     S32 = S32 + boltzmann_wgt*dble( NNNN * NNN )

     S13 = S13 + boltzmann_wgt*dble( NN * NNNN )
     S23 = S23 + boltzmann_wgt*dble( NNN * NNNN )
     S33 = S33 + boltzmann_wgt*dble( NNNN * NNNN )
  end do

  S1=S1/ZZZ
  S2=S2/ZZZ
  S3=S3/ZZZ
  SSB(1,1) = S11 / ZZZ - S1*S1
  SSB(1,2) = S12 / ZZZ - S1*S2
  SSB(1,3) = S13 / ZZZ - S1*S3
  SSB(2,1) = S21 / ZZZ - S2*S1
  SSB(2,2) = S22 / ZZZ - S2*S2
  SSB(2,3) = S23 / ZZZ - S2*S3
  SSB(3,1) = S31 / ZZZ - S3*S1
  SSB(3,2) = S32 / ZZZ - S3*S2
  SSB(3,3) = S33 / ZZZ - S3*S3

  SSB(4,1) = S1
  SSB(4,2) = S2
  SSB(4,3) = S3
else
  SSB(:,:)=-999.
end if

end subroutine rg_NxNp_rlzn
