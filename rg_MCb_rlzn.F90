subroutine rg_MCb_rlzn( NC, NX, B, x1,x2,x3 , SS )

  use renorm_utils, only : blockrg , s4x3ftn
  
  implicit none

  real,    intent(in)  :: x1,x2,x3
  integer, intent(in)  :: NX,NC,B
  real,    intent(out) :: SS(4,3)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: a,c,d,e,f,g,h,i,j,k,l,m,n,o,p,s
  integer :: aa,bb,cc,dd
  integer :: i1,i2,i3,i4,Nspin,Nconf, ii,jj,kk
  integer :: NN,NNN,NNNN
  real    :: Econfig, boltzmann_wgt,w1,w2,w3,ZZZ,ZZp,mm1,mm2,DEE,zzz2,mm1x
  integer :: flipspin,ncorr,isweep
  integer :: CCC(0:Nx+1,0:Nx+1), Bccc( 0:nx/b+1 , 0:nx/b+1 )
  real    :: Delta_Ex( NC ), Delta_Ep( NC ), Energy(NC)
  real    :: MeanMag(NC),ran1,ran2,ran3
  real    :: corr(NX,NX,NX,NX),corrv( 0:NX/2 ),corr3v(3)
  integer, allocatable :: seed(:)
  integer :: clock,countv(0:NX/2)
  real    :: Ranu( NX,NX ),mc1(NC),mc2(NC),mc3(NC),DEEc,DEEx,E_init,E_final
  integer :: IX(NC), JY(NC), iflip, jflip,nxb, ip, jp
  real    :: S1,S2,S3,S11,S12,S13,S21,S22,S23,S31,S32,S33,DEE1,DEE2,DEE3
  logical :: llrenorm

  !! FUNCTIONS
  !!real    :: energy_ccc
  
  !======================================================
  ! Input coupling coefficients are assumed to be in the
  ! form K = -J/(kb*T) where J's come from Hamiltonian
  !
  !       H = SUM_(ij=NN) -J*S_i*S_j
  !
  ! Critical value: K=0.4407... (analytical, Onsager).
  ! Note negative sign in front of J.  Positive J then  
  ! means *lower* Energy with aligned spins.
  !======================================================
  isweep=0

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
  
  nxb=nx/b
  call random_seed(size = n)
  allocate(seed(n))
  call random_seed(get=seed)

    !write(*,*) "T per site=",1./X1," K1=",X1

! SYSTEM_CLOCK call results in a different sequence each time:

  
    CALL SYSTEM_CLOCK(COUNT=clock)     
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put=seed)
    CALL RANDOM_NUMBER(ranu)
    mc1=0.
    mc3=0.

    ccc=0.
    
    !=============================
    ! Initial Spin configuration
    !=============================
#if 0
    !random spins
    ccc(:,:)=0
    where(ranu >=0.5)
       ccc(1:NX,1:NX)=1
    elsewhere
       ccc(1:NX,1:NX)=-1
    end where
#endif
#if 0
       ! cold start
       ccc(:,:)=1
#endif
#if 0
       ! Checker board
       ccc(:,:)=1
       do j=1,NX,2
          do i=1,NX,2
             ccc(i,j)=-1
          end do
       end do
#endif             
#if 0
       ! "Continent"
       ccc(:,:)=1
       do j=1,NX
          do i=1,NX/2
             ccc(i,j)=-1
          end do
       end do
#endif
#if 1
       ! Stripes
       ccc(:,:)=1
       do i=1,NX-4,4
          ccc(i:i+2,:)=-1
       end do
#endif       
          
#if 1
       ! == Periodic BCs ===
       ccc( 1:NX , NX+1 )    = ccc( 1:NX  ,  1  )
       ccc( 1:NX , 0    )    = ccc( 1:NX  ,  NX )
       ccc( NX+1 , 1:NX )    = ccc(  1   , 1:NX )
       ccc( 0    , 1:NX )    = ccc(  NX   , 1:NX )
       ccc( NX+1 , Nx+1 )    = ccc(  1   ,   1 )
       ccc(  0   , Nx+1 )    = ccc(  NX  ,   1 )
       ccc(  Nx+1,  0   )    = ccc(  1   ,  NX )
       ccc(   0  ,  0   )    = ccc( NX   ,  NX )
#endif
    !====================================
    ! Energy of initial spin configuration
    !=====================================   


  mm1=0.
  DEE=0.

  do kk=1,NC-1
#if 1
       call random_number( ran1 )
       jy(kk) = INT( ran1*(NX))+1
       call random_number( ran2 )
       ix(kk) = INT( ran2*(NX))+1

       
     ! Now flip spins
     !===============
     iflip = IX(kk)
     jflip = JY(kk)
#endif
#if 0
     ysweep: do jflip=1,NX
     xsweep: do iflip=1,NX
#endif     
#if 1
     ip=iflip
     jp=jflip
     DEE1 = 2 * X1 *  CCC( ip,   jp  ) * ( CCC( ip+1,jp )+CCC(ip,jp+1)+CCC(ip-1,jp)+CCC(ip,jp-1) ) 

     ! NNN contribution
     DEE2 = 2 * X2 * CCC( ip, jp ) * ( CCC(ip+1,jp+1) + CCC(ip+1,jp-1)+CCC(ip-1,jp+1)+CCC(ip-1,jp-1) ) 

          ! QUAD/NNNN contribution
     DEE3 = 2 * X3 * CCC( ip, jp ) * &    
                  (  CCC( ip-1, jp  )*CCC( ip-1,jp+1  )*CCC( ip,jp+1 )  & ! |-
                  +  CCC( ip, jp+1  )*CCC( ip+1,jp+1  )*CCC( ip+1,jp )  & ! -|
                  +  CCC( ip+1,jp   )*CCC( ip+1,jp-1  )*CCC( ip,jp-1 )  & ! _|
                  +  CCC( ip,jp-1   )*CCC( ip-1,jp-1  )*CCC( ip-1,jp ) )  ! |_
     DEE  = DEE1 + DEE2 + DEE3 
#endif
#if 0
     DEE = 2 * X1 *  CCC( iflip,   jflip  ) * &
                  (  CCC( iflip+1, jflip  )   &
                  +  CCC( iflip-1, jflip  )   &
                  +  CCC( iflip  , jflip+1)   &
                  +  CCC( iflip  , jflip-1)  )

     ! NNN contribution
     DEE = DEE + 2 * X2 * CCC( iflip, jflip ) * &
                  (  CCC( iflip+1, jflip+1  )   &
                  +  CCC( iflip-1, jflip+1  )   &
                  +  CCC( iflip+1, jflip-1)   &
                  +  CCC( iflip-1, jflip-1)  )
     ! QUAD/NNNN contribution
     DEE = DEE + 2 * X3 * CCC( iflip, jflip ) * &
                  (  CCC( iflip+1, jflip+1  )*CCC( iflip,jflip+1  )*CCC( iflip+1,jflip )   &    !ne
                  +  CCC( iflip  , jflip+1  )*CCC( iflip-1,jflip+1)*CCC( iflip-1,jflip )   &    !nw
                  +  CCC( iflip-1, jflip    )*CCC( iflip-1,jflip-1)*CCC( iflip,jflip-1 )   &    !sw
                  +  CCC( iflip  , jflip-1  )*CCC( iflip+1,jflip-1)*CCC( iflip+1,jflip ) )      !se
#endif

     call random_number( ran3 )
     
     if ( ( DEE<0.) .OR. ( ran3 < exp(-DEE) ) ) then
        CCC( iflip,   jflip  ) = -CCC( iflip,   jflip  )
        DEEx=DEE
     else
        DEEx=0.
     end if


  !end do xsweep
  !end do ysweep

#if 1
    ! == Periodic BCs ===
       ccc( 1:NX , NX+1)    = ccc( 1:NX  ,  1  )
       ccc( 1:NX , 0   )    = ccc( 1:NX  ,  NX )
       ccc( NX+1 , 1:NX)    = ccc(  1   , 1:NX )
       ccc( 0    , 1:NX)    = ccc(  NX   , 1:NX )
       ccc( NX+1 , Nx+1)    = ccc(  1   ,   1 )
       ccc(  0   , Nx+1)    = ccc(  NX  ,   1 )
       ccc(  Nx+1,  0  )    = ccc(  1   ,  NX )
       ccc(   0  ,  0  )    = ccc( NX   ,  NX )
#endif

       if (b >=2 ) then
          bccc = blockrg( ccc, nx, b )
       else
          bccc(:,:)=ccc(:,:)
       end if
         !if ( mod(kk,100*nx*nx)==0 ) write(2011) ccc,bccc,kk
       
       if (KK >= NC/3 ) then
               ! Now calculate even <spin-spin> terms for config kk
          NN   = 0.
          NNN  = 0.
          NNNN = 0.
          do j=1,NXb
          do i=1,NXb
             aa = bccc(i,j)
             bb = bccc(i+1,j)
             cc = bccc(i,j+1)
             dd = bccc(i+1,j+1)

             NN   = NN    + aa*bb + aa*cc
             NNN  = NNN   + aa*dd + bb*cc
             NNNN = NNNN  + aa*bb*cc*dd
          end do
          end do

          S1 = S1 +  dble( NN )
          S2 = S2 +  dble( NNN )
          S3 = S3 +  dble( NNNN )

          S11 = S11 +  dble( NN * NN )
          S21 = S21 +  dble( NNN * NN )
          S31 = S31 +  dble( NNNN * NN )

          S12 = S12 +  dble( NN * NNN )
          S22 = S22 +  dble( NNN * NNN )
          S32 = S32 +  dble( NNNN * NNN )

          S13 = S13 +  dble( NN * NNNN )
          S23 = S23 +  dble( NNN * NNNN )
          S33 = S33 +  dble( NNNN * NNNN )
          isweep= isweep+1
          ZZZ=dble( isweep )
       endif

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


  
end subroutine rg_MCb_rlzn
