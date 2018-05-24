subroutine rg_MC_rlzn( NC, NX, B, x1,x2,x3 , S0, S1, lrenorm  ) ! , ZZZc, CCC  )

  use renorm_utils, only : blockrg , s4x3ftn
  
  implicit none

  real,    intent(in)  :: x1,x2,x3
  integer, intent(in)  :: NX,NC,B
  real,    intent(out) :: S0(4,3),S1(4,3)
  logical, optional, intent(in) :: lrenorm
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: a,c,d,e,f,g,h,i,j,k,l,m,n,o,p,s
  integer :: i1,i2,i3,i4,Nspin,Nconf, ii,jj,kk
  real    :: NN,NNN,NNNN
  real    :: Econfig, boltzmann_wgt,w1,w2,w3,ZZZ,ZZp,mm1,mm2,DEE,zzz2,mm1x
  integer :: flipspin,ncorr,isweep
  integer :: CCC(0:Nx+1,0:Nx+1), Bccc( 0:nx/b+1 , 0:nx/b+1 )
  real    :: Delta_Ex( NC ), Delta_Ep( NC ), Energy(NC)
  real    :: MeanMag(NC),ran1,ran2,ran3
  real    :: corr(NX,NX,NX,NX),corrv( 0:NX/2 ),corr3v(3)
  integer, allocatable :: seed(:)
  integer :: clock,countv(0:NX/2)
  real    :: Ranu( NX,NX ),mc1(NC),mc2(NC),mc3(NC),DEEc,DEEx,E_init,E_final
  integer :: IX(NC), JY(NC), iflip, jflip,nxb
  real    :: snm0(4,3),snm1(4,3)
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
  
  if (present(lrenorm)) then
     llrenorm=lrenorm
  else
     llrenorm = .TRUE.
  end if
 
  isweep=0
  S0(:,:)=0.
  S1(:,:)=0.
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

    !write(*,*) " shape bccc ",shape(bccc)
    write(811) NC,NX,B,NXB,X1,1./X1

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
          
    ! == Periodic BCs ===
       ccc( 1:NX , NX+1 )    = ccc( 1:NX  ,  1  )
       ccc( 1:NX , 0    )    = ccc( 1:NX  ,  NX )
       ccc( NX+1 , 1:NX )    = ccc(  1   , 1:NX )
       ccc( 0    , 1:NX )    = ccc(  NX   , 1:NX )
       ccc( NX+1 , Nx+1 )    = ccc(  1   ,   1 )
       ccc(  0   , Nx+1 )    = ccc(  NX  ,   1 )
       ccc(  Nx+1,  0   )    = ccc(  1   ,  NX )
       ccc(   0  ,  0   )    = ccc( NX   ,  NX )

    !====================================
    ! Energy of initial spin configuration
    !=====================================   


  mm1=0.
  DEE=0.


  do kk=1,NC-1
       call random_number( ran1 )
       jy(kk) = INT( ran1*(NX))+1
       call random_number( ran2 )
       ix(kk) = INT( ran2*(NX))+1

       
     ! Now flip spins
     !===============
     iflip = IX(kk)
     jflip = JY(kk)

     DEE = 2 * X1 *  CCC( iflip,   jflip  ) * &
                  (  CCC( iflip+1, jflip  )   &
                  +  CCC( iflip-1, jflip  )   &
                  +  CCC( iflip  , jflip+1)   &
                  +  CCC( iflip  , jflip-1)  )
#if 1
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

    ! == Periodic BCs ===
       ccc( 1:NX , NX+1)    = ccc( 1:NX  ,  1  )
       ccc( 1:NX , 0   )    = ccc( 1:NX  ,  NX )
       ccc( NX+1 , 1:NX)    = ccc(  1   , 1:NX )
       ccc( 0    , 1:NX)    = ccc(  NX   , 1:NX )
       ccc( NX+1 , Nx+1)    = ccc(  1   ,   1 )
       ccc(  0   , Nx+1)    = ccc(  NX  ,   1 )
       ccc(  Nx+1,  0  )    = ccc(  1   ,  NX )
       ccc(   0  ,  0  )    = ccc( NX   ,  NX )



       snm0 = s4x3ftn( ccc  , nx )
       if (llrenorm) then 
          bccc = blockrg( ccc, nx, b )
          snm1 = s4x3ftn( bccc , nxb )
       else
          bccc = -999
          snm1 = -999
       end if
       
       if (KK >= NC/3 ) then
          S0 = S0 + snm0
          S1 = S1 + snm1
          isweep= isweep+1
       endif

       
       if ( (mod(kk-1,(10**3))==0) ) then
          !write(*,*) kk
          write(811) ccc,bccc,kk,snm0,snm1
       !write(*,*) ccc(1:b,1:b)
       !write(*,*) bccc(1,1)
       !pause
       endif

  end do

  S0 = S0/isweep
  S0(1,1) = S0(1,1)-S0(4,1)**2
  S0(1,2) = S0(1,2)-S0(4,1)*S0(4,2)
  S0(1,3) = S0(1,3)-S0(4,1)*S0(4,3)
  S0(2,2) = S0(2,2)-S0(4,2)**2
  S0(2,3) = S0(2,3)-S0(4,2)*S0(4,3)
  S0(3,3) = S0(3,3)-S0(4,3)**2
  S0(2,1) = S0(1,2)
  S0(3,1) = S0(1,3)
  S0(3,2) = S0(2,3)
  
  S1 = S1/isweep
  S1(1,1) = S1(1,1)-S1(4,1)**2
  S1(1,2) = S1(1,2)-S1(4,1)*S1(4,2)
  S1(1,3) = S1(1,3)-S1(4,1)*S1(4,3)
  S1(2,2) = S1(2,2)-S1(4,2)**2
  S1(2,3) = S1(2,3)-S1(4,2)*S1(4,3)
  S1(3,3) = S1(3,3)-S1(4,3)**2
  S1(2,1) = S1(1,2)
  S1(3,1) = S1(1,3)
  S1(3,2) = S1(2,3)

  
end subroutine rg_MC_rlzn
