subroutine ising_MC_rlzn( NC, NX, x1,x2,x3  ) ! , ZZZc, CCC  )

  use diags, only : energy_ccc,corrftn
  
  implicit none

  real, intent(inout)  :: x1,x2,x3
  integer, intent(in) :: NX,NC
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,s
  integer :: i1,i2,i3,i4,Nspin,Nconf, ii,jj,kk
  real    :: NN,NNN,NNNN
  real    :: Econfig, boltzmann_wgt,w1,w2,w3,ZZZ,ZZp,mm1,mm2,DEE,zzz2,mm1x
  integer :: flipspin,ncorr
  integer :: CCC(0:Nx+1,0:Nx+1)
  real    :: Delta_Ex( NC ), Delta_Ep( NC ), Energy(NC)
  real    :: MeanMag(NC),ran1,ran2,ran3
  real    :: corr(NX,NX,NX,NX),corrv( 0:NX/2 ),corr3v(3)
  integer, allocatable :: seed(:)
  integer :: clock,countv(0:NX/2)
  real    :: Ranu( NX,NX ),mc1(NC),mc2(NC),mc3(NC),DEEc,DEEx,E_init,E_final
  integer :: IX(NC), JY(NC), iflip, jflip

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

  E_init = energy_ccc(ccc,x1,x2,x3,NX)

  mm1=0.
  DEE=0.


  ! header for output files
  write(611) NC,NX,X1,1./X1
  write(711) NC,NX,X1,1./X1
  write(811) NC,NX,X1,1./X1

  do kk=1,NC-1
       call random_number( ran1 )
       jy(kk) = INT( ran1*(NX))+1
       call random_number( ran2 )
       ix(kk) = INT( ran2*(NX))+1

       
     ! "propose" a flipped spin
     ! at iflip,jflip
     !========================
     iflip = IX(kk)
     jflip = JY(kk)

     !===================================
     ! Calculate energy change that would
     ! ensue from proposed spin flip
     !====================================
     !   Nearest-neighbor (NN) contribution
     DEE = 2 * X1 *  CCC( iflip,   jflip  ) * &
                  (  CCC( iflip+1, jflip  )   &
                  +  CCC( iflip-1, jflip  )   &
                  +  CCC( iflip  , jflip+1)   &
                  +  CCC( iflip  , jflip-1)  )
#if 0
     !   Next nearest neigh. (NNN) contribution
     DEE = DEE + 2 * X2 * CCC( iflip, jflip ) * &
                  (  CCC( iflip+1, jflip+1  )   &
                  +  CCC( iflip-1, jflip+1  )   &
                  +  CCC( iflip+1, jflip-1)   &
                  +  CCC( iflip-1, jflip-1)  )
     !   QUAD/NNNN contribution
     DEE = DEE + 2 * X3 * CCC( iflip, jflip ) * &
                  (  CCC( iflip+1, jflip+1  )*CCC( iflip,jflip+1  )*CCC( iflip+1,jflip )   &    !ne
                  +  CCC( iflip  , jflip+1  )*CCC( iflip-1,jflip+1)*CCC( iflip-1,jflip )   &    !nw
                  +  CCC( iflip+1, jflip    )*CCC( iflip-1,jflip-1)*CCC( iflip,jflip-1 )   &    !sw
                  +  CCC( iflip  , jflip-1  )*CCC( iflip+1,jflip-1)*CCC( iflip+1,jflip ) )      !se
#endif

     !==============================
     ! pick a uniformly dist. random
     ! number between [0,1]
     !===============================
     call random_number( ran3 )
     mc2(kk) = ran3

     !===========================================
     ! If energy change <0 or if Boltzmann
     ! factor change is bigger than ran3 then
     ! keep proposed flip. (Note to stupid: Since
     ! ran3 must be <1 these conditions are redundant.)
     !============================================
     if ( ( DEE<0.) .OR. ( mc2(KK) < exp(-DEE) ) ) then
        ! KEEP -----
        CCC( iflip,   jflip  ) = -CCC( iflip,   jflip  )
        DEEx=DEE
     else
        ! REJECT ---
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


       Energy(kk) = energy_ccc(ccc,x1,x2,x3,NX)
       mm1x = abs(sum( CCC(1:NX,1:NX) ) )
       mm1  = mm1+ mm1x !abs(sum( CCC(1:NX,1:NX) ) )
       
       Delta_Ex(kk) = DEEx
       Delta_Ep(kk) = DEE
       MeanMag(kk)  = REAL(mm1x)/(NX*NX)


       if ( (mod(kk-1,2*(10**4))==0).AND.(NX>10)) then
          write(*,*) kk,MeanMag(kk)
          write(611) ccc,kk
       endif
       if ( (mod(kk-1,80000*(10**4))==0).AND.(NX>10)) then
          write(*,*) "  --  Calculating correlation functions at ",kk
          corrv = corrftn( ccc, NX,counts=countv)
          write(711) corrv,countv,kk
       end if
  end do

  
  mm1 = mm1/ (NX*NX) /NC
  write(*,*) " Mean value of ABS magnetization scheme 2 :  ",mm1

  E_final = energy_ccc(ccc,x1,x2,x3,NX)

  write(511) NC,NX,X1,1./X1
  write(511) ranu
  write(511) mc1,mc2,mc3,IX,JY
  write(511) E_init,E_final
  write(511) Delta_Ex,Delta_Ep,MeanMag,Energy
  write(511) corrv

  !write(*,*)  " : ",mm1


  !=================================
  ! You need to learn what the
  ! expectation value:
  !
  !          Sum_states [ m* exp(-BE) ]
  !    <m> = --------------------------
  !           Sum_states [ exp(-BE) ]
  !
  ! actually means
  

  corr( :,:,:,: )= 0.

#if 0  
  do kk=1,Nconf
  do j=1,NX
     do i=1,NX
        do jj=1,NX
           do ii=1,NX
              corr(i,j,ii,jj) = corr(i,j,ii,jj) + &
                   ZZZc(kk)* (CCC(i,j,kk)*CCC(ii,jj,kk))/ZZZ
           end do
        end do
     end do
  end do
  end do
#endif

  
end subroutine ising_MC_rlzn
