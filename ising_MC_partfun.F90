subroutine ising_MC_partfun( NC, NX, x1,x2,x3) ! , ZZZc, CCC  )
  implicit none

  real, intent(inout)  :: x1,x2,x3
  integer, intent(in) :: NX,NC
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real    :: z(-1:1,-1:1,-1:1,-1:1)
  integer :: a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,s
  integer :: i1,i2,i3,i4,Nspin,Nconf, ii,jj,kk
  real    :: NN,NNN,NNNN
  real    :: Econfig, boltzmann_wgt,w1,w2,w3,ZZZ,ZZp,mm1,mm2,DEE,zzz2
  integer :: flipspin
  integer :: CCC(0:Nx+1,0:Nx+1,NC ) , CCC_sum(NC) !,CC(Nx*Nx,NC)
  real    :: ZZZc( NC ) , EEc( NC )
  real    :: corr(NX,NX,NX,NX)
  integer, allocatable :: seed(:)
  integer :: clock
  real    :: Ranu( NX,NX ),mc1(NC),mc2(NC),mc3(NC),DEEc(NC),DEEx(NC)
  integer :: IX(NC), JY(NC), iflip, jflip
  integer :: CCs(-1:1,-1:1),flipval
  
  call random_seed(size = n)
  allocate(seed(n))
  call random_seed(get=seed)

    !!write(*,*) "T per site=",1./X1," K1=",X1



! SYSTEM_CLOCK call results in a different sequence each time:

  
    CALL SYSTEM_CLOCK(COUNT=clock)     
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put=seed)
    CALL RANDOM_NUMBER(ranu)

    CALL SYSTEM_CLOCK(COUNT=clock)     
    seed = clock + 197 * (/ (i - 1, i = 1, n) /)
    call random_seed(put=seed)
    CALL RANDOM_NUMBER(mc1)

    CALL SYSTEM_CLOCK(COUNT=clock)     
    seed = clock + 113 * (/ (i - 1, i = 1, n) /)
    call random_seed(put=seed)
    CALL RANDOM_NUMBER(mc2)

    CALL SYSTEM_CLOCK(COUNT=clock)     
    seed = clock + 839 * (/ (i - 1, i = 1, n) /)
    call random_seed(put=seed)
    CALL RANDOM_NUMBER(mc3)

    
    ccc(:,:,:)=0
    where(ranu >=0.5)
       ccc(1:NX,1:NX,1)=1
    elsewhere
       ccc(1:NX,1:NX,1)=-1
    end where

    ! == Periodic BCs ===
       ccc( 1:NX , NX+1, 1)    = ccc( 1:NX  ,  1  ,  1 )
       ccc( 1:NX , 0 ,   1)    = ccc( 1:NX  ,  NX ,  1 )
       ccc( NX+1 , 1:NX, 1)    = ccc(  1   , 1:NX , 1 )
       ccc( 0    , 1:NX, 1)    = ccc(  NX   , 1:NX , 1 )
       ccc( NX+1 , Nx+1, 1)    = ccc(  1   ,   1,    1 )
       ccc(  0   , Nx+1, 1)    = ccc(  NX  ,   1,    1 )
       ccc(  Nx+1,  0  , 1)    = ccc(  1   ,  NX,    1 )
       ccc(   0  ,  0  , 1)    = ccc( NX   ,  NX,    1 )


    
    write(311) NC,NX
    write(311) ranu
    write(311) mc1,mc2,mc3
    

  !do i=1,NC
  !   jy(I) = INT( mc1(i)*(NX-1) )+1
  !   ix(I) = mc1(i)*(NX-1)*(NX-1) - jy(i)*(NX-1) + 1
  !end do
  do i=1,NC
     jy(I) = INT( mc1(i)*(NX))+1
     ix(I) = INT( mc3(i)*(NX))+1
  end do

  do kk=1,NC-1
     ! Now flip spins
     !===============
     CCC(:,:,KK+1) = CCC(:,:,KK)
     iflip = IX(kk)
     jflip = JY(kk)

     DEE = 2 * X1 *  CCC( iflip,   jflip  , KK ) * &
                  (  CCC( iflip+1, jflip  , KK )   &
                  +  CCC( iflip-1, jflip  , KK )   &
                  +  CCC( iflip  , jflip+1, KK )   &
                  +  CCC( iflip  , jflip-1, KK )  )

     if ( ( DEE<0.) .OR. ( mc2(KK) < exp(-DEE) ) ) then
        CCC( iflip,   jflip  , KK+1 ) = -CCC( iflip,   jflip  , KK )
        DEEx(KK)=DEE
     else
        DEEx(KK)=0.
     end if
     !if ( ( DEE<0.) .OR. ( alog(mc2(KK)) < -DEE ) ) CCC( iflip,   jflip  , KK+1 ) = -CCC( iflip,   jflip  , KK )

     DEEc(KK) = DEE
     
    ! == Periodic BCs ===
       ccc( 1:NX , NX+1, KK+1)    = ccc( 1:NX  ,  1  ,  KK+1 )
       ccc( 1:NX , 0 ,   KK+1)    = ccc( 1:NX  ,  NX ,  KK+1 )
       ccc( NX+1 , 1:NX, KK+1)    = ccc(  1   , 1:NX , KK+1 )
       ccc( 0    , 1:NX, KK+1)    = ccc(  NX   , 1:NX , KK+1 )
       ccc( NX+1 , Nx+1, KK+1)    = ccc(  1   ,   1,    KK+1 )
       ccc(  0   , Nx+1, KK+1)    = ccc(  NX  ,   1,    KK+1 )
       ccc(  Nx+1,  0  , KK+1)    = ccc(  1   ,  NX,    KK+1 )
       ccc(   0  ,  0  , KK+1)    = ccc( NX   ,  NX,    KK+1 )

  end do
  DEEc(NC)=0.

     
  write(311) IX,JY,DEEc,DEEx

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
  
  ZZZ = 0.
  
  ! 2x2
  !---------
  ! a b 
  ! c d 
  !---------
  do kk=1,NC
     NN   = 0.
     NNN  = 0.
     NNNN = 0.
     do j=1,NX-1
     do i=1,NX-1
        a = ccc(i,j,kk)
        b = ccc(i+1,j,kk)
        c = ccc(i,j+1,kk)
        d = ccc(i+1,j+1,kk)

        NN   = NN    + a*b + a*c
        NNN  = NNN   + a*d + b*c
        NNNN = NNNN  + a*b*c+d
     end do
   end do

#if 0   
     ! Free Boundaries
     do j=1,NX-1
     i=NX
        a = ccc(i,j,kk)
        c = ccc(i,j+1,kk)
        NN   = NN    + a*c
     end do
     j=NX
     do i=1,NX-1
        a = ccc(i,j,kk)
        b = ccc(i+1,j,kk)
        NN   = NN    + a*b
     end do
#endif
     
#if 1
     ! Periodic Boundaries
     do j=1,NX-1
     i=NX
        a = ccc(i,j,kk)
        b = ccc(1,j,kk)
        c = ccc(i,j+1,kk)
        NN   = NN    + a*c + a*b
     end do
     j=NX
     do i=1,NX-1
        a = ccc(i,j,kk)
        b = ccc(i+1,j,kk)
        c = ccc(i,1,kk)
        NN   = NN    + a*b + a*c
     end do
     j=NX
     i=NX
        a = ccc(i,j,kk)
        b = ccc(1,j,kk)
        c = ccc(i,1,kk)
        NN   = NN    + a*b + a*c
#endif

     Econfig = x1*NN + x2*NNN + x3*NNNN
     !boltzmann_wgt =  w1**nn * w2**nnn * w3**nnnn !   exp( Econfig )
     boltzmann_wgt =  exp( Econfig )

     ZZZ = ZZZ + boltzmann_wgt
     ZZZc(kk)=boltzmann_wgt
     EEc(kk) = -Econfig

#if 0
     ! Now flip spins
     !===============
     CCC(:,:,KK+1) = CCC(:,:,KK)
     iflip = IX(kk)
     jflip = JY(kk)

     CCs(-1:1,-1:1) = CCC( iflip-1:iflip+1, jflip-1:jflip+1, KK )

     flipval = flipspin( CCs, x1,x2,x3, mc2(KK) )

     CCC( iflip, jflip, KK+1 ) = flipval*CCC( iflip, jflip, KK+1 )
#endif
     
  end do


  !write(*,*) "ROutine Ising NxN "
  !write(*,*) "Partition Function ZZZ=",ZZZ


  ! calculate mean magnetization
  mm1=0.
  mm2=0.
  ZZZ2=0.
  do kk=1,NC
     ccc_sum(kk)  =  sum( CCC(1:NX,1:NX,KK) )
     mm2  = mm2+ZZZc(KK) * abs(sum( CCC(1:NX,1:NX,KK) ) )
     mm1  = mm1+abs(sum( CCC(1:NX,1:NX,KK) ) )
  end do

  mm2 = mm2/ (NX*NX) /ZZZ
  mm1 = mm1/ (NX*NX) /NC


  write(311) CCC
  write(311) EEc,ZZZc,ccc_sum
  close(unit=311)



  !=================================
  ! You need to learn what the
  ! expectation value:
  !
  !          Sum_states [ m* exp(-BE) ]
  !    <m> = --------------------------
  !           Sum_states [ exp(-BE) ]
  !
  ! actually means
  
  !write(*,*) " Expectation value of ABS magnetization "
  !write(*,*)  " : ",mm2
  write(*,*) " Mean value of ABS magnetization scheme 1: ",mm1
  !write(*,*)  " : ",mm1

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


  write( 411 ) corr
  
end subroutine ising_MC_partfun
