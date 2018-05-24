program rg_drv
  implicit none
  real    :: x1,x2,x3,y1,y2,y3,T,Tc,dexx
  integer :: niter,i
  !integer, allocatable :: CCC(:,:,:)
  !real, allocatable :: ZZZc(:,:,:)
  integer :: NX,NC,NB,pivot(3), OK,ite,nite,iinr,niinr,noutr,iout
  real :: S0(4,3),S1(4,3),SX(4,3),SXX(4,3),SS0(4,3)
  real :: AA(3,3),delC(3), delX(3),BB(3),KK0(3),KK1(3),delSS,magSS
     
  !write(*,*) " x1,x2,x3 , NC , NX"
  !read(*,*) x1,x2,x3,NC,NX
  !write(*,*) x1,x2,x3 ,NC,NX

  NC = 1*(10**6)  !10**7 !3*(10**6) ! 10**6
  NC = 1500000
  NX = 4
  NB = 2
  x2=0.
  x3=0.
  NC = NC*NX*NX

  nite=10
  niinr = 30 !100
  noutr = 10 !100
  write(911) noutr
  !do i=1,15
  tc=2./( alog(1.+SQRT(2.)))  ! T_c ~ 2.2691853142130221
     !T=  2.5 ! 2.2691853 ! 2.5  ! 3.0 !1.7 ! 2.00 !2.28  !0.5*i
     x1= 0.423 ! 1./T
     T=1/x1
     write(*,*) " T , T_c ",T,Tc
     KK0 = (/ X1, X2, X3 /)
     write(911) kk0
   ! Outer iteration to find fixed point of RG
   do iout=1,noutr   

       call rg_4x4_rlzn( KK0(1),KK0(2),KK0(3) , SXX, S0 ) 
       !!call rg_MCb_rlzn( NC,NX,NB,KK0(1),KK0(2),KK0(3),S0)  !,ZZZc,CCC)
     
       do iinr=1,niinr
       !!write( 6, * ) '                 ---- '

          call rg_NxN_rlzn( NX/NB , KK0(1),KK0(2),KK0(3) , SS0 ) 

          delC(:) = S0(4,:)-SS0(4,:) 
          AA(:,:) = transpose( SS0(1:3,1:3) )
          BB(:)   = delC(:)

          delSS   =  SQRT( sum( (SS0(1:3,1:3)-S0(1:3,1:3) )**2 ) ) / 9.0
          magSS   =  SQRT( sum( (SS0(1:3,1:3)+S0(1:3,1:3) )**2 ) ) / 9.0
          !write(6,*) iinr," SS Matrix difference ",delSS/magSS

          if ( ( delSS / magSS ) > 0.5 ) then
             !write( * , *) 'Thongs giing south with iteration '
             !STOP
          end if

          
               do i=1,3
                 !!write(6, '(4(f12.7,1x)," => ",4(f12.7,1x))') S0(:,I),SS0(:,I)
                end do
   
          call DGESV(3, 1, AA, 3, pivot, BB, 3, ok)

          delX(:) = BB(:)
          KK1 = KK0+delX
                !!write( 6, * ) '      ------ --- -- ',iinr,OK
                do i=1,3
                 !!write(6, '(3(f12.7,1x))') KK0(i),delX(i),KK1(i)
                end do
          KK0 = KK1
          dexx = SQRT( SUM(delX**2) )
          if (dexx < 0.000001 ) EXIT
       END do
       write( 6, * ) ' ---- ',iout
       do i=1,3
         write(6, '((10x,f12.7))') KK1(i)
       end do
       write(911) KK0
    end do
    
 end program rg_drv
