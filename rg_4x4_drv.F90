program rg_drv
  implicit none
  real    :: x1,x2,x3,y1,y2,y3,T,Tc,dexx
  integer :: niter,i,ix2,l
  integer :: NX,NC,NB,pivot(3), OK,itemp,ntemps,iinr,niinr,noutr,iout
  real :: S0(4,3),S1(4,3),SX(4,3),SXX(4,3),SS0(4,3), S0xx(4,3)
  real :: AA(3,3),delC(3), delX(3),BB(3),KK0(3),KK1(3),delSS,magSS

  ntemps=50
  niinr = 30 !100
  noutr = 10 !100
  write(911) noutr,niinr,ntemps


  
  
     tc=2./( alog(1.+SQRT(2.)))  ! T_c ~ 2.2691853142130221
     !T=  2.5 ! 2.2691853 ! 2.5  ! 3.0 !1.7 ! 2.00 !2.28  !0.5*i

x3=0.
x2loop: do ix2=1,2
if (ix2==1) x2=0.
if (ix2==2) x2=0.15   
main: do itemp=1,ntemps
        
     !x1= 0.423 ! 1./T
     x1=.1 + ( dble(itemp)/dble(ntemps) )*0.9
     !x1=.3 + ( dble(itemp)/dble(ntemps) )*0.45
     T=1/x1
     write(*,*) " T , T_c ",T,Tc

     !x2=0.
     !x3=0.
     KK0 = (/ X1, X2, X3 /)
     write(*,*) kk0
     write(911) kk0,0.
     
   ! Outer iteration to find fixed point of RG
   outer: do iout=1,noutr   
       !call rg_4x4_rlzn( KK0(1),KK0(2),KK0(3) , SXX, S0xx ) 
       !write(6,*) " 4x4 "
       !do l=1,3
       !   write( 6, '(4(f12.6,2x))') S0xx(:,l)
       !end do
       !call rg_NxNp_rlzn( 4 , 2, KK0(1),KK0(2),KK0(3) , SXX, S0 ) 
       !write(6,*) " NxNp "
       !do l=1,3
       !   write( 6, '(4(f12.6,2x))') S0(:,l)
       !end do
       NC=2*150000
       call rg_MCb_rlzn( NC, 4 , 2 ,  KK0(1),KK0(2),KK0(3)    , S0 )
     
       inner: do iinr=1,niinr

          call rg_NxNp_rlzn( 2 , 1, KK0(1),KK0(2),KK0(3) , SS0, SXX ) 

          delC(:) = S0(4,:)-SS0(4,:)
          ! SS0 should be symmetric so transpose shouldn't matter.
          ! Left as a reminder
          AA(:,:) = transpose( SS0(1:3,1:3)  )
          BB(:)   = delC(:)

          delSS   =  SQRT( sum( (SS0(1:3,1:3)-S0(1:3,1:3) )**2 ) ) / 9.0
          magSS   =  SQRT( sum( (SS0(1:3,1:3)+S0(1:3,1:3) )**2 ) ) / 9.0

          call DGESV(3, 1, AA, 3, pivot, BB, 3, ok)

          delX(:) = BB(:)
          KK1 = KK0+delX
          KK0 = KK1
          dexx = SQRT( SUM(delX**2) )
          if (dexx < 0.000001 ) EXIT
       END do inner
       write(*,*) KK0
       write(911) KK0,0.
    end do outer
 end do main
end do x2loop

 end program rg_drv
