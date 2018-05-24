#define prob4x4
program minrg_drv
!====================================================================
! Minimalistic 4x4 ==> 2x2 real-space renormalization calculation
! based on Willis et al. (2015, doi:10.3389/fphy.2015.00046)
! Assumes periodic BC's, no external magnetic field.
!====================================================================  
  use renorm_utils, only : zzNxNp,blockrg
  
  implicit none
  real    :: x1,x2,x3,y1,y2,y3,T,Tc,dexx,w1,w2,w3
  integer :: niter,i,ix2,l,ipara
  integer :: OK,c4x4(0:5,0:5, 2**16 )
  integer :: c2x2(0:3,0:3, 2**16 ),cpara(2,2,4)
  real :: b4x4(2**16),zzpara( 4 ),delCC
#ifdef prob4x4
  real :: AA(4,4),BB(4),KK0(4),KK1(4)
  real :: AAx(4,4)
  integer :: pivot(4)
#else  
  real :: AA(3,3),BB(3),KK0(4),KK1(4)
  real :: AAx(3,3)
  integer :: pivot(3)
#endif
  integer :: ntemps,itemp,nconf,noutr,ioutr,kk

  nconf=2**16
  ntemps=50
  noutr=8
  cpara=0

  ! Set up 4 paradigm 2x2 configurations.  In the
  ! absence of external magnetic field these
  ! cover all 16 possibilities
  cpara(:,1,1) = (/  1,  1 /) 
  cpara(:,2,1) = (/  1,  1 /)
  !!!!
  cpara(:,1,2) = (/ -1,  1 /) 
  cpara(:,2,2) = (/  1,  1 /) 
  !!!!
  cpara(:,1,3) = (/ -1, -1 /) 
  cpara(:,2,3) = (/  1,  1 /) 
  !!!!
  cpara(:,1,4) = (/ -1,  1 /) 
  cpara(:,2,4) = (/  1, -1 /) 
   write(6,'(2(i2,1x))') cpara

#ifdef prob4x4
  ! Set up LHS minimal RG matrix for 4x4 problem
  !   AA * K' = ln( Z_para )
  !------------------------------------
  AAx(:,1) = (/  1.,  8.,  8.,  4. /) 
  AAx(:,2) = (/  1.,  0.,  0., -4. /)
  AAx(:,3) = (/  1.,  0., -8.,  4. /)
  AAx(:,4) = (/  1., -8.,  8., -4. /)
#else
  ! Set up LHS minimal RG matrix for 3x3 problem
  !   AA * K' = ln( Z_para(1:3) /Z_para(0) )
  !------------------------------------
  AAx(:,1) = (/  8.,  8.,  8. /)
  AAx(:,2) = (/  8., 16.,  0. /)
  AAx(:,3) = (/ 16.,  0.,  0. /)
#endif
  
     tc=2./( alog(1.+SQRT(2.)))  ! T_c ~ 2.2691853142130221
     x2=0.
     x3=0.
     write(911)noutr,noutr,ntemps


! Loop over SECOND (K2, NNN) coupling coefficient in
! Ising Hamiltonian
x2loop: do ix2=1,2
if (ix2==1) x2=0.
if (ix2==2) x2=0.15

! Loop over non-dim temps (actually over
! reduced coupling K1=beta*V1)
temps: do itemp=1,ntemps
        
     x1=.1 + ( dble(itemp)/dble(ntemps) )*0.9
     T=1/x1
     write(*,*) " T , T_c ",T,Tc
     ! Initialize vector of couplings for current T
     ! Note:  K0=0. and K3=0.
     KK0 = (/ 0.0, X1, X2, X3 /)
     write(*,*) kk0(2:4)
     write(911) kk0(2:4),kk0(1)

     !====================================
     ! Now iterate:
     !     1) Create c4x4 (configurations) and b4x4 (Boltzmann weights) on 4x4 lattice
     !        with given K-vector
     !     2) Coarse grain/block spin to 2x2 lattice
     !     3) Accumulate partial partition functions (ZZpara) for each of 4 2x2 paradigm configs (cpara's)
     !     4) Solve  AA*K = B for new K-vector
     !     5) Return to step 1
     !=====================================
     ! Notes:
     !     *  c4x4 does not actually change as zzNxNp is called. Just wastefully re-calculated
     !        each time.  b4x4 does change as K's change
     !     *  Problem can either be done as 4x4 or 3x3.  With 4x4, K(0) is explicitly calculated
     !        during each iteration.  Value of K(0) does not exhibit fixed point, but does not
     !        interfere with fixed point convergence of K(1:3). For 3x3 problem, RHS is normalized
     !        by ZZpara of config(1,1,1,1).
     !     *  LAPACK solvers want AA matrix as transpose of the way it is naturally in fortran.
     !        Spent many days puzzling over why nothing worked until discovering this ....
     !=====================================     
     outer: do ioutr=1,noutr
     call zzNxNp( 4 , KK0(1), KK0(2),KK0(3),KK0(4) , b4x4 , c4x4 ) 

     zzpara = 0.
     do kk=1,Nconf
        ! Block-spin transform of c4x4
        c2x2(:,:,kk) = blockrg( c4x4(:,:,kk) , 4 , 2 )
        do ipara=1,4
           ! Match c2x2's to paradigm 2x2 configs
           delCC = sum( dble( c2x2(1:2,1:2,kk)-cpara(1:2,1:2,ipara) )**2 )
           if (delCC < 1.e-6 ) then
              ! Accumulate partial partition functions
              zzpara(ipara) = zzpara(ipara)+b4x4(kk)
           end if
        end do
     end do
     !!!!!!!!!!!!!!!!!****************
     !!!!!!!!!!!!!!!!!!!!!************
     !  NOTE that LAPACK subr's seem
     !  to need *TRANSPOSE* of the way A
     !  is naturally written above
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     AA = transpose(AAx)  !!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef prob4x4
          BB(:)   = log( ZZpara )
          call DGESV(4, 1, AA, 4, pivot, BB, 4, ok)
          KK0(:) = BB(:)
#else
          BB(1)   = log( ZZpara(1) ) - log( ZZpara(2) )
          BB(2)   = log( ZZpara(1) ) - log( ZZpara(3) )
          BB(3)   = log( ZZpara(1) ) - log( ZZpara(4) )
          call DGESV(3, 1, AA, 3, pivot, BB, 3, ok)
          KK0(2:4) = BB(1:3)
#endif
#if 0
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !  Various direct solutions of Ax=B
          Y1 = ( log(zzpara(1))-log(zzpara(4)) )/16.
          Y2 = ( log(zzpara(1))-log(zzpara(3)) - 8.*Y1 )/16.
          !Y3 = ( log(zzpara(3))-log(zzpara(2)) + 8.*Y2 )/4.
          Y3 = ( log(zzpara(3))-log(zzpara(2)) + 8.*Y2 )/8.

          !W1 = ( log( ZZpara(1) ) - log( ZZpara(4) ) )/16.
          !W2 = ( log( ZZpara(1) ) - log( ZZpara(3) ) - 8.*W1 )/16.
          !W3 = ( log( ZZpara(1) ) - log( ZZpara(2) ) - 8.*W1 - 8.*W2 )/8.
          W1 = ( log( ZZpara(1) ) - log( ZZpara(4) ) )/aax(1,3)
          W2 = ( log( ZZpara(1) ) - log( ZZpara(3) ) - aax(1,2)*W1 )/aax(2,2)
          W3 = ( log( ZZpara(1) ) - log( ZZpara(2) ) - aax(1,1)*W1 - aax(2,1)*W2 )/aax(3,1)
          KK0 = (/ 0.0, Y1, Y2, Y3 /)
#endif
#ifdef prob4x4          
          write(*,*) KK0(1:4)
#else
          write(*,*) KK0(2:4)
#endif          
          write(911) KK0(2:4),KK0(1)
    end do outer
    end do temps
end do x2loop

end program minrg_drv
