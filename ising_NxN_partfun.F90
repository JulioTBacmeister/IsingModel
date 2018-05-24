subroutine ising_NxN_partfun( NX, x1,x2,x3) ! , ZZZc, CCC  )
  implicit none

  real, intent(inout)  :: x1,x2,x3
  integer, intent(in) :: NX
  !integer, intent(out) :: CCC(Nx,Nx,2**(Nx*Nx) )
  !real, intent(out) :: ZZZc( 2**(Nx*Nx) )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real    :: z(-1:1,-1:1,-1:1,-1:1)
  integer :: a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,s
  integer :: i1,i2,i3,i4,Nspin,Nconf, ii,jj,kk
  real    :: NN,NNN,NNNN
  real    :: Econfig, boltzmann_wgt,w1,w2,w3,ZZZ,ZZp,mm1,mm2
  integer :: iblock
  integer :: CCC(Nx,Nx,2**(Nx*Nx) ),CC(Nx*Nx,2**(Nx*Nx)),ccc_sum(2**(NX*NX))
  real    :: ZZZc( 2**(Nx*Nx) ) , EEc( 2**(Nx*Nx) )
  real    :: corr(NX,NX,NX,NX),tps
  
  w1=exp(x1)
  w2=exp(x2)
  w3=exp(x3)


  ! Input coupling coefficients are assumed to be in the
  ! form K = J/(kb*T) where J's come from Hamiltonian
  !
  !       H = SUM_(ij=NN) J*S_i*S_j
  !
  ! Critical value: K=0.4407... (analytical, Onsager)
  !

  write(*,*) "T per site=",1./X1," K1=",X1

  !NX=4
  Nspin=NX*NX
  Nconf=2**Nspin
  !allocate( CC(Nspin,Nconf ))
  !allocate( CCC(NX,NX,Nconf) )
  !allocate( ZZZc(Nconf) )
  do i=1,Nspin
  do c=1,Nconf
     L=2**(i-1)
     s=(-1)**INT( (c-1)/L )
     cc( i, c ) = s
  end do
  end do

  CCC = RESHAPE( CC, (/ NX,NX,Nconf /) )

  ZZZ = 0.
  
  ! 2x2
  !---------
  ! a b 
  ! c d 
  !---------
  do kk=1,Nconf
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
     EEc(kk) = Econfig
  end do


  write( 211 ) Nconf,Nx
  write( 211 ) EEc,zzzc
  !write( 411 ) CCC
  

  !write(*,*) "ROutine Ising NxN "
  !write(*,*) "Partition Function ZZZ=",ZZZ


  ! calculate mean magnetization
  mm1=0.
  mm2=0.
  do kk=1,nconf
     !write(*,*) sum( CCC(:,:,KK) ), ZZZc(KK)/ZZZ
     ccc_sum(kk) = abs(sum( CCC(:,:,KK) ))
     mm1 = mm1+abs(sum( CCC(:,:,KK) ))
     mm2 = mm2+ZZZc(KK) * abs(sum( CCC(:,:,KK) ) )/ ZZZ
     !mm2 = mm2+ZZZc(KK) / ZZZ
  end do

  mm1 = mm1/ (NX*NX) / Nconf
  mm2 = mm2/ (NX*NX)
  
  write(*,*) " Expectation value ABS magnetization "
  write(*,*)  " : ",mm2
  write(*,*) " Mean value ABS magnetization over all possible states"
  write(*,*)  " : ",mm1

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


  write( 211 ) ccc_sum

end subroutine ising_NxN_partfun
