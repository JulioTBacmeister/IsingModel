subroutine ising_4x4_partfun( x1,x2,x3, ZZZc, CCC  )
  implicit none

  real, intent(in)  :: x1,x2,x3
  !integer, intent(in) :: NX
  integer, intent(out) :: CCC(4,4,2**16)
  real, intent(out) :: ZZZc(2**16)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real    :: z(-1:1,-1:1,-1:1,-1:1)
  integer :: a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,s
  integer :: i1,i2,i3,i4,Nspin,Nconf, ii,jj,kk,NX
  real    :: NN,NNN,NNNN
  real    :: Econfig, boltzmann_wgt,w1,w2,w3,ZZZ,ZZp
  integer :: iblock
  integer, allocatable :: CC(:,:)

  z(:,:,:,:) = 0.
  
  w1=exp(x1)
  w2=exp(x2)
  w3=exp(x3)


  NX=4
  Nspin=NX*NX
  Nconf=2**Nspin
  allocate( CC(Nspin,Nconf ))
  !!allocate( CCC(NX,NX,Nconf) )
  !!allocate( ZZZc(Nconf) )

  do i=1,Nspin
  do c=1,Nconf
     L=2**(i-1)
     s=(-1)**INT( (c-1)/L )
     cc( i, c ) = s
  end do
  end do

  CCC = RESHAPE( CC, (/ NX,NX,Nconf /) )


  ZZZ = 0.
  
  ! 4x4
  !---------
  ! a b e f 
  ! c d g h
  ! i j m n
  ! k l o p
  !---------
  do kk=1,Nconf

     a=ccc(1,1,kk)
     b=ccc(2,1,kk)
     c=ccc(1,2,kk)
     d=ccc(2,2,kk)
  
     e=ccc(3,1,kk)
     f=ccc(4,1,kk)
     g=ccc(3,2,kk)
     h=ccc(4,2,kk)
  
     i=ccc(1,3,kk)
     j=ccc(2,3,kk)
     k=ccc(1,4,kk)
     l=ccc(2,4,kk)
  
     m=ccc(3,3,kk)
     n=ccc(4,3,kk)
     o=ccc(3,4,kk)
     p=ccc(4,4,kk)
  
  ! configs

     NN= a*b + b*e + e*f &
       + c*d + d*g + g*h &
       + i*j + j*m + m*n &
       + k*L + L*o + o*p &
       ! columnwise
       + a*c + c*i + i*k &
       + b*d + d*j + j*L &
       + e*g + g*m + m*o &
       + f*h + h*n + n*p

     NNN= a*d + c*b + b*g + d*e + e*h + g*f &
       +  c*j + i*d + d*m + j*g + g*n + m*h  &
       +  i*L + k*j + j*o + L*m + m*p + o*n

     NNNN= a*b*c*d + b*e*d*g + e*f*g*h &
         + c*d*i*j + d*g*j*m + g*h*m*n &
         + i*j*k*L + j*m*L*o + m*n*o*p    

     Econfig = x1*NN + x2*NNN + x3*NNNN
     !boltzmann_wgt =  w1**nn * w2**nnn * w3**nnnn !   exp( Econfig )
     boltzmann_wgt =  exp( Econfig )

     ZZZ = ZZZ + boltzmann_wgt
     ZZZc(kk)=boltzmann_wgt
     
  end do
     

  write(*,*) "Partition Function ZZZ=",ZZZ

end subroutine ising_4x4_partfun
