subroutine rng_sub2( x1,x2,x3, y1,y2,y3 )
  implicit none

  real, intent(in)  :: x1,x2,x3
  real, intent(out) :: y1,y2,y3

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real    :: z(-1:1,-1:1,-1:1,-1:1)
  integer :: a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,s
  integer :: i1,i2,i3,i4,Nspin,Nconf,NX, ii,jj,kk
  real    :: NN,NNN,NNNN
  real    :: Econfig, boltzmann_wgt,w1,w2,w3,ZZZ,ZZp
  integer :: iblock
  integer, allocatable :: CC(:,:), CCC(:,:,:)

  z(:,:,:,:) = 0.
  
  w1=exp(x1)
  w2=exp(x2)
  w3=exp(x3)


  NX=4
  Nspin=NX*NX
  Nconf=2**Nspin
  allocate( CC(Nspin,Nconf ))
  allocate( CCC(NX,NX,Nconf) )

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

     i1 = iblock( a,b,c,d )
     i2 = iblock( e,f,g,h )
     i3 = iblock( i,j,k,L )
     i4 = iblock( m,n,o,p )

     z(i1,i2,i3,i4) =  z(i1,i2,i3,i4) + boltzmann_wgt 
     ZZZ = ZZZ + boltzmann_wgt
     
  end do
     

  Y1 = log ( z( 1, 1, 1, 1) / z( 1,-1,-1, 1) ) / 8.
  Y2 = log ( z( 1, 1, 1, 1) / z( 1, 1,-1,-1) ) / 4.-Y1
  Y3 = log ( z( 1, 1,-1,-1) / z( 1, 1, 1,-1) ) / 2.+Y2
  write(*,*) "Y1=",y1 ,"Y2=",y2 ,"Y3=",y3
  write(*,*) "Partition Function ZZZ=",ZZZ

  write(*,*) "Renormed ZZ            =",SUM( z )

end subroutine rng_sub2
