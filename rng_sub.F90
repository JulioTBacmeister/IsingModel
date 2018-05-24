subroutine rng_sub( x1,x2,x3, y1,y2,y3 )
  implicit none

  real, intent(in)  :: x1,x2,x3
  real, intent(out) :: y1,y2,y3

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real    :: z(-1:1,-1:1,-1:1,-1:1)
  integer :: a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p
  integer :: i1,i2,i3,i4
  real    :: NN,NNN,NNNN
  real    :: Econfig, boltzmann_wgt,w1,w2,w3
  integer :: iblock

  z(:,:,:,:) = 0.
  
  w1=exp(x1)
  w2=exp(x2)
  w3=exp(x3)

  aloop: do a = -1,1,2
  bloop: do b = -1,1,2
  cloop: do c = -1,1,2
  dloop: do d = -1,1,2
  eloop: do e = -1,1,2
  floop: do f = -1,1,2
  gloop: do g = -1,1,2
  hloop: do h = -1,1,2
  iloop: do i = -1,1,2
  jloop: do j = -1,1,2
  kloop: do k = -1,1,2
  lloop: do l = -1,1,2
  mloop: do m = -1,1,2
  nloop: do n = -1,1,2
  oloop: do o = -1,1,2
  ploop: do p = -1,1,2

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

     
  end do ploop
  end do oloop
  end do nloop
  end do mloop
  end do lloop
  end do kloop
  end do jloop
  end do iloop
  end do hloop
  end do gloop
  end do floop
  end do eloop
  end do dloop
  end do cloop
  end do bloop
  end do aloop


  Y1 = log ( z( 1, 1, 1, 1) / z( 1,-1,-1, 1) ) / 8.
  Y2 = log ( z( 1, 1, 1, 1) / z( 1, 1,-1,-1) ) / 4.-Y1
  Y3 = log ( z( 1, 1,-1,-1) / z( 1, 1, 1,-1) ) / 2.+Y2
  write(*,*) "Y1=",y1 ,"Y2=",y2 ,"Y3=",y3
  

end subroutine rng_sub
