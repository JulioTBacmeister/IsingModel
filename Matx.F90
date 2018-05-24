function Matx(a,b) result (c)
  implicit none
  real, intent(IN) :: a(:,:),b(:,:)
  integer :: n,m,L,xa,ya,xb,yb,i,j

  real :: c( size(a,2) ,size(b,1) )

  xa=size(a,1)
  ya=size(a,2)
  xb=size(b,1)
  yb=size(b,2)

  if ( Yb /= Xa ) then
     c(:,:) = -9e12
  end if
  
  do j=1,Ya
     do i=1,Xb
        c(i,j)=0.
        do n=1,Xa
           c(i,j) = c(i,j)+a(n,j)*b(n,i)
        end do
     end do
  end do
  
end function Matx
