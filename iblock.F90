function iblock(a,b,c,d) result (blockspin)
  implicit none
  integer, intent(IN) :: a,b,c,d
  integer :: blockspin

  integer :: spinsum

  spinsum  = a + b + c + d

  if (spinsum == 0) then
     blockspin = a / abs(a)
  else
     blockspin = spinsum / abs( spinsum )
  endif
  
end function iblock
