close,1
openr,1,/f77,'fort.2011'

nx=6 & nxb=3
nx=nx+2
nxb=nxb+2

ccc=lonarr( nx, nx, 10000L)
bcc=lonarr( nxb, nxb, 10000L)

a=lonarr( nx, nx)
ab=lonarr(nxb,nxb)
window,xs=450,ys=400
i=0L
while not eof(1) do begin

   readu,1,a,ab
   ccc(*,*,i) = a
   bcc(*,*,i) = ab
   ;tv,congrid(ccc(*,*,i),350,350),30,30
   ;wait,.025
   i=i+1L
endwhile
ccc=ccc(*,*,0:i-1)
bcc=bcc(*,*,0:i-1)
close,1

end
