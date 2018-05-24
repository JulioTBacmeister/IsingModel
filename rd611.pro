close,1
openr,1,/f77,'fort.611'

ncx=10L^3
nc=0L & nx=0L
K1=0.d
T1=0.d
readu,1,nc,nx,K1,T1
nx=nx+2

ccc=lonarr( nx, nx, ncx )
ccc0=lonarr(nx,nx)
kk=lonarr(ncx)
kk0=0L
i=0
while not eof(1) do begin

   readu,1,ccc0,kk0
   ccc(*,*,i) = ccc0
   kk(i)=kk0
   i=i+1
endwhile


ccc=ccc(*,*,0:i-1)
ncx=i-1

close,1

window,0,re=2,xs=400,ys=400

for i=0,ncx do begin tv,congrid(ccc(*,*,i),350,350),30,30 & wait,.05 & endfor  


end
  
