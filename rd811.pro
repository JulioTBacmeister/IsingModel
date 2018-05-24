close,1
openr,1,/f77,'fort.811'

ncx=10L^3
nc=0L & nx=0L & nxb=0L & B=0L
K1=0.d
T1=0.d
readu,1,nc,nx,b,nxb,K1,T1
nx=nx+2
nxb=nxb+2

ccc=lonarr( nx, nx, ncx )
ccc0=lonarr(nx,nx)
bccc=lonarr( nxb, nxb, ncx )
bccc0=lonarr(nxb,nxb)
kk=lonarr(ncx)

snml = dblarr(4,3,ncx)
snms = dblarr(4,3,ncx)
snml0 = dblarr(4,3)
snms0 = dblarr(4,3)

kk0=0L
i=0
while not eof(1) do begin

   readu,1,ccc0,bccc0,kk0,snml0,snms0
   ccc(*,*,i) = ccc0
   bccc(*,*,i) = bccc0
   snml(*,*,i) = snml0
   snms(*,*,i) = snms0
   kk(i)=kk0
   i=i+1
endwhile


ccc=ccc(*,*,0:i-1)
bccc=bccc(*,*,0:i-1)
ncx=i-1

close,1
STOP
window,0,re=2,xs=400,ys=400

for i=0,ncx do begin tv,congrid(ccc(*,*,i),350,350),30,30 & wait,.05 & endfor  


end
  
