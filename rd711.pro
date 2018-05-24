close,1
openr,1,/f77,'fort.711'


ncx=10L^3
nx2=0L
nx=0L
readu,1,nc,nx
nx2=nx/2
nx2=nx2+1

c3  =dblarr( 3,ncx)
c30 =dblarr(3 )
cv  =dblarr( nx2 ,ncx)
cv0 =dblarr(nx2 )
cn  =lonarr( nx2 ,ncx)
cn0 =lonarr(nx2 )
kk  =lonarr(ncx)
kk0 =0L
i=0
while not eof(1) do begin

   readu,1,cv0,cn0,kk0,c30
   cv(*,i) = cv0
   cn(*,i) = cn0
   c3(*,i) = c30
   kk(i)=kk0
   i=i+1
endwhile

kk=kk(0:i-1)
cn=cn(*,0:i-1)
cv=cv(*,0:i-1)

end
