if not keyword_set(plotonly) then begin
close,1
openr,1,/f77,'fort.911'
nx2loop=2
noutr=0L&ninnr=0L&ntemps=0L
kk0=dblarr(3)&kk00=1.d
readu,1,noutr,ninnr,ntemps
kk=dblarr( 3,noutr+1,ntemps*nx2loop)
kk_0x=dblarr( noutr+1,ntemps*nx2loop)

for j=0,nx2loop*ntemps-1 do begin
for i=0,noutr do begin

   ;print,'reading',i,j
   readu,1,kk0,kk00
   kk(*,i,j)=kk0
   kk_0x(i,j) = kk00
   ;print,'              got',i,j
endfor
endfor
endif

 plot,kk(0,*,50),kk(1,*,50),xr=[.2,.5],yr=[0,.2],/xst,/yst
 for n=0,nx2loop*ntemps-1 do oplot,kk(0,*,n),kk(1,*,n)                  

end
  
