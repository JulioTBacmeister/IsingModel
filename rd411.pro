close,1
openr,1,/f77,'fort.411'

nc=0L
nx=0L


readu,1,nc,nx

eec=dblarr( nc )
ccc=lonarr(nx,nx,nc)
cor=dblarr(nx,nx,nx,nx)

dummy=1.d
readu,1,eec
readu,1,ccc
readu,1,cor


end
  
