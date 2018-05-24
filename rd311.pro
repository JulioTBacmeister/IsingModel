close,1
openr,1,/f77,'fort.311'

nc=0L
nx=0L


readu,1,nc,nx

ccc=lonarr(nx+2,nx+2,nc)
ranu=dblarr(nx,nx)
mc1=dblarr(nc)
mc2=dblarr(nc)
mc3=dblarr(nc)
deec=dblarr(nc)
deex=dblarr(nc)
eec=dblarr(nc)
zzzc=dblarr(nc)
ix=lonarr(nc)
jy=lonarr(nc)
ccc_sum = lonarr(nc)

readu,1,ranu
readu,1,mc1,mc2,mc3
readu,1,ix,jy,deec,deex
readu,1,ccc
readu,1,eec,zzzc,ccc_sum


close,1


ccy = total( total( ccc(1:nx-2,1:nx-2,*),1) ,1 )
end
  
