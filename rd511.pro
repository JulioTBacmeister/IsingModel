close,1
openr,1,/f77,'fort.511'

nc=0L
nx=0L
K1=1.d
T=1.d

readu,1,nc,nx,K1,T
print,K1,T

;ccc=lonarr(nx+2,nx+2,nc)
ranu=dblarr(nx,nx)
mc1=dblarr(nc)
mc2=dblarr(nc)
mc3=dblarr(nc)
dE_x=dblarr(nc)
dE_p=dblarr(nc)
MnMag=dblarr(nc)
Energy=dblarr(nc)
ix=lonarr(nc)
jy=lonarr(nc)
ccc_sum = lonarr(nc)

E_init=1.d
E_final=1.d

readu,1,ranu
readu,1,mc1,mc2,mc3,ix,jy
readu,1,E_init,E_final
readu,1,dE_x,dE_p,MnMag,Energy


close,1

end
  
