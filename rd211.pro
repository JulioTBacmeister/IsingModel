close,1
openr,1,/f77,'fort.211'

nc=0L
nx=0L


readu,1,nc,nx

eec_5x5=dblarr( nc )
zzzc_5x5=dblarr( nc )
;ccc=lonarr(nx,nx,nc)
ccc_sum_5x5=lonarr(nc)
;cor=dblarr(nx,nx,nx,nx)

dummy=1.d
readu,1,eec_5x5,zzzc_5x5
readu,1,ccc_sum_5x5


end
  
