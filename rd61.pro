
N=3
a=dblarr(N,N) & b=dblarr(N) & x=dblarr(N)

  close,1
  openr,1,'fort.61',/f77_u

  readu,1,a,b
  readu,1,x

  close,1
  end
