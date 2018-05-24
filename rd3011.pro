close,1
openr,1,/f77,'fort.3011'

nconf=0L & cpara=lonarr(2,2,4)

readu,1,nconf,cpara

c2x2=lonarr(4,4,nconf) & paraci = lonarr( nconf )

readu,1,c2x2,paraci
close,1
end
