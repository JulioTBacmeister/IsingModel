pro corrftn,ccc,corr=corrs,count=count,smean=smean,svar=svar

  s=size(ccc)

  minn=min( s(1:2) )
  corrs = dblarr(minn)
  count = lonarr(minn)
  x2=s(1)/2
  y2=s(2)/2

  smean = total( ccc(1:s(1)-2,1:s(2)-2) ) /( (s(1)-2)*(s(1)-2) )
  svar  = total( (ccc(1:s(1)-2,1:s(2)-2) - smean )^2 ) /( (s(1)-2)*(s(1)-2) )

  for y2=0,s(2)-2 do begin
  for x2=0,s(1)-2 do begin
  for j=0,s(2)-2 do begin
     for i=0,s(1)-2 do begin
         r=sqrt( 1.*(i-x2)^2 + 1.*(j-y2)^2 )
         ir=fix(r)
         if ir le minn-2 then begin
         ss= (ccc(i,j)-smean)*(ccc(x2,y2)-smean)
         corrs(ir) = corrs(ir)+double(ss)
         count(ir) = count(ir)+1
         endif
  endfor
  endfor
  endfor
  endfor


  
  return
end
