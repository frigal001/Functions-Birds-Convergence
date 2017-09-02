Evar<-function(x)   {
  v<-as.numeric(ncol(x)-1)
  for(k in 1:nrow(x)) {
    y<-x[k,2:ncol(x)]
    a<-as.numeric(ncol(x)-1); b<-a; d<-a
    for(i in 1:(ncol(x)-1)) {
      a[i]<-if(y[i]!=0) log(y[i]) else 0 }
    S<-sum(y>0)
    for(j in 1:(ncol(x)-1))    {
      b[j]<-a[[j]]/S}
    c<-sum(b)
    for(i in 1:(ncol(x)-1)) {
      d[i]<-if(y[i]!=0) ((a[[i]]-c)^2/S) else 0 }
    f<-sum(d)
    v[k]<-(1-2/pi*atan(f))   }
  v }
