library(ggplot2)
22
mutationrate<- function (n,snp) {
  a=0
  for(i in 1:(n-1)) {
    a= a+ 1/i
  }
  print(snp/a)
}

mutationrate(dim(ceu)[1],dim(ceu)[2])


--
  hwe.calc=function(x)
  {
    geno.X=apply(cbind(x[seq(1,length(x),2)],x[seq(2,length(x),2)]),1,sum)
    p_X=length(x[x==0])/length(x)
    expected.vals=c(p_X^2,2*p_X*(1-p_X),(1-p_X)^2)*length(geno.X)
    observed.vals=c(length(geno.X[geno.X==0]),length(geno.X[geno.X==1]),length(geno.X[geno.X==2]))
    pearson.chi.sq=sum(((observed.vals-expected.vals)^2)/expected.vals)
    p.val=1-pchisq(pearson.chi.sq,1)
    #return(pearson.chi.sq)
    return(p.val)
  }
  num.snpsceu= dim(ceu)[2]
  hweceu=rep(NA,num.snpsceu)
  for(i in 1: num.snpsceu) hweceu[i]= hwe.calc(ceu[,i])
  hist(hweceu, main="15 chr CEU",xlab="HWE2 scores across snps")
  
  num.snpsceu20= dim(ceu20)[2]
  hweceu20=rep(NA,num.snpsceu20)
  for(i in 1: num.snpsceu20) hweceu20[i]= hwe.calc(ceu20[,i])
  hist(hweceu20, main="20 chr CEU",xlab="HWE2 scores across snps")
  
  num.snpspop15= dim(pop15)[2]
  hwepop15=rep(NA,num.snpspop15)
  for(i in 1: num.snpspop15) hwepop15[i]= hwe.calc(pop15[,i])
  hist(hwepop15, main="15 chr POP",xlab="HWE2 scores across snps")
  
  num.snpspop20= dim(pop20)[2]
  hwepop20=rep(NA,num.snpspop20)
  for(i in 1: num.snpspop20) hwepop20[i]= hwe.calc(pop20[,i])
  hist(hwepop20, main="20 chr ",xlab="HWE2 scores across snps")
  
  
  hwe.snps.ceu=hwe.snps.ceu20=hwe.snps.pop15=hwe.snps.pop20=rep(NA,num.snps)
  for (i in 1:num.snps) hwe.snps.ceu[i]=hwe.calc(ceu[,i])
  for (i in 1:num.snps) hwe.snps.ceu20[i]=hwe.calc(ceu20[,i])
  for (i in 1:num.snps) hwe.snps.pop15[i]=hwe.calc(pop15[,i])
  for (i in 1:num.snps) hwe.snps.pop20[i]=hwe.calc(pop20[,i])
  par(mfrow=c(2,2)) ## plots 2 figures per row, and 2 per column
  hist(hwe.snps.ceu,main="CEU",xlab="HWE scores across SNPs")
  hist(hwe.snps.ceu20,main="CEU20",xlab="HWE scores across SNPs")
  hist(hwe.snps.pop15,main="POP15",xlab="HWE scores across SNPs")
  hist(hwe.snps.pop20,main="POP20",xlab="HWE scores across SNPs")
  par(mfrow=c(1,1)) ## goes back to one figure total
  boxplot(hwe.snps.ceu,hwe.snps.ceu15,hwe.snps.pop15,hwe.snps.pop20,
          ylab="HWE p-values",names=c("CEU","CEU20","POP15","POP20"))
  
  
  
  ---
  #Linkage disequilibrium question
    
    d.prime.calc=function(x,y)
    {
      D.00=length(x[x==0 & y==0])/length(x)-(length(x[x==0])/
                                               length(x))*(length(y[y==0])/length(y))
      D.minus=min((length(x[x==1])/length(x))*(length(y[y==1])/length(y)),
                  (length(x[x==0])/length(x))*(length(y[y==0])/length(y)))
      D.plus=min((length(x[x==1])/length(x))*(length(y[y==0])/length(y)),
                 (length(x[x==0])/length(x))*(length(y[y==1])/length(y)))
      if (D.00>=0) D.prime=D.00/D.plus
      if (D.00<0) D.prime=D.00/D.minus
      return(abs(D.prime))
    }
  
  #D' Linkage disequilibrium for ceu (chr15)
  
  D_ceu=matrix(NA,nrow=dim(ceu)[2],ncol=dim(ceu)[2]) 
  for (i in 1:(dim(ceu)[2]-1)){
    for(j in (i+1):dim(ceu)[2]){
      D_ceu[i,j]=d.prime.calc(ceu[,i],ceu[,j])
    }
  }
  
  #### plot scatter D link 
  distance_ceu15 = c()
  d_prime_ceu15= c() 
  for (i in 1:(dim(ceu)[2]-1)){
    for (j in  (i+1):dim(ceu)[2]){
      distance_ceu15 = c(distance_ceu15,j-i)
      d_prime_ceu15 = c(d_prime_ceu15,D_ceu[i,j] )
    }
  }
  par(mfrow=c(1,1))
  Data <- data.frame(distance_ceu15, d_prime_ceu15) 
  ggplot(Data,aes(distance_ceu15, d_prime_ceu15))+geom_point()+geom_smooth()
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #r2 Linkage disequilibrium for ceu (chr15)
  
  r2_ceu=matrix(NA,nrow=dim(ceu)[2],ncol=dim(ceu)[2]) 
  for (i in 1:(dim(ceu)[2]-1)){
    for(j in (i+1):dim(ceu)[2]){
      r2_ceu[i,j]=cor(ceu[,i],ceu[,j])^2
    }
  }
  
  
  #### plot scatter r2 link ceu15 
  distance_r2ceu15 = c()
  r2ceu15= c() 
  for (i in 1:(dim(ceu)[2]-1)){
    for (j in  (i+1):dim(ceu)[2]){
      distance_r2ceu15 = c(distance_r2ceu15,j-i)
      r2ceu15 = c(r2ceu15,r2_ceu[i,j] )
    }
  }
  par(mfrow=c(1,1))
  Data <- data.frame(distance_r2ceu15, r2ceu15) 
  ggplot(Data,aes(distance_r2ceu15, r2ceu15))+geom_point()+geom_smooth()
  
  
  
  
  #D' Linkage disequilibrium for ceu20 (chr20)
  
  D_ceu20=matrix(NA,nrow=dim(ceu20)[2],ncol=dim(ceu20)[2]) 
  for (i in 1:(dim(ceu20)[2]-1)){
    for(j in (i+1):dim(ceu20)[2]){
      D_ceu20[i,j]=d.prime.calc(ceu20[,i],ceu20[,j])
    }
  }
  
  
  #### plot scatter D link 
  distance_ceu20 = c()
  d_prime_ceu20= c() 
  for (i in 1:(dim(ceu20)[2]-1)){
    for (j in  (i+1):dim(ceu20)[2]){
      distance_ceu20 = c(distance_ceu20,j-i)
      d_prime_ceu20 = c(d_prime_ceu20,D_ceu20[i,j] )
    }
  }
  par(mfrow=c(1,1))
  Data <- data.frame(distance_ceu20, d_prime_ceu20) 
  ggplot(Data,aes(distance_ceu20, d_prime_ceu20))+geom_point()+geom_smooth()
  
  
  
  #r2 Linkage disequilibrium for ceu20 (chr20)
  
  r2_ceu20=matrix(NA,nrow=dim(ceu20)[2],ncol=dim(ceu20)[2]) 
  for (i in 1:(dim(ceu20)[2]-1)){
    for(j in (i+1):dim(ceu20)[2]){
      r2_ceu20[i,j]=cor(ceu20[,i],ceu20[,j])^2
    }
  }
  
  
  #### plot scatter r2 link ceu20 
  distance_r2ceu20 = c()
  r2ceu20= c() 
  for (i in 1:(dim(ceu20)[2]-1)){
    for (j in  (i+1):dim(ceu20)[2]){
      distance_r2ceu20 = c(distance_r2ceu20,j-i)
      r2ceu20 = c(r2ceu20,r2_ceu20[i,j] )
    }
  }
  par(mfrow=c(1,1))
  Data <- data.frame(distance_r2ceu20, r2ceu20) 
  ggplot(Data,aes(distance_r2ceu20, r2ceu20))+geom_point()+geom_smooth()
  
  
  
  
  #D' Linkage disequilibrium for pop15 (pop15)
  
  D_pop15=matrix(NA,nrow=dim(pop15)[2],ncol=dim(pop15)[2]) 
  for (i in 1:(dim(pop15)[2]-1)){
    for(j in (i+1):dim(pop15)[2]){
      D_pop15[i,j]=d.prime.calc(pop15[,i],pop15[,j])
    }
  }
  
  #### plot scatter D link 
  distance_pop15 = c()
  d_prime_pop15= c() 
  for (i in 1:(dim(pop15)[2]-1)){
    for (j in  (i+1):dim(pop15)[2]){
     distance_pop15 = c(distance_pop15,j-i)
     d_prime_pop15 = c(d_prime_pop15,D_pop15[i,j] )
    }
  }
  par(mfrow=c(1,1))
  Data <- data.frame(distance_pop15, d_prime_pop15) 
  ggplot(Data,aes(distance_pop15, d_prime_pop15))+geom_point()+geom_smooth()
  
  #r2 Linkage disequilibrium for ceu20 (chr20)
  
  r2_pop15=matrix(NA,nrow=dim(pop15)[2],ncol=dim(pop15)[2]) 
  for (i in 1:(dim(pop15)[2]-1)){
    for(j in (i+1):dim(pop15)[2]){
      r2_pop15[i,j]=cor(pop15[,i],pop15[,j])^2
    }
  }
  
  
  
  #### plot scatter r2 link pop15 
  distance_r2pop15 = c()
  r2pop15= c() 
  for (i in 1:(dim(pop15)[2]-1)){
    for (j in  (i+1):dim(pop15)[2]){
      distance_r2pop15 = c(distance_r2pop15,j-i)
      r2pop15 = c(r2pop15,r2_pop15[i,j] )
    }
  }
  par(mfrow=c(1,1))
  Data <- data.frame(distance_r2pop15, r2pop15) 
  ggplot(Data,aes(distance_r2pop15, r2pop15))+geom_point()+geom_smooth()
  
  
  
  
  #D' Linkage disequilibrium for pop20 (pop20)
  
  D_pop20=matrix(NA,nrow=dim(pop20)[2],ncol=dim(pop20)[2]) 
  for (i in 1:(dim(pop20)[2]-1)){
    for(j in (i+1):dim(pop20)[2]){
      D_pop20[i,j]=d.prime.calc(pop20[,i],pop20[,j])
    }
  }
  
  #### plot scatter D link 
  distance_pop20 = c()
  d_prime_pop20= c() 
  for (i in 1:(dim(pop20)[2]-1)){
    for (j in  (i+1):dim(pop20)[2]){
      distance_pop20 = c(distance_pop20,j-i)
      d_prime_pop20 = c(d_prime_pop20,D_pop20[i,j] )
    }
  }
  par(mfrow=c(1,1))
  Data <- data.frame(distance_pop20, d_prime_pop20) 
  ggplot(Data,aes(distance_pop20, d_prime_pop20))+geom_point()+geom_smooth()
  
  
  
  #r2 Linkage disequilibrium for ceu20 (chr20)
  
  r2_pop20=matrix(NA,nrow=dim(pop20)[2],ncol=dim(pop20)[2]) 
  for (i in 1:(dim(pop20)[2]-1)){
    for(j in (i+1):dim(pop20)[2]){
      r2_pop20[i,j]=cor(pop20[,i],pop20[,j])^2
    }
  }
  
  
  #### plot scatter r2 link pop20 
  distance_r2pop20 = c()
  r2pop20= c() 
  for (i in 1:(dim(pop20)[2]-1)){
    for (j in  (i+1):dim(pop20)[2]){
      distance_r2pop20 = c(distance_r2pop20,j-i)
      r2pop20 = c(r2pop20,r2_pop20[i,j] )
    }
  }
  par(mfrow=c(1,1))
  Data <- data.frame(distance_r2pop20, r2pop20) 
  ggplot(Data,aes(distance_r2pop20, r2pop20))+geom_point()+geom_smooth()
  
  
  
    
  #D' Linkage disequilibrium for ceu&ceu20 (chr15&20)
  
  D_ceu1520=matrix(NA,nrow=dim(ceu)[2],ncol=dim(ceu20)[2]) 
  for (i in 1:(dim(ceu)[2]-1)){
    for(j in (i+1):dim(ceu20)[2]){
      D_ceu1520[i,j]=d.prime.calc(ceu[,i],ceu20[,j])
    }
  }
  
  hist(D_ceu1520, main="Linkage disequilibrium for ceu&ceu20",xlab="Disequilibrium range")
  
  
  #r2 Linkage disequilibrium for ceu&ceu20 (chr15&20)
  
  r2_ceu1520=matrix(NA,nrow=dim(ceu)[2],ncol=dim(ceu20)[2]) 
  for (i in 1:(dim(ceu)[2]-1)){
    for(j in (i+1):dim(ceu20)[2]){
      r2_ceu1520[i,j]=cor(ceu[,i],ceu20[,j])^2
    }
  }
  
  #D' Linkage disequilibrium for pop15&pop20 (pop15&20)
    
  D_pop1520=matrix(NA,nrow=dim(pop15)[2],ncol=dim(pop20)[2]) 
  for (i in 1:(dim(pop15)[2]-1)){
    for(j in (i+1):dim(pop20)[2]){
      D_pop1520[i,j]=d.prime.calc(pop15[,i],pop20[,j])
    }
  }
  
  hist(D_ceu1520, main="Linkage disequilibrium for pop15&pop20",xlab="Disequilibrium range")
  
  #r2 Linkage disequilibrium for pop15&pop20 (pop15&20)
  
  r2_pop1520=matrix(NA,nrow=dim(pop15)[2],ncol=dim(pop20)[2]) 
  for (i in 1:(dim(pop15)[2]-1)){
    for(j in (i+1):dim(pop20)[2]){
      r2_pop1520[i,j]=cor(pop15[,i],pop20[,j])^2
    }
  }  
    
  



