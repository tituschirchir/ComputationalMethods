Option.Binomial.Price<- function(K,tau,S,r,N,u,d)
{
  dt <- tau/N
  p <- (exp(r*dt)-d)/(u-d)
  disc <- exp(-r*dt)
  AssetPrices <- c(S*u^N)
  for(i in 2:(N+1))
  {
    AssetPrices[i] <- AssetPrices[i-1] * d / u
  }
  #OptionPricesAtMaturity
  op<-c(max(AssetPrices[1]-K,0))
  optionPrices<- matrix(nrow = 2*N+1, ncol = N+1)
  for(i in 1:(N+1))
  {
    optionPrices[2*i-1, N+1]<-max(0, AssetPrices[i]-K)
  }
  #Cascade Down
  rows <- c()
  for(i in N:1){
    j<-N+2-i
    for(m in 1:i)
    {
      optionPrices[j,i]<- disc*(p*optionPrices[j-1,i+1] + (1-p)* optionPrices[j+1,i+1])
      j<-j+2
    }
  }
  optionPrices
}

op<-Option.Binomial.Price(K=100,tau=1,S=100,r=0.06,N=30,u=1.1,d=0.9091)
print(op)