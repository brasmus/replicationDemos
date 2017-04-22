shake <- function(x) {
    x <- x[order(rnorm(length(x)))]
}

demoConf <- function(m=3,s=2,N=10000,
                     n=c(10,50,100,200,500,1000,5000)) {
  
  X <- rnorm(N,mean=m,sd=s)
  #n <- round(n/100)*100
  l <- length(n)

  par(bty="n",yaxt="n",xaxt="n",xpd=NA)
  plot(c(0,l)*150,m+s*c(-1,1),type="n",
       main="Estimate of mean and confidence levels for the estimate",
       xlab="size of subset",ylab="mean estimate")
  lines(range(n),rep(m,2),lwd=5,col="grey")
  lines(range(n),rep(m,2)+s,lty=2,col="grey")
  lines(range(n),rep(m,2)-s,lty=2,col="grey")

  legend(l*100,m+0.97*s,c(expression(paste(bar(x),"=")),m,
                          expression(paste(sigma,"=")),s),bty="n",
         horiz=TRUE,x.intersp=0.25)  

  #print(n)

  SS <- rep(NA,l*100); dim(SS) <- c(l,100)
  for (i in 1:l) {
    im <- i - 1
    #print(n[i])
    for (ii in 1:100) {
    # Shake the sample: sort the data randomly...
      X <- shake(X)
      M <- mean(X[1:n[i]])
      S <- sd(X[1:n[i]])
      SS[i,ii] <- 2*S/sqrt(n[i]-1)
      #print(c(M,S))
      lines(im*150+rep(ii,2),M+c(1,-1)*2*S/sqrt(n[i]-1),col="grey30")
      lines(im*150+ii+c(-2,2),rep(M+2*S/sqrt(n[i]-1),2),col="grey30")
      lines(im*150+ii+c(-2,2),rep(M-2*S/sqrt(n[i]-1),2),col="grey30")
      points(im*150+ii,M,pch=19,cex=0.5,col="grey30")
    }
    text(im*150+50,m-1.1*s,n[i])    
  }

  lines( (0:(l-1))*150+50, m+rowMeans(SS),col="red",lty=2)
  lines( (0:(l-1))*150+50, m-rowMeans(SS),col="red",lty=2)
  text(l*150,m-1.1*s,"sample st.dev.",col="grey",pos=3)
  text(l*150,m+1.1*s,"sample st.dev.",col="grey")
  text(l*150,m,"sample mean",col="grey",pos=3)
}

demoRange <- function(m=3,s=1) {
  # Add more and more data ...

  par(bty="n",yaxt="n",xaxt="n",xpd=NA)
  plot(c(0,21),c(0,7),type="n",
       main="Increasing the sample size",
       xlab="sample size",ylab="value")
  lines(c(0,21),rep(m,2),lwd=5,col="grey70")
  lines(c(0,21),rep(m,2)+s,lty=2,col="grey70")
  lines(c(0,21),rep(m,2)-s,lty=2,col="grey70")

  X <- NA
  for (i in 1:20) {
    x <- rnorm(5,mean=m,sd=s)
    X <- c(X,x)
    M <- mean(X,na.rm=TRUE)
    S <- sd(X,na.rm=TRUE)
    n <- sum(is.finite(X))
    points(rep(i,length(X)),X,pch=19,cex=0.7,col="grey")
    points(i,M,pch=19,cex=0.9)
    lines(rep(i,2),M+c(1,-1)*2*S/sqrt(n-1))
    lines(i+c(-0.2,0.2),rep(M+2*S/sqrt(n-1),2))
    lines(i+c(-0.2,0.2),rep(M-2*S/sqrt(n-1),2))
    text(i,0,n,cex=0.5)    
  }

  legend(2,7,c(expression(paste(bar(x),"=")),m,
               expression(paste(sigma,"=")),s),bty="n",
         horiz=TRUE,x.intersp=0.25)  

}
