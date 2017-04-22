# Rasmus E. Benestad
# R-script to simulate a forced-damped oscllator using 4th order
# Runge-Kutta integration.

resonance <- function(x0=1,v0=0,t0=0,t1=1000,N=1000,
                      LHS=cos(2*pi*(1:1000)/75),
                      f=0.001,m=0.1,w0=0.07) {
  FOscillator(x0=x0,v0=v0,t0=0,t1=10,N=N,LHS=LHS,f=f,m=m,w0=w0) -> a
  invisible(a)
}

FOscillator<- function(x0=1,v0=0,t0=0,t1=1000,N=1000,main=NULL,
                       LHS="cos(2*pi*t/75)",f=0.001,m=0.1,w0=0.07,scal=TRUE) {
# It's a long time since I last used 4th order RungeKutta to simulate an ODE.  
# But I think I still remember the essence, and I used the information about
# forced damped damped oscillator from URL
# http://physics.tamuk.edu/~suson/html/4390/DiffEq.html

  require(deSolve)
  # To avoide warining:  no visible global function definition for ‘F.ext’
  F.ext <- 0


  if (is.vector(LHS)) {
    lhs <- LHS[1]
    for (i in 2:length(LHS)) lhs <- paste(lhs,LHS[i],sep=",")
    funccode <- paste("F.ext= function(t) { y=c(",lhs, 
                      "); x=seq(",t0,",",t1,",length=",length(LHS),
                      "); Y=approx(x,y,t,n=",N,",rule=2)$y; Y}")
  }
  if (is.character(LHS)) funccode <-
    paste("F.ext= function(t) { y=",LHS,"; y}")
  #print(funccode)

  eval(parse(text=funccode))
  
  dXdt <- function(t,X,parms) {
    x <- X[1]; v <- X[2]
    with(as.list(parms), {
      dvdt <- -f/m*v - w0^2*x  + F.ext(t)/m
      dxdt <- v
      dF <- F.ext(t)
      Y <- c(x=dxdt,dxdt=dvdt,forcing=dF)
      list(Y)
    } )

  }


  
# rk4(y, times, func, parms)

  parms <- c(f=f,w0=w0,m=m)
  X <- c(x0,v0,0)
  time <- seq(t0,t1,length=N)
  h <- diff(time)[1]
  attr(X,'names') <- c("x","dxdt","forcing")
  out <- as.data.frame(rk4(X,time,dXdt,parms))

  #plot(time,F.ext(time),type="l"); dev.new()
  #print(summary(out$x)); print(F.ext(time))
  
  if (scal) x <- out$x/sd(out$x) else x <- out$x
        
  par(bty="n",xaxt="n",yaxt="n",pty="m",xpd=NA)
  if(is.null(main)) main <- "Forced damped oscillator"
  plot(time,x,type="l",lwd=4,main=main,col="grey",
       xlab="t",ylab=expression(paste(x(t),", ",F[ext](t))))
  lines(time,rep(0,length(time)),lty=2)
  xs <- 0.5*time*pi/w0; xs <- xs[xs < max(time)]
  points(xs,rep(0,length(xs)),cex=0.5)
  if (length(time)==length(F.ext(time))) {
    lines(time,F.ext(time)/sd(F.ext(time)),col="red",lty=2)

    legend(0,min(x),c("x",expression(F[ext])),lty=c(1,2),lwd=c(4,1),
                      col=c("grey","red"),bty="n")
  }  else
    legend(0,min(x),"x",lty=1,lwd=4,col="grey",bty="n")
  legend(time[round(4*N/5)],max(x),
           c(expression(m==phantom(0)),
             expression(omega[0]==phantom(0)),expression(f==phantom(0)),
             m,w0,h),ncol=2,bty="n",cex=0.8)

  text(time[round(5*N/7)],min(x),
  expression(frac(d^2*x,d*t^2) + frac(f,m)*frac(dx,dt) + omega[0]^2 == F[ext]))

  k <- m*w0^2; c <- f*m
  attr(out,'time') <- time
  attr(out,'N') <- N
  attr(out,'f') <- f
  attr(out,'w0') <- w0
  attr(out,'m') <- m
  attr(out,'h') <- h
  attr(out,'damping.ratio') <- c/(2*sqrt(m*k))
  attr(out,'LHS') <- LHS
  invisible(out)
}

resonanceTest <- function(N=1000) {
  dev.new()
  FOscillator(x0=1,v0=0,N=N,main="Unforced oscillator",
           LHS="0",f=0,m=0.1,w0=0.1,scal=FALSE) -> a
  f <- attr(a,'w0')*attr(a,'h')
  lines(cos( f*attr(a,'time') ),col="blue",lty=2,lwd=1)
  legend(0,-1.1,"analytical",lty=2,lwd=1,col="blue",bty="n")
  
  dev.new()
  FOscillator(x0=1,v0=0,main="Unforced damped oscillator",
          LHS="0",f=0.001,m=0.1,w0=0.1,scal=FALSE) -> a
  f <- attr(a,'w0')*attr(a,'h')
  damp <- exp( -attr(a,'damping.ratio')*attr(a,'time') )
  lines(cos( f*attr(a,'time') )*damp,
        col="blue",lty=2,lwd=1)
  legend(0,-1.01,"analytical",lty=2,lwd=1,col="blue",bty="n")
  stop()

  dev.new()

  FOscillator(x0=1,v0=0,main="Damped oscillator forced w. sinusoid",
          LHS="sin(2*pi*t/100)",f=0.003,m=0.1,w0=0.4)

  dev.new()          
  F.ext <- decomposeFT(N=10000)
  
  dev.new()
  FOscillator(x0=1,v0=0,main="Damped oscillator forced w. random signal",
              LHS=F.ext,f=0.3,m=0.1,w0=1)

  dev.new()
  y <- FOscillator()

  dev.new()
  spectrum(y)

  data("co2",envir=environment())
  F.ext <- log(co2)
  spectrum(F.ext,main="log|CO2|")
  dev.new()
  
  FOscillator(x0=0,v0=0,main="Damped oscillator forced w. ln|CO2|",
              LHS=F.ext,f=3,m=10,w0=1)
}
