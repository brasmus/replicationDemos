
require(replicationDemos)

am <- function(x) {
  if (class(x)=="zoo") {
    y <-  aggregate(x, by= as.Date(cut(time(x), "year")), mean, na.rm = FALSE)
  } else if (class(x)=="data.frame") 
    y <- rowMeans(svalbard[,3:14])
  
  y <- y[is.finite(y)]
  y <- y - mean(y,na.rm=TRUE)
  invisible(y)
}

showACF <- function(x0,x1,x2=NULL) {
  n <- min(length(x0),length(x1))
  x0 <- x0[1:n]; x1 <- x1[1:n]
  par(bty="n",mfcol=c(2,1))
  plot(coredata(am(x1)),type="l",lwd=5,xlab="years",ylab="T(2m)")
  lines(coredata(am(x0)),col="grey",lwd=5)
  if (!is.null(x2)) lines(coredata(am(x2)),col="red",lwd=2,lty=2)

  ar.0 <- acf(coredata(am(x0)),plot=FALSE,na.action=na.pass)
  ar.1 <- acf(coredata(am(x1)),plot=FALSE,na.action=na.pass)
  if (!is.null(x2)) ar.2<- acf(coredata(am(x2)),plot=FALSE)
  plot(ar.0$lag,ar.0$acf,type="l",col="grey",lwd=5, log="x",
       main=attr(x0,'description'),xlab="lag (years)",ylab="ACF")
  lines(ar.1$lag,ar.1$acf,lwd=5)
  if (!is.null(x2)) lines(ar.2$lag,ar.2$acf,col="red",lwd=2,lty=2)
}

testLTP <- function() {
  
# Need to add the observed series for Svalbard, etc.

  data("echam5.0",envir=environment())
  data("echam5.1",envir=environment())
  data("svalbard.0",envir=environment())
  data("svalbard.1",envir=environment())
  data("svalbard",envir=environment())
  
  showACF(echam5.0,echam5.1)

  dev.new(width=4,height=8)
  showACF(svalbard.0,svalbard.1,svalbard)
}





