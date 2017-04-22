require(zoo)


stand <- function(x) {
  y <- (x - mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
  return(y)
}

#HC4nh <- read.table('http://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/time_series/HadCRUT.4.2.0.0.monthly_nh.txt')
HC4nh <- read.table('http://www.metoffice.gov.uk/hadobs/hadcrut4/data/4.2.0.0/time_series/HadCRUT.4.2.0.0.monthly_nh.txt')
HC4sh <- read.table('http://www.metoffice.gov.uk/hadobs/hadcrut4/data/4.2.0.0/time_series/HadCRUT.4.2.0.0.monthly_sh.txt')
HC430 <- read.table('http://www.metoffice.gov.uk/hadobs/hadcrut4/data/4.2.0.0/time_series/HadCRUT.4.2.0.0.monthly_30S_30N.txt')

yr <- substr(HC4nh$V1,1,4)
mo <- substr(HC4nh$V1,6,7)
dates <- as.Date(paste(yr,mo,'01',sep="-"))

# Create zoo objects for satandardised anomalies:
NH <- zoo(stand(HC4nh$V2),order.by=dates)
SH <- zoo(stand(HC4sh$V2),order.by=dates)
Tr <- zoo(stand(HC430$V2),order.by=dates)

# detrend:
dev.new(width=5,height=9)
par(mfrow=c(3,2),bty="n")
gl.tr <- 0.5*(NH + SH)
plot(gl.tr,main="NH + SH",ylab="",xlab="")
acf(coredata(gl.tr),lag.max = 240,main="Global")

gl.res <- NH - SH
plot(gl.res,main="NH - SH",ylab="",xlab="")
acf(coredata(gl.res),lag.max = 240,main="Hemisphere diff.")

tr.res <- Tr - gl.tr
plot(tr.res,main="TR - 0.5(NH + SH)",ylab="",xlab="")
acf(coredata(tr.res),lag.max = 240,main="Tropical diff.")
dev2bitmap('ltpforc.png')



