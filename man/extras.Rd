\name{Extra functions}
\alias{hfit}
\alias{diff12}
\alias{showACF}
\alias{am}
\alias{DJF}
\alias{shake}
\alias{obs2tab1}
\alias{Walker.test}
\alias{instring}
\alias{dT}
\alias{ma.filt}
\alias{gauss.filt}
\alias{reverse}
\title{Demonstrations which debunk some methods and analytical set-ups.}
\description{
This set of function provide demonstrations showing how different
methods or models work. This computer code is also used to do the analysis on
which the paper 'Agnotology: learning from mistakes' by Benestad et al. is based.

Some of these analyses are indeed very trivial, but carried out
nevertheless. The objective with this R-package is partly to show that
code sharing and open-source can be an effective means of resolving
differences. Methods and analytical set-up should be tested with
surrogate data for which the anwers are known a priori (a kind of method
calibration and evaluation). The spirit of this is very much like the
replication carried out by Benestad and Schmidt (2009)
\url{http://pubs.giss.nasa.gov/abs/be02100q.html}


A regression to harmoincs uses \code{hfit} fits harmonics according to
eq. (3)in Scafetta (2011). 

\code{Walker.test()} is the function for the Walker test to test
the significance when many tests are made (problem of multiplicity/field
significance). 

A comparison between auto-correlation functions (\code{showACF}).

The annual mean values are computed with (\code{am}). 


}

\usage{
hfit(x,t,T1=60,T2=20)
diff12(x,wfl=NULL)
showACF(x0,x1,x2=NULL)
Walker.test(N,alpha=0.05)
am(x)
DJF(obs=NULL,ssh2011.tab1=NULL,yr1c24=2009)
shake(x)
instring(c,target,case.match=TRUE)
dT(y,maxhar=NULL,plot=FALSE,chk.conf=1)
ma.filt(x,n)
gauss.filt(x,n)
reverse(x)
}
\arguments{
  \item{obs}{observations: station data.}
  \item{ssh2011.tab1}{Table as in \code{ssh2011.tab1} for defining
      solar cycle epochs, correponding to Table 1 in Solheim et al. (2011).}
    \item{x0}{time series for control}
    \item{x1}{time series from model results}
    \item{x2}{time series from observations}
    \item{x}{curve to which harmonics are fitted or time series - a
      vector of numbers; or a time series}
    \item{t}{times in year}
    \item{T1}{periodicity for first harmonic}
    \item{T2}{periodicity for second harmonic}
    \item{N}{length of time series; Size of data sample}
    \item{wfl}{window filter length}
    \item{alpha}{Level of statistical significance}
    \item{yr1c24}{Year 1 of solar cycle 24}
    \item{c}{character}
    \item{target}{string in which to look for character}
    \item{case.match}{TRUE: case sensitive}
    \item{y}{time series - a vector}
    \item{maxhar}{Maximum number of harmonics to use in the fit/Fourier truncation.}
    \item{plot}{TRUE: plot}
    \item{chk.conf}{Include coefficient associated with a given confidence level}
    \item{n}{length of time window.}
}
\value{A table or lists containing the relevant data.}
\author{R.E. Benestad}
\examples{
# Plot the winter temperature for Svalbard:
DJF()

# Show the autocorrelation functions to see the effect of LTP:
data(echam5.0)
data(echam5.1)
showACF(echam5.0,echam5.1)


}
\keyword{manip}


