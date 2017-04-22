\name{Demonstrations}
\alias{decomposeFT}
\alias{ENSO.example}
\alias{resonance}
\alias{resonanceTest}
\alias{diffdemo}
\alias{diff12demo}
\alias{testLTP}
\alias{demoConf}
\alias{demoRange}
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


\code{ENSO.example()} shows how the strategy adopbed in Scafetta (2012)
fails when applied to ENSO (the NINO3.4 index). \code{resonance()} shows
that a system with resonance will pick up the resonant frequency from
any noisy forcing - this is analogous to the whine from the wind blowing
aroundcorners, and how musical pipes/trumpets work (the function
\code{resonanceTest} provides some
demonstrations). A regression to harmoincs uses \code{hfit} fits harmonics according to
eq. (3)in Scafetta (2011). \code{decomposeFT()} 
shows how any curve can be represented as a sum of harmonics - Fourier
series. REference: Scafetta,N.,Testing an astronomically based
decadal-scale empirical harmonic climate model versus the
IPCC (2007) general...., Journal of Atmospheric and Solar-Terrestrial
Physics (2011) doi:10.1016/j.jastp.2011.12.005.
Also on \url{http://arxiv.org/abs/1201.1301} 

\code{Walker.test()} is the function for the Walker test to test
the significance when many tests are made (problem of multiplicity/field
significance). 

\code{testLTP} provides a demonstration/test showing that sophisticated
trend fits will also be affected by the trend in data. This function
compares the auto-correlation functions (\code{showACF}) of time series
from a global climate model: both global mean and interpolated to
Svalbard, where one uses constant boundary conditions (e.g. no trend)
and the other includes the 20th century greenhouse gas forcings
(non-zero trends). This analysis is applied to the annual mean values
(\code{am}). These tests apply to the papers
\url{http://dx.doi.org/10.1029/2005GL024476} and
\url{http://dx.doi.org/10.1029/2012GL054244}, where sophisticated trend
models were used for hypothesis testing. It is important to keep in mind
that these trend models then take the trends to be part of the noise,
and hence are unsuitable for testing whether the trends are
statistical significant.

\code{demoConf} and \code{demoRange} demonstrate how the confidence
interval gets narrower for the estimate of the mean value as the sample
size increases. This is done by generating a large set of random numbers
("original data") with known mean and standard deviation (indicator of the data range),
and then a set of random subsets of the data. The mean values and standard deviation
are computed for each of these subsets. As the sample size increases,
the estimates converge towards the prescribed ("true") value. However,
the range of the original master sample is not affected by the
subsampling, and the range of the subsamples converge towards that of
the original data as their sample size approach that of the original
master sample. The data range and the confidence interval for the
estimate of the mean were mixed up in \url{http://dx.doi.org/10.1002/joc.1651}.


}

\usage{
diffdemo(x=0.7*cos(seq(0,10*pi,length=1000))+0.4*rnorm(1000),
                     y=0.9*cos(seq(0,10*pi,length=1000))+0.3*rnorm(1000))
diff12demo(x=0.5*cos(seq(0,10*pi,length=1200))+rnorm(1200),
                     y=0.7*cos(seq(0,10*pi,length=1200))+rnorm(1200),
                     wfl=12)
decomposeFT(N=1000)
ENSO.example(interval=1980:1989)
resonance(x0=1,v0=0,t0=0,t1=1000,N=1000,
          LHS=cos(2*pi*(1:1000)/75),
          f=0.001,m=0.1,w0=0.07)
resonanceTest(N=1000)
testLTP()
demoConf(m=3,s=2,N=10000,n=c(10,50,100,200,500,1000,5000))
demoRange(m=3,s=1) 
}
\arguments{
    \item{interval}{Time interval to analyse}
    \item{x0}{Initial condition for x of oscillator; or time series for control}
    \item{v0}{Initial condition for v = dx/td}
    \item{t0}{Initial time index}
    \item{t1}{Finalt time index}
    \item{LHS}{Left hand side of the equation}
    \item{f}{friction term for damped oscillator}
    \item{m}{inertia term for oscillator; mean value}
    \item{w0}{frequency term}
    \item{x}{curve to which harmonics are fitted or time series - a
      vector of numbers; or a time series}
    \item{y}{time series - a vector of numbers}
    \item{N}{length of time series; Size of data sample}
    \item{wfl}{window filter length}
    \item{s}{standard deviation}
    \item{n}{size of subsets}
}
\value{A table or lists containing the relevant data.}
\author{R.E. Benestad}
\examples{
\dontrun{

# Demonstrate the limitations of the diff-operator for two noisy signals
# (red and black in the upper panel respectively) with similar long-term
# harmonics. The lower panel shows the lagged correlation for the
# diff-operated series. 
diffdemo()
diff12demo()

# Demonstration: show that a noise consists of many Fourier components/harmonics
decomposeFT()


# Test the assumption about on good cycle-fit for a curve-fit toanother
# cycle, as done in Scafetta (2011)     
ENSO.example()
# Test the Runge-Kutta integration of a forced damped oscillator to test
# the claim about resonance made by Scafetta.
resonanceTest()



}
}
\keyword{manip}


