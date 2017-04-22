\name{Replication}
\alias{forskning.no}
\alias{Humlum.et.al.2011}
\alias{LoehleScafetta2011}
\alias{Scafetta2010}
\alias{Scafetta2012}
\alias{Scafetta.tab1}
\alias{Solheim.et.al.2011}
\alias{do.vardo}
\alias{Scafetta2006}
\alias{check.table1}
\alias{MM2004}
\alias{DCPS2007}
\alias{paleaoproxy}
\alias{FL1991}
\alias{Humlum.et.al.2013}
\title{Demonstrations which debunk some methods and analytical set-ups.}
\description{
This set of function provide demonstrations showing why the choice of
methods and analytrical set-ups in a number of papers are
inappropriate. This computer code is also used to do the analysis on
which the paper 'Bad science papers' by Benestad et al. is based.

Some of these analyses are indeed very trivial, but carried out
nevertheless. The objective with this R-package is partly to show that
code sharing and open-source can be an effective means of resolving
differences. Methods and analytical set-up should be tested with
surrogate data for which the anwers are known a priori (a kind of method
calibration and evaluation). The spirit of this is very much like the
replication carried out by Benestad and Schmidt (2009)
\url{http://pubs.giss.nasa.gov/abs/be02100q.html}

\code{Humlum.et.al.2011} exmines the results of Humlum et al. (2011)
\url{http://www.sciencedirect.com/science/article/pii/S0921818111001457}
and extends their analysis to the data that they cut off.

\code{Humlum.et.al.2012} exmines the results of Humlum et al. (2012)
\url{http://www.sciencedirect.com/science/article/pii/S0921818112001658}.
The URLs for the data are given by \code{URLs4HSS2012}

\code{LoehleScafetta2011}. Reference: Loehle and Scafetta (2011) \url{http://benthamscience.com/open/toascj/articles/V005/74TOASCJ.htm}

\code{Scafetta2010} Some replication and some demonstrations showing
how different tempertures change results, and how different noise
processes may seem to contain long-term cycles. 


\code{Solheim.et.al.2011} replicates the analysis carried out by
solheim et al. (2011), and shows how their bootstrapping is biased, and
that their results are not really significant after
all. \code{\link{Walker.test}} is the function for the Walker test to test
the significance when many tests are made (problem of multiplicity/field
significance). \code{check.table1} and \code{obs2tab1} are
supporting functions for this replication. \code{do.vardo} repeats the
analysis for Vardo. Reference: Solheim et al (2011)
\url{http://arxiv.org/pdf/1112.3256}. \code{DJF} extracts the
December-February mean.

\code{Scafetta2006} is a modification of the R-scripts used in
Benestad & Schmidt (2008) replicating and testing the analysis in the
papers Scafetta (2006).

\code{MM2004} The R-script used to carry out the analysis of Benestad
(2004) \url{http://www.int-res.com/articles/cr2004/27/c027p171.pdf}

\code{DCPS2007} shows how the conficende intervals used in Douglas
et al. (2007) excludes many of the GCM results, on which it was
based. Reference: 
Douglas et al. (2007) \url{http://www.pas.rochester.edu/~douglass/papers/Published\%20JOC1651.pdf}

\code{paleaoproxy} examines some assertions made by Veizer (2005):
Comparison between cosmogenic Be-10 isotope proxies, CO2, and
temperature from the Vostoc ice core. \url{http://www.realclimate.org/index.php/archives/2005/05/on-veizers-celestial-climate-driver/}

\code{forskning.no} replicates a figure (again, by Solheim et al.) shown in
\url{http://www.forskning.no/artikler/2012/mars/316178} that claims that
the warming has stopped. The figure gives a wrong impression, as it's
shown on monthly scales where the monthly anomalies swamp the signal. By
plottingthe monthly mean values instead, the warming trends becomes
visible. 



These papers are also discussed on \url{http://www.RealClimate.org}: \url{http://www.realclimate.org/index.php/archives/2011/12/curve-fitting-and-natural-cycles-the-best-part/} and \url{http://www.realclimate.org/index.php/archives/2011/08/an-exercise-about-meaningful-numbers-examples-from-celestial-attribution-studies/}
}

\usage{
Humlum.et.al.2011()
Humlum.et.al.2013(wfl=12,forcing=FALSE,HadCRUT4=FALSE,HadSST3=FALSE)
forskning.no(uah=TRUE,rss=TRUE,giss=TRUE,ncdc=TRUE,hadcrut3=TRUE,
             lwd=1,ylim=c(-0.4,1),xlim=c(1995,2012),
             base.period=1981:2010,type="s")
LoehleScafetta2011()
Scafetta2010()
Scafetta2012()
Scafetta.tab1()
check.table1(ssh2011.tab1=NULL)
Solheim.et.al.2011(obs=NULL,ssh2011.tab1=NULL,N.tests=30000)
do.vardo()
Scafetta2006(GISS.temp = TRUE, do.MonteCarlo = TRUE, test.bp = TRUE, 
    lag = 0, stepwise = TRUE, interval = 1958:2000, same.interval = TRUE, 
    all.data = FALSE, SW06.coefs.only = FALSE, wavelets.only = FALSE, 
    bivariate = TRUE, figures = TRUE, tables = TRUE, wavelet = TRUE, 
    boundary = "reflection")
MM2004()
DCPS2007()
paleaoproxy()
FL1991(dmi=FALSE,
   url="http://web.dmi.dk/fsweb/solarterrestrial/sunclimate/SCL.txt")
}
\arguments{
    \item{obs}{observations: station data.}
    \item{GISS.temp}{TRUE: use GISTEM, otherwise hadCRUT3}
    \item{test.bp}{pand-pas test}
    \item{lag}{lag between forcing and temperature response (years)}
    \item{stepwise}{TRUE: step-wise regression}
    \item{interval}{Time interval to analyse}
    \item{same.interval}{TRUE: same.interval as SW2006a}
    \item{all.data}{TRUE: use all the data data}
    \item{SW06.coefs.only}{TRUE: only estimate the coefficients in SW2006}
    \item{wavelets.only}{TRUE: only do the wavelet}
    \item{bivariate}{TRUE: only include soalr forcing and GHG in regression}
    \item{figures}{TRUE: plot figures}
    \item{tables}{TRUE: print tables}
    \item{ssh2011.tab1}{Table as in \code{ssh2011.tab1} for defining
      solar cycle epochs, correponding to Table 1 in Solheim et al. (2011).}
    \item{wavelet}{TRUE: do wav}
    \item{boundary}{Boundary for the wavelet: mra}
    \item{N.tests}{Number of Monte-Carlo simulations}
    \item{do.MonteCarlo}{TRUE: Carry out Monte-Carlo simulations}
    \item{url}{URL address of FL(1991) solar cycle length data.}
    \item{dmi}{Flag: if TRUE read from DMI web site}
    \item{uah}{Univ. Alabama Huntsville}
    \item{rss}{REmote Sensing Systems}
    \item{giss}{NASA/GISS}
    \item{ncdc}{NOAA NCDC}
    \item{hadcrut3}{Hadley Centre and Climate Research Unit, U.K.}
    \item{lwd}{see \code{\link{plot}}}
    \item{ylim}{see \code{\link{plot}}}
    \item{xlim}{see \code{\link{plot}}}
    \item{base.period}{Reference for anomalies.}
    \item{type}{see \code{\link{plot}}}
    \item{wfl}{window filter length}
    \item{forcing}{TRUE: use log(CO2) rather than CO2}
    \item{HadCRUT4}{Use HadCRUT4}
    \item{HadSST3}{Use HadSST3}
    \item{s}{standard deviation}
}
\value{A table or lists containing the relevant data.}
\author{R.E. Benestad}
\examples{
\dontrun{
# Produce stretched-out temperature graph to hide the trend: Demonstrate
# how a stretched x-axis for plots of temperature, showing
# monthly values hide the long-term trend readily visible in the annual
# mean values.
forskning.no()


# Replicat the results from Humlum et al. (2011)
Humlum.et.al.2011()


# Replicat the results from Humlum et al. (2012)
Humlum.et.al.2012()
# Demonstrate the limitations of the diff-operator for two noisy signals
# (red and black in the upper panel respectively) with similar long-term
# harmonics. The lower panel shows the lagged correlation for the
# diff-operated series. 
diffdemo()
diff12demo()


# Replication of the analysis of Solheim et al. (2011)
Solheim.et.al.2011()


# Replicate the results from Loehle and Scafetta (2011)
LoehleScafetta2011()
# Demonstration: show that a noise consists of many Fourier components/harmonics
decomposeFT()


# Replicate the analysis done by Benestad & Schmidt (2009) repeating the
# work of Scafetta (2005, 2006a, 2006b)
Scafetta2006()


# Replicate the results from Scafetta (2010)
Scafetta2010()


# Replicate the results from Scafetta (2011)
Scafetta2011()
# Plot the results shown in Table 1 in Scafetta (2011) and the results
# obtained with an unbiased method.
Scafetta.tab1()
# Test the assumption about on good cycle-fit for a curve-fit toanother
# cycle, as done in Scafetta (2011)     
ENSO.example()
# Test the Runge-Kutta integration of a forced damped oscillator to test
# the claim about resonance made by Scafetta.
resonanceTest()


# Replicate the results from McKitrick & Michaels (2004)
MM2004()


# Replicate the results from Douglass et al. (2007)
Douglass2007()


# Replicate the results from Veizer (2005)
paleaoproxy()


# Replicate the solar cycle length from Friis-Christensen and Lassen (1991)
FL1991()
}
}
\keyword{manip}


