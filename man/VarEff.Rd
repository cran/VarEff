\name{VarEff}
\alias{VarEff}

\title{
Function VarEff from the package called also VarEff
}

\description{
\code{VarEff} (Variation of effective size) is a model to estimate the evolution of effective population size with coalescent approach. The estimation is done on simulated demographies modelled by steps of constant size for which the posterior probabilities are derived using an approximation of likelihood.

The user is asked the number of different steps he wishes to use to model the variation of effective size in the past.

Here, the function VarEff is given with all parameters to solve the model \code{VarEff}. 

It produces 2 output files: summary statistics (.Theta) and the list of detailed states simulated by the MCMC chain (.Batch). 

Further, files are produced by the additional functions that allows results to be visualized.\cr

For an overview of the package, see \link{VarEff-package}. \cr

Caution: Source the MCMC package of C.J. Geyer (2009): version 0.7-3.
}

\usage{
VarEff(infile = NULL, parafile =NULL, NBLOC = NULL, JMAX = NULL, 
MODEL = NULL, MUTAT = NULL, NBAR = NULL, VARP1 = NULL, RHOCORN = NULL, 
GBAR = NULL, VARP2 = NULL, DMAXPLUS = NULL, Diagonale = NULL, 
NumberBatch = NULL, LengthBatch = NULL, SpaceBatch = NULL)
}

\arguments{
  \item{parafile}{[matrix]:Output data.
}
  \item{infile}{[matrix]:Input data.
}
  \item{NBLOC}{[integer]:Number of markers.
}
  \item{JMAX}{[integer]:Number of times when the effective size has changed.
}
  \item{MODEL}{[numeric]:Mutation model.
}
  \item{MUTAT}{[numeric]:Mutation rate.
}
  \item{NBAR}{[integer]:Global prior mean of effective size.
}
  \item{VARP1}{[integer]:Logarithmic variance of effective size.
}
  \item{RHOCORN}{[numeric]:Coefficient of correlation between effective sizes in successive steps.
}
  \item{GBAR}{[integer]:Number of generations since the origin of the population.
}
  \item{VARP2}{[integer]:Logarithmic variance of time interval.
}
  \item{Diagonale}{[numeric]:Smoothing parameter. 
}
  \item{DMAXPLUS}{[integer]:Range of allele distances to be analysed.
}
  \item{NumberBatch}{[integer]:Number of batch.
}
  \item{LengthBatch}{[integer]:Length of batch.
}
  \item{SpaceBatch}{[integer]:Space of batch
}
}

\details{
The package named \code{VarEff} works with the main function called VarEff. 

This function is the corpus of coalescent estimation of effective sizes from the time of sampling to ancestor population size.

The user can get the posterior distribution of effective population size at any time in the past, as asked in further functions.
}

\section{Summary}{
To summarize the main part of this package go thought these four section:
\enumerate{
		\item INFILE
		\item INPUT
		\item BATCH FILE
		\item OUTPUTS
		\item PLOT
		}
}

\section{1. INFILE}{
The infile is the microsatellite data of a population.
It is near of MSVAR (Beaumont 1999) infile.
Go to https://qgp.jouy.inra.fr to convert a MSVAR file in VarEff file.

In this package an exemple can be found in the directory data called InputTest.
}

\section{2. INPUT}{

- parafile (Name that you give to the job and to the output files created by the model).

- infile (Name of the data file).

- NBLOC (Number of Loci).

- JMAX (Number of times when the effective size has changed, used to generate step functions simulating the past demography. 
Ex: JMAX=2, if you think that the population took 3 different effective sizes in the past).

- MODEL (choose one mutation model in: S = Single Step Model, T = Two Phase Model, G = Geometric Model, and provide an additional coefficient (C) for T and G models).

- MUTAT (Mutation rate, assumed the same for all loci).

- NBAR (Prior value for the effective size).

- VARP1 (Variance of the prior log-distribution of effective sizes. Ex: VARP1=3 allows for searches with 20- to 40-fold relative variations of effective size).

- RHOCORN (Coefficient of correlation between effective sizes in successive intervals).

- GBAR (Number of generations since the assumed origin of the population).

- VARP2 (Variance of the prior log-distribution of time intervals during which the population is assumed of constant size).

- DMAXPLUS = DMAX+1 (DMAX is the maximal distance between alleles (number of microsatellite motifs) that is used in the estimation algorithm).

- Diagonale (A smoothing parameter to balance the observed covariance structure with a theoretical diagonal variance matrix and avoid numerical instability. Diagonale = 0.5 is a robust choice).

- NumberBatch (number of batch (nbatch) for metrop in MCMC).

- LengthBatch (length of batch (blen) for metrop in MCMC).

- SpaceBatch (space of batch (nspac) for metrop in MCMC).
}

\section{3. BATCH FILE}{
Two files are produced by the model \code{VarEff}:

-job.Theta provides summary statistics (see below).

-job.Batch which is the core of the estimate. \cr


It reports a list of demographic evolutions described by step functions. Each line includes:

Column 1: the number  i of the simulated state (from 1 to Numberbatch).

Column 2: quadratic deviation of data from the  i-th  simulated state.

Column 3: natural logarithm of the prior probability the i-th state.

Columns 4 to JMAX+4: the JMAX+1 population sizes in the i-th state.

Columns JMAX+5 to 2 JMAX+4: the JMAX times when the size changed in the i-th state.

Columns 2 JMAX + 5: value of the C parameter of the mutation model. \cr
 

Results are kept in the job.Batch file in reduced scales: 

Theta's for population sizes, products of generation numbers times mutation rate for times of size changes. 

The additional C parameter is a constant set to 0 for the Single Step Mutation Model, to a positive value for geometrical model or a negative value for the Two Phase Model.
}

\section{4. OUTPUTS}{

Firstly, the VarEff() function produces summary statistics in a file job.Theta.\cr

The first 3 lines are the same as produces by the Theta() function:

Line 1: global Theta0, Theta1 and Theta2 estimates.

Line 2: imbalance indices ln(Theta1/Theta0) and ln(Theta2/Theta0).

Line 3: expected range of Ne values, from the minimum and maximum global Theta estimates.

An the 4th line provides means and standard deviations over simulations of the quadratic deviations of data from simulated state and of natural logarithm of the prior probabilities of the simulated states. \cr
  

To obtain the distributions of effective size at a number of generations in the past, 
from the time of sampling to an ancestral time, use the function called NatSizeDist(). \cr

This function provides 2 files with the results in the Ne scale:

-job.Nstat

-job.Ndist \cr

To obtain the distributions of logarithm of effective size at a number of generations in the past, 
from the time of sampling to an ancestral time, use the function called LogSizeDist(). \cr


This function provide 2 files with the results given in the Log(Ne) scale:

-job.Lstat

-job.Ldist \cr


Format of Nstat or Lstat file:

Column  1: Time in generations (if MUTAT is not 0) or the corresponding reduced time.

Columns 2: Arithmetic Mean of Ne or Log(Ne).

Columns 3: Harmonic means of Ne (not provided for Log(Ne)).

Columns 4: Mode of Ne or Log(Ne).

Columns 5: Median of Ne or Log(Ne).

Columns 6: Quantile 5 percent of Ne or Log(Ne).

Columns 7: Quantile 95 percent of Ne or Log(Ne. \cr


Format of Ndist or Ldist file: 

Posterior densities of Ne or Log(Ne) at past times (fitted distribution using the density R function).

File with (Nbinstants+1) lines and 514 (Ndist) or 1025 (Ldist) columns.

Lines: Instants when the distribution of N(T(i-1)) was calculated.\cr


File Ndist:

Columns in line i :

Column  1 : Values of  T(i-1).
 
Columns 2 : Size of each of the intervals (=TMAX/511) in the abcsissa (Ne scale).
 
Columns 3 to 514 : Ordinates (densities of Ne at 512 points). \cr


File Ldist:

Columns in line i :

Column  1 : Values of  T(i-1).
 
Columns 2 to 513 : Abscissa ( Log(Ne) values).
 
Columns 514 to 1025 : Ordinates (densities of these Log(Ne)). \cr

}

\section{5. PLOT}{


- Two functions are provided to plot posterior distributions: \cr

1) NTdist():	Graphical summary of the posterior distribution of estimates of Log(Ne) in the past (a 2D plot), 
using the job.Batch file and the length of time during which distributions are retrieved. \cr

2) plotNdistrib():	Plots the posterior distributions of the estimates of Ne (or Log(Ne)) at a number of times in the past, 
using a job.Ndist or job.Ldist file as previously calculated with NatSizeDist() or LogSizeDist(). \cr


- Plot Harmonic Mean from the .Nstat file. In addition to the proposed functions,
 we illustrate how to plot specific results from a file produced by NatSizeDist():

dat2=read.table("job.Nstat")

x2=dat2[,1]

y2=dat2[,3]

x2=as.numeric(x2)

y2=as.numeric(y2)

maxX2=max(x2)

maxY2=max(y2)

plot(x2,y2,type='l', lwd="2",ylim=c(0,maxY2),xlab="Time T in the past (generations)",ylab="Ne(T)")

}

\value{
\item{job.Theta}{Summary statistics}
\item{job.Batch}{Core of program to produce all files .Nstat, Ndist, etc.}
}

\references{
Nikolic N., Chevalet C. 2011. Detecting past changes of effective population size. \cr

Chevalet C., Nikolic N. 2010. Distribution of coalescent times and distances between microsatellite alleles with changing effective population size. Theoretical Population Biology, 77(3): 152-163.  
}

\author{
Natacha Nikolic <documents_57@hotmail.com> and Claude Chevalet <claude.chevalet@toulouse.inra.fr>
}

\note{
More details on the model can be found on the website: \url{https://qgp.jouy.inra.fr}
This package needs mcmc package.
The example is done with few batches to test the model.
The minimum is NumberBatch =10000, LengthBatch=10, SpaceBatch=10.
}

\seealso{
Overview: \code{\link{VarEff-package}} \cr
Exemple: \code{\link{InputTest}} \cr
HelpData: \code{\link{HelpData}} 
}

\examples{
\donttest{
VarEff(infile=system.file("data/InputTest.txt", package = "VarEff"),
                          parafile = 'job',
                          NBLOC=20,
                          JMAX=3,
                          MODEL = 'S',
                          MUTAT=0.01,
                          NBAR=1000,
                          VARP1=3,
                          RHOCORN=0,
                          GBAR=5000,
                          VARP2=3,
                          DMAXPLUS=12,
                          Diagonale=0.5,
                          NumberBatch = 100,
                          LengthBatch = 10,
                          SpaceBatch = 10)
 }
  }
\keyword{ package }
\keyword{ math }
\keyword{ methods }
\keyword{ misc }
\keyword{ models }
\keyword{ survey }
