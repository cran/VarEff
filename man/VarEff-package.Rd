\name{VarEff-package}
\alias{VarEff-package}
\docType{package}
\title{
Overview: Estimation of effective sizes from present to ancestral time 
}

\description{
This package is a model called \code{VarEff} to estimate the evolution of effective population size with coalescent approach. 

The estimation is done on simulated demographies modelled by steps of constant size for which the posterior probabilities are derived using an approximation of likelihood.
}

\details{
\tabular{ll}{
Package: \tab VarEff\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-02-07\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
Depends: \tab mcmc\cr
URL: \tab URL: https://qgp.jouy.inra.fr\cr
	}
}

\section{Overview}{
This package depends on package MCMC (C.J. Geyer, 2009, version 0.7-3), 
so you have to load the library mcmc: library(mcmc). \cr

To use this package and explore the results go through four steps:

	\enumerate{
		\item Data preparation
		\item Variable input
		\item Output files
		\item Explore the results
		}
}

\section{1. Data preparation}{

The data file describes the genotypes of a population at microsatellite markers.

The model assumes alleles defined by their lengths (number of microsatellite repeats).

The format of the file is close to MSVAR infile (Beaumont 1999). To convert a MSVAR file in \code{VarEff} file go to the Web site: 
\url{https://qgp.jouy.inra.fr} \cr


Infile: 

Each markers is described by 2 lines.  

The first line gives the number of alleles (allelic classes) at the locus.

The second line gives the numbers of alleles at each corresponding length of the locus. \cr


Caution: 

You have to specify all potential alleles between those of  minimum and maximum lengths.

It means that if you have a locus with 2 types of alleles at the lengths 10 and 12 (number of repeat motifs), 
you have to mention the unobserved allele with 11 motifs. \cr

So if the alleles 10 and 12 have observed at the frequencies 24 and 6 respectively, you have to describe the locus by:

3

24 0 6  \cr

In this package the infile test is called \code{\link{InputTest}}.
}

\section{2. Variable input}{

Because this model follows a Bayesian approach, you have to give priors on effective sizes (current and ancestral) and age of the population,
specifying means and variances on the logarithmic scale. 

Estimations are concerned with reduced population sizes (on the Theta = 4 * N * u scale) 
and reduced time (product of generation times (T) and mutation rate (u)). 

The mutation rate is not estimated in this package. The package used the mutation rate as a scale parameter to recover actual census size and actual times (generation numbers) from the results.

In the case that you want estimation not in Theta (4*Ne*u) but in effective size Ne), you have to estimate previously the mutation rate (u) with existing method (Ex: MSVAR (Beaumont 1999). 

The other parameters use in the package to visualize results concern the time (generations) that you want to go back and the times you want to watch. \cr


Call the package \code{VarEff} then answer the questions: \cr

- parafile (Name that you give to the job and to the output files created by the model).

- infile (Name of the data file).

- NBLOC (Number of Loci).

- JMAX (Number of times when the effective size has changed, used to generate step functions simulating the past demography. Ex: JMAX=2, if you think that the population took 3 different effective sizes in the past).

- MODEL (choose one mutation model in: S = Single Step Model, T = Two Phase Model, G = Geometric Model, and provide an additional coefficient (C) for T and G models).

- MUTAT (Mutation rate, assumed the same for all loci).

- NBAR (Global prior mean of effective size).

- VARP1 (Variance of the prior log-distribution of effective sizes. Ex: VARP1=3 allows for searches with 20- to 40-fold relative variations of effective size).

- RHOCORN (Coefficient of correlation between effective sizes in successive intervals).

- GBAR (Number of generations since the assumed origin of the population).

- VARP2 (Variance of the prior log-distribution of time intervals during which the population is assumed of constant size).

- DMAXPLUS = DMAX+1 (DMAX is the maximal distance between alleles (number of microsatellite motifs) that is used in the estimation algorithm).

- Diagonale (A smoothing parameter to balance the observed covariance structure with a theoretical diagonal variance matrix and avoid numerical instability. Diagonale = 0.5 is a robust choice).

- NumberBatch (number of batch (nbatch) for metrop in MCMC).

- LengthBatch (length of batch (blen) for metrop in MCMC).

- SpaceBatch (space of batch (nspac) for metrop in MCMC). \cr



You can also directly give the parameters into R console. \cr


Exemple with data \code{\link{InputTest}}: 

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
                          NumberBatch = 2,
                          LengthBatch = 1,
                          SpaceBatch = 1)


}

\section{3. Output files}{

At the end of the calculations, VarEff() returns global theta values, the summaries of adjustment criteria of data to model, and the distribution of posterior probabilities, which are added to the .Theta file.\cr 


The 4 lines includes in .Theta file contain:

Line 1: global Theta0, Theta1 and Theta2 estimates.

Line 2: imbalance indices ln(Theta1/Theta0) and ln(Theta2/Theta0).

Line 3: expected range of Ne values, from the minimum and maximum global Theta estimates.

Line 4: means and standard deviations over simulations of the quadratic deviations of data from simulated state and of natural logarithm of the prior probabilities of the simulated states.\cr


The main result of VarEff() is the .Batch file, which reports a list of demographic evolutions described by step functions. Each line includes:

Column 1: the number  i of the simulated state (from 1 to Numberbatch).

Column 2: quadratic deviation of data from the  i-th  simulated state.

Column 3: natural logarithm of the prior probability the i-th state.

Columns 4 to JMAX+4: the JMAX+1 population sizes in the i-th state.

Columns JMAX+5 to 2 JMAX+4: times of size changes in the i-th state.

Columns 2 JMAX + 5: value of the C parameter of the mutation model. \cr


Results are kept in the .Batch files in reduced scales: 

Theta's for population sizes, products of generation numbers times mutation rate for times of size changes. 

The additional C parameter is set to 0 for the Single Step Mutation Model, positive for geometrical model or negative for the Two Phase Model. \cr


\code{\link{NatSizeDist}}  \cr


To obtain the distributions of effective size at a number of generations in the past, from the time of sampling to an ancestral time, use the function called NatSizeDist().\cr

This function provides 2 files with the results in the Ne scale:

-job.Nstat

-job.Ndist \cr

Format of Nstat or Lstat file

Column  1: Time in generations (if MUTAT is not 0) or the corresponding reduced time.

Columns 2: Arithmetic Mean of Ne or Log(Ne).

Columns 3: Harmonic means of Ne (not provided for Log(Ne), set to 0 in .Lstat file).

Columns 4: Mode of Ne or Log(Ne).

Columns 5: Median of Ne or Log(Ne).

Columns 6: Quantile 5 percent of Ne or Log(Ne).

Columns 7: Quantile 95 percent of Ne or Log(Ne). \cr


\code{\link{LogSizeDist}} \cr

To obtain the distributions of logarithm of effective size at a number of generations in the past, from the time of sampling to an ancestral time, use the function called LogSizeDist(). \cr

This function provides 2 files with the results in the Log(Ne) scale:

-job.Lstat

-job.Ldist \cr


Format of Ndist or Ldist file: 

Posterior densities of Ne or Log(Ne) at past times (fitted distribution using the density R function).

File with (Nbinstants+1) lines  and 514 or 1025 columns.

Lines: Instants when the distribution of N(T(i-1)) was calculated (1<i<Nbinstants+1; 0<T(i-1)<Tempsmax). \cr


File Ndist:

Columns in line i :

Column  1 : Value of  T(i-1).
 
Columns 2 : Size of each of the intervals (=TMAX/511) in the abcsissa (Ne scale).
 
Columns 3 to 514 : Ordinates (densities of  Ne at 512 points). \cr


File Ldist:

Columns in line i :

Column  1 : Value of  T(i-1).
 
Columns 2 to 513 : Abscissa ( Log(Ne) values).
 
Columns 514 to 1025 : Ordinates (densities of these Log(Ne)).

}

\author{
Natacha Nikolic <documents_57@hotmail.com> and Claude Chevalet <claude.chevalet@toulouse.inra.fr>

Maintainer: Who to complain to <documents_57@hotmail.com>
Natacha Nikolic
}

\references{
Nikolic N, Chevalet C. 2011. Detecting the evolution of coalescent effective population size. In preparation. \cr

Chevalet C & Nikolic N. 2010. Distribution of coalescent times and distances between microsatellite alleles with changing effective population size. Theoretical Population Biology, 77(3): 152-163.
}

\keyword{ package }
\keyword{ math }
\keyword{ methods }
\keyword{ misc }
\keyword{ models }
\keyword{ survey }
\keyword{ htest }

\seealso{
Summary: \code{\link{VarEff}} \cr
Exemple: \code{\link{InputTest}} \cr
HelpData: \code{\link{HelpData}} \cr
Functions to built output files: \code{\link{NatSizeDist}} \cr or \code{\link{LogSizeDist}} \cr
Functions to visualize and plot the results: \code{\link{plotNdistrib}} \cr and  \code{\link{NTdist}} \cr
Web site: \url{https://qgp.jouy.inra.fr} 
}
