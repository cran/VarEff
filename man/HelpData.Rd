\name{HelpData}
\alias{HelpData}
\docType{data}
\title{
Help: Infile to test the model
}

\description{
The infile called InputTest provides an example used to give an overview of possible results.
  
It is a simulation of a population that underwent a huge bottleneck 200 generations ago. The past effective population was  1000 and the present size is 100.  
}

\usage{data(InputTest)}
\format{
The format is:\cr


 num [1:240] 9 4 12 2 0 14 21 15 3 9 ...
  }

\details{
The data file \code{InputTest.txt} is an example which can be used directly in the R console.
}
\section{How to use the example}{
Read the parts \code{\link{VarEff-package}} \cr
and
\code{\link{VarEff}} \cr
}

\section{Meaning of input variables}{
The variables to run one of the two functions (Theta or VarEff) correspond to the input arguments: \cr

- parafile (Name that you give to the job and to the output files created by the model).

- infile (Name of the data file).

- NBLOC (Number of Loci).

- JMAX (Number of times when the effective size has changed, used to generate step functions simulating the past demography. 

Ex: JMAX=2, if you think that the population took 3 different effective sizes in the past).

- MODEL (choose a mutation model in: S = Single Step Model, T = Two Phase Model, G = Geometric Model, and provide an additional coefficient (C) for T and G models).

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

- SpaceBatch (space of batch (nspac) for metrop in MCMC). \cr  


The input example corresponds to a population genotyped with 20 microsatellite markers (NBLOC): \cr
 
We supposed 3 huge events (JMAX) which have affected this population. JMAX has no impact on the present effective size but it can affect the distribution of past effective size. Hence check different JMAX (from 1 to 10) and keep the one in which the mean, median and mode are the nearest.

Generally the mutation rate (MUTAT) is around 0.01 to 0.0001. 0.01 was used in the simulation which built the data. 

Concerning the other priors, we set the priors of actual effective size (NBAR) at 1000 with a large variance (VARP1=3). 

We also assumed no correlation (RHOCORN=0) between the successive effective sizes from present to ancestral time. 

We supposed a time since the ancestral population (GBAR) of 5000 generations with a large variance (VARP2=3).

Even if the range of alleles length differences (number of repeat motifs) in InputTest is equal at 18, we restricted DMAXPLUS to 12 because it is representative of more than 95 percent of data (above the red line in graph Fk).\cr 

Here, the metrop (method of Monte Carlo Markov Chain) parameters are 1000 for the number of batches (NumberBatch), 10 for the length of batches (LengthBatch) and 10 for the space between batches (SpaceBatch) for a quick run.
Recommended parameters are 10000 (NumberBatch), 10 (LengthBatch) and 10 (SpaceBatch) i.e. 1 million of iterations. \cr                       

From the core file (.Batch) use the function \code{\link{NatSizeDist}} \cr or \code{\link{LogSizeDist}} \cr to estimate the effective size (Ne) (or Log(Ne)) at a number of times from 0 to a certain time ago (given by the user).


To visualize the results use the functions:
\code{\link{NTdist}} \cr and \code{\link{plotNdistrib}} \cr
  
}

\source{
This data is a simulated set obtained by the model DemoDivMs from Nikolic and Chevalet \url{https://qgp.jouy.inra.fr}.
}

\references{
Nikolic N., Butler J., Bagliniere JL., Laughton R., McMyn I.A.G, Chevalet C. 2009. An examination of genetic diversity and effective population size in Atlantic salmon populations, and applications for conservation and management. Genetics Research. 91: 1-18.  
}

\keyword{datasets}
