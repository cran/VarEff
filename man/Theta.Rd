\name{Theta}
\alias{Theta}
\title{
Preliminary analysis of data to choose Theta priors and DMAX from microsatellites data
}
\description{
This function allows to get Theta (4*Ne*u) using 3 estimators correlated to the mean of
Present (Theta_0), Intermediate (Theta_1) and Ancestral (Theta_2) population size.\cr

The results are written to a .Theta file. \cr

It also provides Imbalance indices ln(Theta_1/Theta_0) and ln(Theta_2/Theta_0),
the minimum of census Ne (Ne0=Theta_0/4*u) and the maximum Ne (Ne_2=Theta_2/4*u).\cr



Values of .Theta file:

- Theta_0, Theta_1, Theta_2.

- Imbalance indices ln(Theta_1/Theta_0) and ln(Theta_2/Theta_0).

- MinNe, MaxNe.

- Means and standard deviations of the quadratic deviations of data. \cr



These results provide an overview of priors necessary for VarEff function.

The programm also allows to built a file .R to process the model VarEff with the function VarEff()
}
\usage{
Theta(infile = NULL, parafile = NULL, NBLOC = NULL, MUTAT = NULL, 
JMAX = NULL, MODEL = NULL, NBAR = NULL, VARP1 = NULL, RHOCORN = NULL, 
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
  \item{JMAX}{[integer]:Number of time that the effective size has changed.
}
  \item{MODEL}{[numeric]:Mutation model.
}
  \item{MUTAT}{[numeric]:Mutation rate.
}
  \item{NBAR}{[integer]:Global mean of effective size.
}
  \item{VARP1}{[integer]:Logarithmic variance of effective size.
}
  \item{RHOCORN}{[numeric]:Coefficient of correlation between successive effective size.
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
  \item{SpaceBatch}{[integer]:Space of batch.
}
}
\details{
It provides three global theta: Theta0 (linked to current effective size N0, Theta=4*N0*u), 
Theta1 (linked to intermediate size), Theta2 (linked to ancestral size), 
imbalance indices, and a range of population size.\cr

It also builds a script to run VarEff (job.R) and a plot showing the distribution of allele distances in the sample.
}
\value{
\item{job.Theta}{File with the first 3 lines: the 3 theta values, the 2 imbalances indices, the minimal and maximal suggested population sizes.}
\item{job.R}{A file containing a script to run the VarEff function.}
\item{R.plot}{Figure showing the distribution of allele distances in the sample.}
}
\references{
Nikolic N., Chevalet C. 2011. Detecting past changes of effective population size. \cr

Chevalet C., Nikolic N. 2010. Distribution of coalescent times and distances between microsatellite alleles with changing effective population size. 
Theoretical Population Biology, 77(3): 152-163.  
}
\author{
Natacha Nikolic <documents_57@hotmail.com> and Claude Chevalet <claude.chevalet@toulouse.inra.fr>
}
\note{
More details on the model can be found on the website: \url{https://qgp.jouy.inra.fr}
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
Theta (infile=system.file("data/InputTest.txt", package = "VarEff"),
                          parafile = 'job',
                          NBLOC=20,
                          MUTAT=0.01,
                          JMAX=3,
                          MODEL = 'S',
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
\keyword{ priors }