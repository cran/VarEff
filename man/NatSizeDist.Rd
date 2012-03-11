\name{NatSizeDist}
\alias{NatSizeDist}
\title{
Reveal the estimates of Ne
}
\description{
Calculates estimates of effective size (Ne), at a number of times from 0 to a certain time ago (given by the user), plots and saves results on files.
}
\usage{
NatSizeDist(NameBATCH = NULL, MUTAT = NULL, TMAX = NULL, NBT = NULL)
}
\arguments{

  \item{NameBATCH}{[matrix]:File .Batch create by the function VarEff.
}
  \item{MUTAT}{[numeric]:Mutation rate.
}
  \item{TMAX}{[integer]:Length of the period for which the distributions of Ne in the past are generated.
}
  \item{NBT}{[integer]:Number of time intervals. 
  
Example: If TMAX=1000 generations, and NBT=100, the estimates are calculated every 10th generation until 1000 generations ago.
}
}
\details{
To obtain the distribution of effective size (Ne) at a number of generations in the past, run the function  NatSizeDist(). \cr

The function makes use of the job.Batch file created by VarEff(). Results can be shown in natural census values of population sizes, or in reduced theta scale if mutation rate is set to 0 in these functions.\cr


Nstat file:

Column  1: Time in generations (if MUTAT is not 0) or the corresponding reduced time.

Columns 2: Arithmetic Mean of Ne.

Columns 3: Harmonic means of Ne.

Columns 4: Mode of Ne.

Columns 5: Median of Ne.

Columns 6: Quantile 5 percent of Ne.

Columns 7: Quantile 95 percent of Ne. \cr


Ndist file: 

Posterior densities of Ne or Log(Ne) at past times (fitted distribution using the density R function).

File with  (Nbinstants+1) lines and 514 columns.

Lines: Instants when the distribution of N(T(i-1)) was calculated: 1<i<NBT+1; 0<T(i-1)<TMAX.
TMAX corresponds to the length of the period for which the distributions of Log(Ne) in the past are generated).
NBT corresponds to the number of time intervals. 
Example: If TMAX=1000 generations, and NBT=100, estimates are calculated every 10th generation until 1000 generations ago).


Columns in line i :

Column  1 : Values of  T(i-1).
 
Columns 2 : Size of each of the intervals (=TMAX/511) in the abcsissa (Ne scale).
 
Columns 3 to 514 : Ordinates (densities of  Ne at 512 points).\cr

}
\value{
\item{job.Nstat}{Summary statistics of posterior effective population size. Matrix(ncol=7,nrow=NBT+1):
 Time in generations (or reduced time (Theta) if MUTAT was set to 0), Arithmetic mean of Ne, Harmonic mean of Ne, Mode of Ne, Median of Ne, Quantile 5 percent of Ne, Quantile 95 percent of Ne.}
\item{job.Ndist}{Containing the detailed densities of the posterior distributions of Ne. }
}
\references{
Nikolic N., Chevalet C. 2011. Detecting past changes of effective population size.\cr

Chevalet C., Nikolic N. 2010. Distribution of coalescent times and distances between microsatellite alleles with changing effective population size. Theoretical Population Biology, 77(3): 152-163.  
}
\author{
Natacha Nikolic <documents_57@hotmail.com> and Claude Chevalet <claude.chevalet@toulouse.inra.fr>
}
\note{
More details on the model can be found on the website: \url{https://qgp.jouy.inra.fr}
}
\seealso{
Overview: \code{\link{VarEff-package}} \cr
Exemple: \code{\link{InputTest}} \cr
HelpData: \code{\link{HelpData}} 
}
\examples{ 
NatSizeDist (NameBATCH = system.file("extdata/job.Batch", package = "VarEff") 
, MUTAT = 0.01, TMAX = 200, NBT = 10)
 }
\keyword{ models }

