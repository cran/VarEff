\name{LogSizeDist}
\alias{LogSizeDist}
\title{
Reveal the estimates of Log(Ne)
}
\description{
Calculates estimates of effective size in Logarithm (Log(Ne)), at a number of times from 0 to a certain time ago (given by the user), plots and saves results on files.
}
\usage{
LogSizeDist(NameBATCH = NULL, MUTAT = NULL, TMAX = NULL, NBT = NULL)
}
\arguments{
  \item{NameBATCH}{[matrix]:File .Batch created by the function VarEff.
}
  \item{MUTAT}{[numeric]:Mutation rate.
}
  \item{TMAX}{[integer]:Length of the period for which the distributions of Log(Ne) in the past are generated.
}
  \item{NBT}{[integer]:Number of time intervals.
 
Example: If TMAX=1000 generations, and NBT=100, the estimates are calculated every 10th generation until 1000 generations ago.
}
}
\details{
To obtain the distribution of Log(Ne) at a number of generations in the past, 
run the function LogSizeDist(). \cr

The function uses the job.Batch file created by VarEff(). 
Results can be shown in natural census values of population sizes, or in reduced theta scale if mutation rate is set to 0 in this function. \cr


Lstat file:

Column  1: Time in generations (if MUTAT is not 0) or the corresponding reduced time.

Columns 2: Arithmetic Mean of Log(Ne).

Columns 3: None (value set to 0).

Columns 4: Mode of Log(Ne).

Columns 5: Median of Log(Ne).

Columns 6: Quantile 5 percent of Log(Ne).

Columns 7: Quantile 95 percent of Log(Ne).\cr


Ldist file: 

Posterior densities of Log(Ne) at past times (fitted distribution using the density R function).

File with (Nbinstants+1) lines and 514 or 1025 columns.

Lines: Instants when the distribution of N(T(i-1)) was calculated: 1<i<NBT+1; 0<T(i-1)<TMAX.

TMAX corresponds to the length of the period for which the distributions of Log(Ne) in the past are generated)

NBT corresponds to the number of time intervals. Example: If TMAX=1000 generations, and NBT=100, 
estimates are calculated every 10th generation until 1000 generations ago)

Columns in line i :

Column  1 : Values of  T(i-1) 

Columns 2 to 513 : Abscissa (Log(Ne) values) 

Columns 514 to 1025 : Ordinates (densities of these Log(Ne)) \cr

}

\value{
\item{job.Lstat}{Summary statistics of posterior effective population size in logarithm (Log(Ne)). 
Matrix(ncol=7,nrow=NBT+1): Time in generations or in in reduced scale (Theta) (if MUTAT set to 0), 
Arithmetic mean of Log(Ne), Mode of Log(Ne), Median of Log(Ne), 
Quantile 5 percent of Log(Ne), Quantile 95 percent of Log(Ne).}
\item{job.Ldist}{Densities of the posterior distributions of Ne 

(Matrix(ncol=1025,nrow=NBT+1).
}
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
\donttest{  
LogSizeDist (NameBATCH = system.file("extdata/job.Batch", package = "VarEff") 
, MUTAT = 0.01, TMAX = 200, NBT = 10)

}
}
\keyword{ models }
