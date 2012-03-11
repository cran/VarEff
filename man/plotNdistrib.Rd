\name{plotNdistrib}
\alias{plotNdistrib}
\title{
Plots the posterior distributions of Ne or Log(Ne))
}
\description{
Plots the posterior distributions of the estimates of Ne (or Log(Ne)) at a number of times in the past.
}
\usage{
plotNdistrib(infile = NULL, nbcases=NULL)
}
\arguments{
  \item{infile}{[matrix]:job.Ndist created by the function NatSizeDist
or
job.Ldist created by the function LogSizeDist.
}
\item{nbcases}{[integer]: number of times (suggestion: <=5) you wish to plot the distribution.
}
}
\details{
The function plotNdistrib() makes use of files .Ndist or .Ldist previously built by NatSizeDist() or LogSizeDist(). 

It allows the user to exhibit the density of the posterior distribution of Ne (or of Log(Ne)) at several times in the past. 

Compared to NTdist() that gives a rough but global view of these densities, plotNdistrib() gives a precise view of these densities at a small number of times. 

The global 2D plot given by NTdist() may help choosing the times when plotNdistrib() is used. \cr

Caution! In the example below, we used nbcases=0 to get out. It allowed to built the package because the function is interactive.

If you want to use it just call the function 

plotNdistrib () 

}
\value{
\item{job_Ndist}{Figure of the posterior distributions of Ne at R console.}
\item{job_Ldist}{Figure of the posterior distributions of Log(Ne) at R console.}
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
}
\seealso{
Overview: \code{\link{VarEff-package}} \cr
Exemple: \code{\link{InputTest}} \cr
HelpData: \code{\link{HelpData}} 
}
\examples{
\donttest{
plotNdistrib (NameBATCH = system.file("extdata/job.Ndist", package = "VarEff") 
, nbcases=0)
plotNdistrib (infile = "../VarEff/extdata/job.Ndist",nbcases=0) 
}
}
\keyword{ models }
