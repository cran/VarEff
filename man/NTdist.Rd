\name{NTdist}
\alias{NTdist}
\title{
Graphical summary
}
\description{
Graphical summary of the posterior distribution of estimates of Log(Ne) in the past (a 2D plot).
}
\usage{
NTdist(NameBATCH = NULL, MUTAT = NULL, TMAX = NULL)
}
\arguments{
  \item{NameBATCH}{[matrix]: File .Batch create by the function VarEff.
}
  \item{MUTAT}{[numeric]:Mutation rate.
}
  \item{TMAX}{[integer]:Length of the period for which the distributions of Log(Ne) in the past are generated.
}
}
\details{
The function NTdist() uses the job.Batch file and makes a figure to summarise the posterior distribution of 
Log(Ne) in a certain period (TMAX) given by the user. 

The given mutation rate is used to choose natural scales (census size and generations), 
or reduced scales (theta and generations times mutation rate) if the mutation rate is set to 0. \cr


The graph is plotted and saved in a text file job.2D of 101 lines: \cr

Line 1, column 1: rate of coverage of the plotted distributions.

Line 1, columns 2 to 101: abscissa (times from 0 to TMAX when the Log(Ne) distribution is estimated).\cr

Lines 2 to 101, column 1: ordinates (values of Log(Ne) in the selected range).

Lines 2 to 101, columns 2 to 101: densities of Log(Ne) (integer values).

}
\value{
\item{R plot}{In R console: Graphical summary of the posterior distribution of estimates of Log(Ne) in the past.}
\item{job.2D}{Text file, 101 lines and 101 columns.}
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
NTdist (NameBATCH = system.file("extdata/job.Batch", package = "VarEff") 
, MUTAT = 0.01, TMAX = 200)
}
\keyword{ models }
