\name{InputTest}
\alias{InputTest}
\docType{data}
\title{
Infile: Population simulated
}

\description{
This dataset corresponds to a population simulated with a size of 1000 individuals that underwent a bottleneck 200 generations ago to become a population with effective size 100.

The genetic data of this population are samples at 20 markers microsatellites.
}

\format{
The format is:

 num [1:240] 9 4 12 2 0 14 21 15 3 9 ...
 
  }

\details{
The infile is called InputTest and it provides the microsatellite data at 20 markers of one population.
}

\section{Infile format}{
Each markers is describe by 2 lines and the alleles not present are mentionned by zero:

9 (number of alleles at the first locus)

4  12   2   0  14  21  15   3   9 (counts of number of alleles with the same length)

12 (number of alleles at second locus)
 
.

.

.

14 (number of alleles at last locus)

1   4   0   0   0  38  18   7   9   1   0   0   1   1 \cr



It means that if you have a locus with 2 types of alleles of lengths 10 and 12 (number of motifs) you have to mention the unobserved allele the 11 repeat motifs!\cr

So if the alleles 10 and 12 have frequencies 24 and 6 respectively, you have to describe the locus by:

3

24 0 6 \cr

This file is close to MSVAR file and you can convert a MSVAR file to VarEff file at the website \url{https://qgp.jouy.inra.fr} 
               
}

\source{
This data is a simulated set obtained by the model DemoDivMs  from Nikolic and Chevalet.
\url{https://qgp.jouy.inra.fr}
}

\references{
Nikolic N., Butler J., Bagliniere JL., Laughton R., McMyn I.A.G, Chevalet C. 2009. An examination of genetic diversity and effective population size in Atlantic salmon populations, and applications for conservation and management. Genetics Research. 91: 1-18.
}

\keyword{datasets}
