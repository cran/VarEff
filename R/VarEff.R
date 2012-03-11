# MCMC estimation of past effective size from microsatellite data
# August / december 2011
#
# integration of functions:
# Theta() preparation of MCMC calculations
# VarEff(): processing estimation for a single set of markers
# VarEffM(): processing successive calculations for several sets of markers
#            from the same population, or analysing repeated simulated data
# visualisation functions:
#
#
#
#
# Three mutation models are proposed: 
#  Single Step Model (S)
#  Two Phase Model   (T)
#  Geometric Model   (G)
#  in the last two cases, an additional parameter C must be given:
#       the rate of mutations with steps of 2 microsatellite motifs (T)
#    or the ratio of the geometric series fixing the distribution of step lengths (G).
#   (this parameter is assumed fixed and it is NOT estimated in the present version)
#
# A mutation rate u (MUTAT parameter) is used to translate standardized units 
# (theta= 4Nu, and reduced time T = Generation * u) into natural values of
# effective size (census of the diploid population), and of time (number of generations).
# However, since the theory makes only use of reduced units, the actual mutation rate
# (parameter MUTAT) is not estimated from data, it must be provided by the user. 
#
# The following functions include:
#
# Theta(): Preliminary analysis of data, building a parameter file to run VarEff
#
# VarEff(): the program to run to generate the posterior distribution of effective sizes
#            in the past
#
#           Results: 
#           if the function is called as  Result = VarEff(),
#           Result$Theta  is the set of 3 global theta estimates, and means and standard deviation
#                         of Quadratic deviations and of prior probabilities of simulated states
#           Result$Batch  is the matrix of the simulated states
#
#           Both lists are also written to files  .Theta and .Batch
#
# and three functions allowing the results of VarEff to be analysed:
#
#    Note that if the parameter MUTAT is not provided in these functions,
#    results are shown using reduced scales (theta's, and G*u units).
#    Any intermediate result is saved under this reduced format.
# 
# Functions NatSizeDist(...) and LogSizeDist(BATCHfile,TMAX,NBT,MUTAT)
# Extracts Ne (or log_10(Ne) ) posterior distributions at a number of past times
# Input:    the name of the parafile.Batch file produced by VarEff
#           the interval of observation of past times TMAX 
#           the number of observed instants NBT
#           the MUTAT mutation rate (if set to 0, TMAX must be given in reduced time,
#                                    and results are given in reduced scales)
# Results:  global stats (mean of posterior distributions of Ne at past times)
#              (plot and file parafile.Nstat)
#           detailed posterior distribution of Ne at the analysed times
#              (file parafile.Ndist)
#
# Function NTdist(BATCHfile,TMAX,MUTAT)
# Graphical summary of the 2D joint distribution of N and of T
# Input:    the name of the parafile.Batch file produced by VarEff
#           the interval of observation of past times TMAX 
#           the MUTAT mutation rate (if set to 0, TMAX must be given in reduced time,
#                                    and results are given in reduced scales)
# Results:  the 2D joint distribution of T and log_10( N(T) )
#                  plot using scales as defined by parameter MUTAT 
#                        (reduced scales if MUTAT set to 0)
#              and file parafile.2D using the reduced scales
#
# Function plotNdistrib(NDISTfile)
# Plots posterior distributions of N(T) at given times
# Input:    the name of a parafile.Ndist file produced by the function SizeDist()
#           (the choice of scales is made in function SizeDist)
# Results:  plots of the posterior distributions of population size at different
#           past times. Choosing these times and  the range of population sizes
#           to plot is under user's control.
#
###############################################################################################
#  Main function
#
VarEff <- function(infile=NULL,
                   parafile=NULL,
                         NBLOC=NULL,
                         JMAX=NULL,
                         MODEL=NULL,
                         MUTAT=NULL,
                         NBAR=NULL,
                         VARP1=NULL,
                         RHOCORN=NULL,
                         GBAR=NULL,
                         VARP2=NULL,
                         DMAXPLUS=NULL,
                         Diagonale=NULL,
                         NumberBatch=NULL,
                         LengthBatch=NULL,
                         SpaceBatch=NULL)
{
# Load mcmc
library(mcmc)
# Load R Graphics dev2
library(grDevices)
# 
print(" Estimation of effective sizes Ne(T) in the past",quote="FALSE")
print(" from microsatellites data",quote="FALSE")
print(" Note that this R code makes use of the mcmc library ",quote="FALSE")
print(" ",quote="FALSE")
#
#############################################
# Parameters for the model, given by the user
# 
print(" *** Data ***",quote="FALSE")
if(!is.null(infile)){
print (paste(" Name of the data file= ",infile,sep=""),quote="FALSE")}else{
infile = readline(" Name of the data file= ")}
#
if(!is.null(NBLOC)){
print (paste(" Number of microsatellite markers= ",NBLOC,sep=""),quote="FALSE")}else{ 
NBLOC = as.numeric(readline(" Number of microsatellite markers= "))}
#
# If not given, give a name for the job
if(!is.null(parafile)){
print(paste(" Name of the job = ",parafile,sep=""),quote="FALSE") }else
{parafile = readline(" Give a name for this job = ")}
#
print(" *** Demographic model ***",quote="FALSE")
print(" Choose the number of transitions between periods with constant Ne ",quote="FALSE")
if(!is.null(JMAX)){
print (paste(" Value of JMAX = ",JMAX,sep=""),quote="FALSE")}else{
JMAX <- as.numeric(readline(" Value of JMAX = "))}
#
print(" *** Mutation model ***",quote="FALSE")
print(" Fix the mutation model (S, T or G) ",quote="FALSE")
print("   (S = Single Step Model, T = Two Phase Model, G = Geometric Model)",quote="FALSE")
print("   and the prior value of its C parameter,  0 < C < 1",quote="FALSE")
print("   e.g. answer ' T   0.15 ' for the Two Phase Model with a proportion 15% of two motifs steps",quote="FALSE")
if(!is.null(MODEL)){
print (paste(" MODEL = ",MODEL),quote="FALSE")}else{
MODEL = readline(" Mutation model and C coefficient? ")
}
CoGeo = .MutMod(MODEL)
#
print(" Give the mutation rate MUTAT",quote="FALSE")
if(!is.null(MUTAT)) {
print (paste(" MUTAT = ",MUTAT),quote="FALSE")}else{
MUTAT <- as.numeric(readline(" MUTAT = "))}
#
print(" *** Priors about population sizes ***",quote="FALSE")
print(" ",quote="FALSE")
#
print(" Prior expected value of Effective Population Size Ne:",quote="FALSE")
if(!is.null(NBAR)){
print (paste(" Prior Ne= ",NBAR,sep=""),quote="FALSE")}else{
NBAR = as.numeric(readline(" Prior Ne = "))}
#
if(!is.null(VARP1)){
print (paste(" Variance of log(Ne)= ",VARP1,sep=""),quote="FALSE")}else{
VARP1 = as.numeric(readline(" Variance of log(Ne)= "))}
#
print(" Prior correlation between successive population sizes:",quote="FALSE")
#
if(!is.null(RHOCORN)){
print (paste(" correlation= ",RHOCORN,sep=""),quote="FALSE")}else{
RHOCORN = as.numeric(readline(" correlation= "))}
#
print(" *** Priors about time and time intervals ***",quote="FALSE")
print(" ",quote="FALSE")
#
print(" Prior time since the Origin (number of generations))",quote="FALSE")
if(!is.null(GBAR)){
print (paste(" Origin time (generations)= ",GBAR,sep=""),quote="FALSE")}else{
GBAR =as.numeric(readline(" Origin time (generations) = "))}
# 
print(" Prior about time intervals:",quote="FALSE")
#
if(!is.null(VARP2)){
print (paste(" Variance of log(time interval)= ",VARP2,sep=""),quote="FALSE")}else{
VARP2 = as.numeric(readline(" Variance of log(time interval)= "))}
#
# MCMC
#
print(" *** Smoothing the Covariance matrix ***",quote="FALSE")
print(" ",quote="FALSE")
print("    Calculation of an approximate likelihood assumes the theoretical",quote="FALSE") 
print("    variance-covariance matrix V of observations is known.",quote="FALSE") 
print("    It is in fact replaced by a combination of the sample estimate Vs",quote="FALSE")
print("    and of the diagonal matrix Vd made of theoretical variances:",quote="FALSE")
print("    V = (1-D) * Vs + D * Vd",quote="FALSE")
print("     D  must be > 0 and less than 1",quote="FALSE")
print("     D = 1  means that only theoretical variances are considered",quote="FALSE") 
print("            and correlations are ignored ",quote="FALSE")
print("     D = 0  means that only observed data are considered, a choice ",quote="FALSE")
print("            that may raise numerical instability",quote="FALSE")
print("     For more detail read the help and article",quote="FALSE")
print("     A medium choice, D = 0.5 , is generally robust",quote="FALSE")
#
if(!is.null(Diagonale)){
print (paste(" D = ",Diagonale,sep=""),quote="FALSE")}else{
Diagonale= as.numeric(readline(" D = "))}
#
Diagonale = min(Diagonale,1)
if( Diagonale <= 0) 
{
print(" Diagonale must be > 0, it is reset to 0.1",quote="FALSE")
Diagonale = 0.1
}
#
print(" *** Length of the MCMC chain (look at the metrop specifications) ***",quote="FALSE")
print(" Values recommended: nbatch = 10000, blen=10, nspac=10 ",quote="FALSE")
print(" Enter the three values you wish to use",quote="FALSE")
if(!is.null(NumberBatch)){
print (paste(" nbatch = ",NumberBatch,sep=""),quote="FALSE")}else{
NumberBatch = as.numeric(readline(" nbatch = "))}
if(!is.null(LengthBatch)){
print (paste(" blen = ",LengthBatch,sep=""),quote="FALSE")}else{
LengthBatch = as.numeric(readline(" blen = "))}
if(!is.null(SpaceBatch)){
print (paste(" nspac = ",SpaceBatch,sep=""),quote="FALSE")}else{
SpaceBatch = as.numeric(readline(" nspac = "))}
#
#############################################
# Fixing constants of the prior distributions
#
# checking JMAX, a parameter that must be > 0
# if JMAX is entered as 0 (the model with a single constant population
# size, it is changed to  JMAX = 1, with RHOCORN = 0.999 so as to
# simulate a single population size
if( JMAX == 0 )
{
JMAX = 1
RHOCORN = 0.999
}
   JMAX1 <- JMAX+1
   JMAX2 <- JMAX+2
   COEF2 <- - 1/(2*VARP2)
   COEF1 <- - 1/(2*VARP1*(1-RHOCORN^2))
#
# add priors for extensions of the mutation model
   NBPARAM <- 2*JMAX+2
#
# in the present version of the program, the additional parameter C is 
# assumed fixed throughout the calculations
      Prod3 = 0
      COEF3 = 0 
      VARP3 = 0
#      if( !(CoGeo == 0 ) )
#   {
#      Prod3 = 0.5
#      VARP3 = 0.69
#      COEF3 = - 1/( 2* VARP3)
#    }
#
#############################################
# Initialisation of the set of parameters to be searched
#
   parametre = rep(0,NBPARAM)
   THETATAUZERO <- vector(mode="numeric",length=NBPARAM)
   the <- 4*MUTAT*NBAR 
   gba <- GBAR/(2*NBAR*JMAX)
   THETATAUZERO[1:JMAX1] <- the
   THETATAUZERO[JMAX2:(NBPARAM-1)] <- gba
#
#  Priors mean for the model of mutation 
#
   THETATAUZERO[NBPARAM] = CoGeo
#
priors = list(TZ=THETATAUZERO,C1=COEF1,C2=COEF2,C3=COEF3,RO=RHOCORN)
#
##############################################
# Definition of the proposal distribution (matrix scalecal)
#
   SCALE11 <- t(chol(.vdmdiretendue(RHOCORN, JMAX)))
  scalecal <- SCALE11*.prodEE(JMAX,0.5,0.5,Prod3)
  paraminit <- parametre
#
#############################################
# Reading the data from infile
#
datex = .DataExtract(infile,NBLOC)
PkMoy  = datex$P
VarCov = datex$V
NbMaxAll   = datex$A
NELOC      = datex$L
#
print(c(" Effective Number of markers = ",NELOC)  ,quote="FALSE" )
#
# Showing the distribution of P_k from data to choose DMAXPLUS
#
print(" Plot of the (mean) distribution of distances between alleles",quote="FALSE")
###
abscissa = " Distance k between alleles"
ordinates= " Frequencies f_k of pairs of alleles at distance k"
title = paste(" Data file: ", infile)
par(adj=0.5)
plot(c(-0.5,NbMaxAll   +0.5),c(0,max(PkMoy)),xlab=abscissa,ylab=ordinates,main=title,type='n')
lines(-0.5+c(0:(NbMaxAll   -1)),PkMoy,type='h')
esfplus = c(PkMoy,0)
lines(-0.5+c(0:NbMaxAll   ),esfplus,type='s')
#
lines(c(-0.5,NbMaxAll   -0.5),c(0.005,0.005),col="red")
#
# save the plot of f_k frequencies to file
#dev.copy2pdf(file=paste(infile,"_Pk.pdf",sep=""))
#
# Choose DMAXPLUS
#
print(" Choose the range 0:DMAX of allele distances to be analysed",quote="FALSE")
print(" e.g. so that f_k < 0.005 for k > DMAX (under the red line)",quote="FALSE")
if(!is.null(DMAXPLUS)){
print (paste(" DMAX = ",(DMAXPLUS-1),sep=""),quote="FALSE")}else{
DMAXPLUS = 1 + as.numeric(readline(" Enter DMAX = "))}
#
# check that DMAXPLUS is smaller than the largest allele distance
if(DMAXPLUS > NbMaxAll)
{
print(c(" DMAX is reduced to ",(NbMaxAll-1)),quote="FALSE")
DMAXPLUS = NbMaxAll
}
#############################################
#  Save this set of parameters to file 
#
.parDiagwrite(parafile,infile,NBLOC,JMAX,MODEL,MUTAT,NBAR,VARP1,RHOCORN,GBAR,VARP2,DMAXPLUS,Diagonale,NumberBatch,LengthBatch,SpaceBatch)
#
#############################################
# Mathematical resolutions
#############################################
# Markov Chains with a priori and the parameters 
# NbMaxAll, PkMoy, VarCov, DMAXPLUS, 
#
# Global estimates of a single theta value
#
# theta_0, theta_1 et theta_2 are put in thetasim[1:3]
#
thetasim <- vector(mode="numeric",length=7)
#
# Estimation of theta_0, theta_1 et theta_2
#
# Theta_1
theta1 = 2 * sum(c(1:(NbMaxAll-1))*PkMoy[2:NbMaxAll])
theta1 = theta1*(theta1+sqrt(theta1**2+1))
thetasim[2]=theta1
# Theta_2
theta2 = 2*sum((c(1:(NbMaxAll-1))**2)*PkMoy[2:NbMaxAll])
thetasim[3]= theta2
# Theta0
theta0 = PkMoy[1]
theta0 = theta0**2
theta0 = ( 1/theta0 )-1
theta0 = theta0/2
thetasim[1]=theta0
#
# Outputs
#
# Saving the results on files (start)
# Global theta estimates
print(" Mean values of estimates Theta_0, Theta_1, Theta_2",quote="FALSE")
print(thetasim[1:3],quote="FALSE")
print(" Imbalance indices ln(Theta_1/Theta_0) and ln(Theta_2/Theta_0)",quote="FALSE")
ib1 = log(theta1/theta0) 
ib2 = log(theta2/theta0) 
print(c( ib1 , ib2 ),quote="FALSE")
MinNe = as.integer(min(thetasim[1:3])/(4*MUTAT))
MaxNe = as.integer(max(thetasim[1:3])/(4*MUTAT))
print(paste("order of magnitude of Ne from ",MinNe," to ",MaxNe) ,quote="FALSE")
#
write(thetasim[1:3],file=paste(parafile,"Theta",sep="."))
write(c(ib1,ib2),file=paste(parafile,"Theta",sep="."),append="TRUE")
#
write(c(MinNe,MaxNe),file=paste(parafile,"Theta",sep="."),append="TRUE")
#
# Likelihood model
#
# "Diagonale" is entered as a parameter
#
FREQEXPE <- PkMoy[1:DMAXPLUS]
VCDMX = VarCov[1:DMAXPLUS,1:DMAXPLUS]
# The approximate theoretical Diagonal 
diagt = .diagtheo(theta1,DMAXPLUS)
VCDMX = Diagonale*diag(diagt)+ (1-Diagonale)*VCDMX 
VARCOVINV = solve(VCDMX)
#
# Get the constant term of quadratic deviations
XLKZ = 0
XLKZ = t(FREQEXPE) %*% ( VARCOVINV %*% FREQEXPE )
#
# Get the constant term of the prior log-probability
PZZERO = 0
# If(VARP3 > 0) PZZERO = log(VARP3)
PZZERO = PZZERO + NBPARAM * log(2*pi) + JMAX*log(VARP2)
PZZERO = PZZERO + (JMAX+1)*log(VARP1) + JMAX*log(1-RHOCORN^2)
PZZERO = - PZZERO / 2
#
# Keep observed values in the list  "observed"
#
observed = list(DM=DMAXPLUS,NL=NELOC,FR=FREQEXPE,IV=VARCOVINV)
# 
#############################################
# Adjust "scalecal" to get the chosen acceptation rate of MCMC 
REFACC=0.25
#
scalecal = .Refa(.Laposter,paraminit,JMAX,observed,priors,REFACC,scalecal)
# 
#############################################
# Burn-in
out <- metrop(.Laposter,paraminit,nbatch=10000,blen=1,nspac=1,scale=scalecal,,,JMAX,observed,priors)
#
#############################################
# Chain 
out <-  metrop(out,nbatch=NumberBatch,blen=LengthBatch,nspac=SpaceBatch,,,,,,JMAX,observed,priors)
#
#############################################
# Recover Posterior probabilities and likelihoods of the simulated states
#
distlp = .distriL(out$batch,JMAX,observed,priors)
distLike = distlp$L
distPost = distlp$P
#
# Quadratic deviation from expectation (f_obs-f_exp)' Varcovinv (f_obs-f_exp)
Quaddev = - 2 *  distLike + NELOC * XLKZ
#
# Record Quadev and distPost stats in thetasim
thetasim[4] = mean(Quaddev)
thetasim[5] = sd(Quaddev)
thetasim[6] = mean(distPost-distLike) + PZZERO
thetasim[7] = sd(distPost-distLike)
#
# Recover parameters to their natural values (reduced scales theta and G*u times))
tt = .ExtracNatural(out$batch,JMAX,priors$TZ,MUTAT)
# 
#############################################
# Outputs
#
# Saving the results on files
#
print(" Mean and standard deviation of Quadratic deviations from data",quote="FALSE")
print(thetasim[4:5],quote="FALSE")
print(" Mean and standard deviation of log prior probabilities",quote="FALSE")
print(thetasim[6:7],quote="FALSE")
write(thetasim[4:7],file=paste(parafile,"Theta",sep="."),append="TRUE")
#
# The detailed series of simulations
# 
# File .Batch    the full set of simulated values 
#           column 1 = Number of the simulation
#           column 2 = Quaddev
#           column 3 = Prior proba of parameters
#           columns 3+1 to 3+JMAX+1: theta_j
#           columns 4+JMAX+1 to 4+JMAX+JMAX: times t_j in generation*mutation units
#           column  5 + 2*JMAX: the current C value (extended mutation models)
#                      for Geometric model, C is encoded with positive values
#                      for the Two Phase Model, C is encoded as its opposite
#
fullout = matrix(c(1:NumberBatch,Quaddev,distPost-distLike+PZZERO,tt),nrow=NumberBatch)
write(t(fullout),ncolumns=(NBPARAM+3),file=paste(parafile,"Batch",sep="."))
return()
}
##############################################################################
############################################################################## 
# Visualization of the results of VarEff
# 
#    Note that if the parameter MUTAT is set to 0 in these functions,
#    results are shown using reduced scales (theta's, and G*u units).
#    Any intermediate result is saved under this reduced format.
#
##############################################################################
# Effective size distribution functions:
#    natural values (Ne or theta scales):              Function  NatSizeDist(BATCHfile,TMAX,NBT,MUTAT)
# or decimal logarithms (log_10(Ne) or log_10(theta)): Function  LogSizeDist(BATCHfile,TMAX,NBT,MUTAT)
#
# Extracts Ne distribution at a number of past times
# Input:    the name of the parafile.Batch file produced by VarEff
#           the interval of observation of past times TMAX 
#           the number of observed time intervals NBT (hence NBT+1 times from 0 to TMAX)
#           the MUTAT mutation rate (if set to 0, TMAX must be given in reduced time,
#                                    and results are given in reduced scales)
#
# Results:  Global stats (mean of posterior distributions of (log)Ne at past times),
#           given as a table with (NBT+1) lines and 7 columns containing each:
#           time of observation, arithmetic mean, harmonic mean, mode, median, 5% and 95% quantiles.
#           The table is written to file parafile.Nstat, and the results plotted.
#              
#           In addition the densities of the posterior distribution of (log)Ne at the analysed times
#           are written to file parafile.Ndist (to be visualized with the plotNdistrib() function, below).
#
#      technical comment on the code:
#      The values extrait$tK, -$moy, -$eca, -$Q are in units theta = 4 * N * MUTAT 
#      The values  extrait$tech  are in units Generation * MUTAT
#
#   parametre ",NatSizefile=NULL"  SUPPRIME
#
NatSizeDist <- function(NameBATCH=NULL,MUTAT=NULL,TMAX=NULL,NBT=NULL)
{
# get parameters
print(" VarEff - View past effective sizes Ne(T)",quote="FALSE")
print("          and save posterior distributions",quote="FALSE")
if(!is.null(NameBATCH)){
print (paste(" Name of the batch table= ",NameBATCH),quote="FALSE")}else{
NameBATCH = readline(" Name of the batch file= ")}
#
if(is.null(MUTAT))
{
print(" Mutation rate: Enter 0 to use reduced scales Theta's and G*u",quote="FALSE")
MUTAT = as.numeric(readline(" mutation rate = "))
}
#
if(!is.null(TMAX)){
print (paste(" Length of the period = ",TMAX),quote="FALSE")}else{
print(" Enter the length of the period to be analysed",quote="FALSE")
if(MUTAT > 0) print(" give time in number of generations ",quote="FALSE")
if(MUTAT == 0)print(" give time in the reduced time scale G*u ",quote="FALSE")
TMAX = as.numeric(readline(" Length of the period = "))}
#
if(!is.null(NBT)){
print (paste(" Number of time intervals= ",NBT),quote="FALSE")}else{
NBT <- as.numeric(readline(" Number of time intervals= "))}
#
#
if(!(MUTAT == 0))
{
print (paste(" Assumed mutation rate= ",MUTAT),quote="FALSE")
TMAX = TMAX * MUTAT 
}
else
{
print(" Results are expressed in reduced scales Theta's and G*u",quote="FALSE") 
}
#
#  PARTIE SUPPRIMEE  #  PARTIE SUPPRIMEE  #  PARTIE SUPPRIMEE  #  PARTIE SUPPRIMEE
#
# If not given, give a name for this set of parameters 
# if(is.null(NatSizefile)){
# NatSizefile = readline(" Give a name for this set of parameters = ")
# }
# Write the parameters     
#
# info=matrix(c(TMAX,NBT,MUTAT)) # A CE MOMENT DU PROGRAMME, TMAX EST EN UNITES REDUITES
#
# write(info,file=paste(NatSizefile,"NatSizeInfo",sep="."))
#
#  PARTIE SUPPRIMEE  #  PARTIE SUPPRIMEE  #  PARTIE SUPPRIMEE  #  PARTIE SUPPRIMEE
#
.LSD(NameBATCH,TMAX,NBT,MUTAT,FALSE)
return()
}
##############################################################################
# Function LogSizeDist(BATCHfile,TMAX,NBT,MUTAT)
# 
# same as SizeDist() for log-sizes
#
LogSizeDist <- function(NameBATCH=NULL,MUTAT=NULL,TMAX=NULL,NBT=NULL)
{
# get parameters
print(" VarEff - View past effective sizes Ne(T)",quote="FALSE")
print("          and save posterior distributions",quote="FALSE")
if(!is.null(NameBATCH)){
print (paste(" Name of the batch table= ",NameBATCH),quote="FALSE")}else{
NameBATCH = readline(" Name of the batch file= ")}
#
if(is.null(MUTAT))
{
print(" Mutation rate: Enter 0 to use reduced scales Theta's and G*u",quote="FALSE")
MUTAT = as.numeric(readline(" mutation rate = "))
}
if(!is.null(TMAX)){
print (paste(" Length of the period = ",TMAX),quote="FALSE")}else{
print(" Enter the length of the period to be analysed",quote="FALSE")
if(MUTAT > 0) print(" give time in number of generations ",quote="FALSE")
if(MUTAT == 0)print(" give time in the reduced time scale G*u ",quote="FALSE")
TMAX = as.numeric(readline(" Length of the period = "))}
if(!is.null(NBT)){
print (paste(" Number of time intervals= ",NBT),quote="FALSE")}else{
NBT <- as.numeric(readline(" Number of time intervals= "))}
if(!(MUTAT == 0))
{
print (paste(" Assumed mutation rate= ",MUTAT),quote="FALSE")
TMAX = TMAX * MUTAT 
}
else
{
print(" Results are expressed in reduced scales Theta's and G*u",quote="FALSE") 
}
.LSD(NameBATCH,TMAX,NBT,MUTAT,TRUE)
return()
}
#
##############################################################################
# Function NTdist(BATCHfile,TMAX,MUTAT)
# Graphical summary of the 2D joint distribution of N and of T
# Input:    the name of the parafile.Batch file produced by VarEff
#           the interval of observation of past times TMAX 
#           the MUTAT mutation rate (if set to 0, TMAX must be given in reduced time,
#                                    and results are given in reduced scales)
# Results:  the 2D joint distribution of T and log_10( N(T) )
#                  plot using scales as defined by parameter MUTAT 
#                        (reduced scales if MUTAT set to 0)
#              and file parafile.2D using the reduced scales
#
NTdist <- function(NameBATCH=NULL,MUTAT=NULL,TMAX=NULL)
{
# get parametersBATCH
print(" VarEff - View joint distribution of T and Ne(T)",quote="FALSE")
print("          and save this joint posterior distributions",quote="FALSE")
if(!is.null(NameBATCH)){
print (paste(" Name of the batch file= ",NameBATCH),quote="FALSE")}else{
NameBATCH = readline(" Name of the batch file= ")}
#
if(is.null(MUTAT))
{
print(" Mutation rate: Enter 0 to use reduced scales Theta's and G*u",quote="FALSE")
MUTAT = as.numeric(readline(" mutation rate = "))
}
if(!is.null(TMAX)){
print (paste(" Length of the period = ",TMAX),quote="FALSE")
}else{
print(" Enter the length of the period to be analysed",quote="FALSE")
if(MUTAT > 0) print(" give time in number of generations ",quote="FALSE")
if(MUTAT == 0)print(" give time in the reduced time scale G*u ",quote="FALSE")
TMAX = as.numeric(readline(" Length of the period = "))
}
if(MUTAT == 0) MUTAT=NULL
if(!is.null(MUTAT))
{
print (paste(" Assumed mutation rate= ",MUTAT),quote="FALSE")
TMAX = TMAX * MUTAT 
}
else
{
print(" Results are expressed in reduced scales Theta's and G*u",quote="FALSE") 
}
parafile = substring(NameBATCH,1,nchar(NameBATCH)-6)
BATCH = read.table(NameBATCH)
Lenfull = length(BATCH[1,])
NumberBatch = length(BATCH[,1])
JMAX = (Lenfull-3)/2 - 1 
tt = BATCH[,4:Lenfull]
#
Ntime = 100
Nsize = 100
#
map = .LogMap(tt,TMAX,Ntime, Nsize,JMAX)
#
MinLN = map$sm 
MaxLN = map$SM 
Cover = 100*sum(map$M)/(NumberBatch*Ntime)
# save to file, adding 1st line for times, 1st column for log_10(theta)'s
abscissa = TMAX*c(1:Ntime)/Ntime
ordinates = (MinLN+c(1:Nsize)*(MaxLN-MinLN)/Nsize)/log(10)
mapww = matrix(0,ncol = (Nsize+1),nrow=(Ntime+1) )
mapww[1,1] = Cover
mapww[2:(Nsize+1),1] = abscissa
mapww[1,2:(Ntime+1)] = ordinates
mapww[2:(Nsize+1),2:(Ntime+1)] = map$M
 write(t(mapww),ncolumns=(Nsize+1),file=paste(parafile,"2D",sep=".")) 
# 
MinLN = ( map$sm )/log(10)
MaxLN = ( map$SM )/log(10)
#
if(is.null(MUTAT))
{
print( " Plotting the 2D ( T,log_10(Theta) ) distribution", quote="FALSE")
print( c(" Time range: from 0 to ",TMAX," (G*u units)"),quote="FALSE")
print( " Minimum and Maximum values of log_10(theta) in the plot:",quote="FALSE")
print( c("       Minimum = ",MinLN),quote="FALSE")
print( c("       Maximum = ",MaxLN),quote="FALSE")
}
else
{
print( " Plotting the 2D ( g ,log_10( Ne(g) ) ) distribution", quote="FALSE")
print( c(" Time range: from 0 to ",TMAX/MUTAT," (generations)"),quote="FALSE")
delta = log(4*MUTAT)/log(10)
print( " Minimum and Maximum values of log_10( Ne ) in the plot:",quote="FALSE")
print( c("       Minimum = ",MinLN-delta),quote="FALSE")
print( c("       Maximum = ",MaxLN-delta),quote="FALSE")
}
print( c("       Coverage  = ",Cover, " %"),quote="FALSE")
#
# plot
titre = paste("N-T joint distribution (coverage = ",Cover," % )")
if(is.null(MUTAT))
{
legord=" log_10( Theta(T) )"
legabs= " Time T in the past (unit = G*u)"
TMT = TMAX / Ntime
MinOrdo = MinLN
} 
else
{
legord=" log_10( N(T) )"
legabs= " Time T in the past (generations)"
TMT = TMAX/(Ntime*MUTAT)
MinOrdo = MinLN - delta
}
par(adj=0.5)
image((TMT*c(1:Ntime)),MinOrdo+ c(1:Nsize)*(MaxLN-MinLN)/Nsize,(-map$M),xlab=legabs,ylab=legord)
par(adj=0)
title(main=titre,sub=parafile,cex.main=1,cex.sub=1)
#
# save the graph to pdf file
#if(!is.null(parafile)) { dev.copy2pdf(file=paste(parafile,"_2D.pdf",sep="") ) }
#else {dev.copy2pdf(file=paste(parafile,"_2D.pdf",sep="")}
#return() 
# save the graph to pdf file
#dev.copy2pdf(file=paste(parafile,"_2D.pdf",sep="") )
return()
}
#
##############################################################################
# Function plotNdistrib(NDISTfile)
# Plots posterior distributions of N(T) at given times
# Input:    the name of a parafile.Ndist file produced by the function functions 
#           NatSizeDist() and LogSizeDist()
#           (the choice of scales is made in function SizeDist)
# Results:  plots of the posterior distributions of population size at different
#           past times. Choosing these times and  the range of population sizes
#           to plot is under user's control.
#
plotNdistrib <- function(infile=NULL, nbcases=NULL)
{
print(" *** Plot posterior distributions of N(T) ***",quote="FALSE")
if(!is.null(infile)){
print (paste(" Name of the .Ndist file or .Ldist file= ",infile),quote="FALSE")}else{
infile = readline(" Name of the 'parafile.Ndist' file or 'parafile.Ldist' = ") }
#
parafile = substring(infile,1,nchar(infile)-6)
#
dat=read.table(infile)
# check number of times and show intermediate times
Nbtimes = length(dat[,1])
# check the format of .Ndist file
# 1025 corresponds to log option
XAvail= FALSE
if( length(dat[1,]) == 1025 ) XAvail=TRUE
#
print("Instant       past time ",quote="FALSE")
for (time in 1:Nbtimes) print(paste(time,"                 ",dat[time,1]),quote="FALSE")
# choose the number of instants you wish to plot the distribution
#
#nbcases = 1
#
#while (nbcases > 0) {
#
if(!is.null(nbcases)){
print("Choose the number of times (suggestion: <=5) you wish to plot the distribution",quote="FALSE")}else{
nbcases = as.integer(readline("  number of times (stop on 0): "))
}
if (nbcases > 0) {
if(nbcases > 11) {
return( print( paste(nbcases," is too large a number. Please retry !"),quote="FALSE") )
}
if(nbcases == 0) {
# save the last graph to pdf file
if(XAvail) { dev.copy2pdf(file=paste(parafile,"_Ldist.pdf",sep="") ) }
else { dev.copy2pdf(file=paste(parafile,"_Ndist.pdf",sep="") )  }
return()
                 }
cases = c(1:nbcases)
for (ti in 1:nbcases)
{
cases[ti] = readline(paste(" case # ",ti,", instant number = : "))
}
# define colours
colours = c("black","blue","red","green","grey","purple","orange","pink","brown","yellow","turquoise")
time = cases[1]
if(XAvail) 
{
x = as.numeric(dat[time,2:513])
y = as.numeric(dat[time,514:1025])
}
else
{  
x=as.numeric(dat[time,2])
x= x*c(0:511)
y=as.numeric(dat[time,3:514])
}
z=as.integer( dat[time,1] )
#
if( nbcases == 1) {
legabs =  " Ne  values at time  -"
if(XAvail) legabs =  " log_10(Ne)  values at time  -"
par(adj=0.5)
plot(x,y,xlab=paste(legabs,z,"   "), ylab="Density",type="l")
par(adj=0)
title(sub= parafile)
}
else {
# with several distributions, check limits
if(XAvail)
{
xmax = max( dat[cases,513] )
xmin = min( dat[cases,2] )
}
else
{
xmax = 511 * max( dat[cases,2] )
xmin = 0
}
# choice of the range of (log)N values
print(" Choosing the range of Ne values ",quote="FALSE" )
if(XAvail) print("  note that a log_10 scale of sizes is used here",quote="FALSE" )
print(paste("  maximal  Ne  value is ",xmax,". Give an upper value:"),quote="FALSE" )
xmax= readline("  enter maximum  Ne  value= ")
print(paste("  minimal  Ne  value is ",xmin,". Give a lower value:"),quote="FALSE" )
xmin= readline("  enter minimum  Ne  value= ")
if(XAvail)
{
ymax = max(dat[cases,514:1025])
legy = "density of posterior distribution of log_10(N(T)) "
legx = " log_10(Ne)  values "
}
else
{
ymax = max(dat[cases,3:514])
legy = "density of posterior distribution of N(T) "
legx = " Ne  values "
}
#
TITRE = paste(colours[1:nbcases],": at T= ",(dat[cases,1])," in the past")
# TITRE = paste("Time in the past \n", colours[1:nbcases], "\n", (dat[cases,1])
par(adj=0.5)
plot(c(xmin,xmax),c(0,ymax),type='n',xlab=legx,ylab=legy)
par(adj=0)
title(sub= parafile,main=TITRE,cex.main=1,cex.sub=1)
#
for (ti in 1:nbcases) {
time = cases[ti]
if(XAvail) 
{
x = as.numeric(dat[time,2:513])
y = as.numeric(dat[time,514:1025])
}
else
{  
x=as.numeric(dat[time,2])
x= x*c(0:511)
y=as.numeric(dat[time,3:514])
}
lines(x,y,col=colours[ti])
}
     }
#  end of while loop
               }
return()
}
##############################################################################
##############################################################################
# # Functions with explicit transfer of parameters
##############################################################################
# Estimation of the likelihood of thetatau (parameters N_i and differences  t_(j) - t_(j-1)  )
.Logvraisembl <- function(thetatau,JMAX,observed)
{
DMAXPLUS = observed$DM
espfreqdis <- vector(mode="numeric",length=DMAXPLUS)
# Means of theorics fk
  for (dista in 1:DMAXPLUS)
 { 
   intint = integrate(.fonvarjmax2 ,0,pi,dista-1,JMAX,thetatau,stop.on.error = FALSE)
   if(intint$message == "OK"){espfreqdis[dista] <- intint$value/pi}else{return(-Inf)}
  }
# Vector of VarCovInv * espfreqdis
vciprod <- observed$IV %*% espfreqdis
return(observed$NL*sum(vciprod*(2*observed$FR-espfreqdis))/2)
}
############################################################################## 
# Posterior probability of parameters
.Laposter <- function(parametre,JMAX,observed,priors)
{
COEF1 = priors$C1
COEF2 = priors$C2
COEF3 = priors$C3
RHOCORN = priors$RO
thetatau <- priors$TZ * exp(parametre)
# check for c < 1
if( abs(thetatau[2*JMAX+2]) > .99) return(-Inf)
valeur <- .Logvraisembl(thetatau,JMAX,observed)
if(valeur == -Inf) return(-Inf)
# 
 ap1 <- COEF2*sum(parametre[(JMAX+2):(2*JMAX+1)]^2)
 ap2 <- sum(parametre[1:(JMAX+1)]^2)
 ap2 = ap2 + RHOCORN *(RHOCORN * sum(parametre[2:JMAX]^2)
     - 2 * sum(parametre[1:JMAX] * parametre[2:(JMAX+1)]))
#
# add the term for "C"
return(valeur + ap1 + COEF1 * ap2 + COEF3 * parametre[2*JMAX+2]^2)
}
##############################################################################
# Function submitted to the Fourier inversion, to get expected values of the fk
.fonvarjmax2  <- function(x,dista,JMAX,thetatau)
{
        qmuc <- 1 - .geom(x,thetatau[2*JMAX+2])
        jplus <- JMAX + 1
        valeur <- 1/(1 + thetatau[jplus] * qmuc)
        j <- JMAX 
        jplus <- jplus + JMAX 
        {
                while(j > 0) {
                        ajet <- 1 + thetatau[j] * qmuc
                        beta <- exp( - thetatau[jplus] * ajet)
                        valeur <- valeur * beta + (1 - beta)/ajet
                        j <- j - 1
                        jplus <- jplus - 1
                }
        }
        return(valeur * cos(dista * x))
}
##############################################################################
#Goemetric function
#TPM, C < 0
.geom <- function(x,c)
{
cosx = cos(x)
if( c == 0 ) return(cosx)
#
if( c < 0 ) return((1+c)*cosx - c * cos(2*x))
#
a = cosx - c
ac= c*a
return((a-ac)/(1-ac-ac-c*c))
}
##############################################################################
# Matrix function to set the scale matrix
.prodEE <- function(j,V1,V2,V3)
{
j1 <- j +1
j2 <- j +2
Nb <- 2*j1
 param1 = rep(0,Nb)
 param2 <- param1
 param3 <- param1
param1[1:j1]      = 1
param2[j2:(Nb-1)] = 1
param3[Nb]        = 1
return( V1* param1 %*% t(param1)
        + V2 * param2 %*% t(param2) 
        + V3 * param3 %*% t(param3))
}
##############################################################################
# Matrix function to set the scale matrix
.vdmdiretendue <- function(rho,jm)
 {
 djmp <- 2*jm+2
vander = diag(rep(1,djmp))
for (i in 1:jm) {
  for (j in (i+1):(jm+1)) {
  vander[j,i] <- rho^(j-i) 
  vander[i,j] = vander[j,i]
  }
}
 return(vander) 
}
##############################################################################
# Calculates the likelihoods and posterior probabilities of the series (Iterations) 
# of results tt = out$batch
.distriL <- function(tt,JMAX,observed,priors)
{
Nbsim <- length(tt[,1])
distriLike <- vector(mode="numeric",length=Nbsim)
distriPost <- vector(mode="numeric",length=Nbsim)
   for (simu in 1:Nbsim) {
    Lapos2 = .Laposter2(tt[simu,],JMAX,observed,priors)
    distriLike[simu] <- Lapos2$L
    distriPost[simu] <- Lapos2$P
                         }
return(list(L=distriLike,P=distriPost))
}
##############################################################################
# Return the likelihoods and posterior probabilities of parameters
.Laposter2 <- function(parametre,JMAX,observed,priors)
{
COEF1 = priors$C1
COEF2 = priors$C2
COEF3 = priors$C3
RHOCORN = priors$RO
thetatau <- priors$TZ * exp(parametre)
# check for c < 1
if( thetatau[2*JMAX+2] > .99 ) return( -Inf )
valeur <- .Logvraisembl(thetatau,JMAX,observed)
# Fixe the COEF2 <- 1/(2*VARP2)  and  COEF1 <- - 1 /( 2* VARP1 *(1- RhoCorN^2))
ap1 <- COEF2*sum(parametre[(JMAX+2):(2*JMAX+1)]^2)
 ap2 <- sum(parametre[1:(JMAX+1)]^2)
 ap2 = ap2 + RHOCORN*(RHOCORN*sum(parametre[2:JMAX]^2)
     - 2*sum(parametre[1:JMAX]*parametre[2:(JMAX+1)]))
#
# add term for "C"
valpost = valeur + ap1 + COEF1 * ap2 + COEF3 * parametre[2*JMAX+2]^2 
return(list(L=valeur,P=valpost))
}
##############################################################################
# Theoretical approximate value of the diagonal of the variance covariance matrix of the fk
.diagtheo <- function(theta,kmax)
{
Expa = 0.05307618
Bthe = 1.1407126
coefth = exp(-2/sqrt(1+2*theta))
dth = c(1:kmax)
dth[1] = Expa/(theta**Bthe)
if( kmax > 1 ) {
for (k in 2:kmax)
   {
   dth[k] = dth[k-1] * coefth
   }
               }
return(dth)
}
##############################################################################
# From "out$batch" of metrop, extraction of natural parameters
# input : tt = out$batch
#         prio = prior values of theta and tau parameteres 
#                (units 4Nu and g / 2N and  c_0 )
# output
# "tt" is the transformed result from "out$batch" where
#    log parameters (1:JMAX+1) are returned to their natural scale (theta)
#    and parameters (JMAX+2:2*JMAX+2) are transformed from log-intervals
#                                     to tansitions times (scale g_i u)
#
# suppression du N_0 si t_1 < une generation, ie t_1 < MUTAT
#             en imposant theta_0 = theta_1
#
.ExtracNatural <- function(tt,JMAX,prio,MUTAT)
{
 Nbsim <- length(tt[,1])
for (sim in 1:Nbsim) {tt[sim,] <- prio * exp( tt[sim,])}
# Units = generation * mutation
 for (sim in 1:Nbsim) {
   gmu <- 0
            for  (j in 1:JMAX) {   
            gmu <- gmu + 0.5 * tt[sim,(JMAX+1+j)] * tt[sim,j]
            tt[sim,(JMAX+1+j)] <- gmu 
                               
# overall correction N_(i+1)  -->  N_(i)  if  T_(i) < MUTAT
for (jt in (2*JMAX+1):(JMAX+2))
{
if(tt[sim,jt] < MUTAT) tt[sim,(jt-JMAX-1)] = tt[sim,(jt-JMAX)] 
}                              }
                      }
return(tt)
}
##############################################################################
# From "out$batch" of metrop, extraction of the distribution of N(t ) 
# for a set of  analysed times (NbInstants)
# "tt" is the transformed result from "out$batch" where
#    log parameters (1:JMAX+1) are returned to their natural scale (theta)
#    and parameters (JMAX+2:2*JMAX+2) are transformed from log-intervals
#                                     to tansitions times (scale g_i u)
# This transformation is performed by function ExtracNatural(...)
#           
.ExtracTheta <- function(tt,TempsMax,Nbinstants,JMAX)
{
 Nbsim <- length(tt[,1])
#Distribution of theta in function of analysed times (NbInstants)
Longint <- TempsMax / Nbinstants
Kmax <- Nbinstants + 1
ThetaK <- matrix(nrow=Nbsim,ncol=Kmax)
ThetaK[,1] <-  tt[,1]
deltak <- 0
{ for (k in 2:Kmax) 
  {  deltak <- deltak + Longint
      { for (sim in 1:Nbsim) 
         {  
Nbdep <-1
{ for (j in 1:JMAX) 
  { if( deltak >=  tt[sim,(JMAX+1+j)] ) Nbdep <- Nbdep + 1 }
}
             ThetaK[sim,k] <- tt[sim,Nbdep]   }
       }
  }
}
# Mean 
# July 2011 "ecartype" is used to keep harmonic means
moyennes <- vector(mode="numeric",length=Kmax)
 ecartype <- vector(mode="numeric",length=Kmax)
 for  (k in 1:Kmax)
{ 
moyennes[k]  <- mean( ThetaK[,k])
ecartype[k]  <- 1 / mean( 1 / ThetaK[,k]) 
}
# Sample times 
TempEch <- vector(mode="numeric",length=Kmax)
  deltak <- 0
  TempEch[1] <- 0
 { for (k in 2:Kmax){ deltak <- deltak + Longint
 TempEch[k] <- deltak }}
# Quantiles
 QuantP3 <- matrix(nrow=Kmax,ncol=5)
 {for (k in 1:Kmax) {
  QuantP3[k,] <- quantile(ThetaK[,k],probs=c(0.01,0.05,0.5,0.95,0.99),type=4)
 }}
return( list(tK=ThetaK,moy=moyennes,eca=ecartype,tech=TempEch,Q=QuantP3) )
}
##############################################################################
# From trajectories  {t - N} in natural scales, table tt,
# "tt" is the transformed result from "out$batch" where
#    log parameters (1:JMAX+1) are returned to their natural scale (theta)
#    and parameters (JMAX+2:2*JMAX+1) are transformed from log-intervals
#                                     to tansitions times (scale g_i u)
# This transformation is performed by function ExtracNatural(...)
#  
# build a map, a Ntime * Nsize  matrix in which element (i,j) is
# the density of the joint (T,N) distribution
#
# range is  0-TempsMax for times
# the range for sizes is adjusted, censoring highest and lowest sizes of proba < 0.10
#
# data :  tt, TempsMax, Ntime, Nsize, JMAX
# output: map density, in log scale (natural log) for theta's
#         
.LogMap <- function(tt,TempsMax,Ntime, Nsize,JMAX)
{
Nbsim <- length( tt[,1] )
map = matrix(0,nrow=Ntime,ncol=Nsize)
deltaTime = TempsMax/Ntime
# calculate upper size
maxiz = vector(mode="numeric",length=Nbsim)
for (sim in 1:Nbsim) { maxiz[sim] = max( tt[sim,1:(JMAX+1)] ) }
miniz = vector(mode="numeric",length=Nbsim)
for (sim in 1:Nbsim) { miniz[sim] = min( tt[sim,1:(JMAX+1)] ) }
SizeMax = as.numeric( quantile( maxiz, probs=c(0.95), type=4 ) )
SizeMin = as.numeric( quantile( miniz, probs=c(0.05), type=4 ) )
#
deltaSize = log( SizeMax/SizeMin )/Nsize
#
Limit = vector(mode="numeric",length=JMAX)
#
proba = 1
#
for (sim in 1:Nbsim) {
# identify times of transition 
  for (j in 1:JMAX)  {
     Lim      = as.integer( tt[sim,(JMAX+1+j)] / deltaTime ) 
     Limit[j] = min( Lim , Ntime )
                    }
#
# fill map 
# first interval
  if( Limit[1] >0 )
    {
    ordo = as.integer( log( tt[sim,1]/SizeMin ) / deltaSize )
    if( ((ordo+1)>0) & (ordo<Nsize) ) {
    ordo = ordo + 1
    for (ell in 1:Limit[1] ) 
       { map[ell,ordo] = map[ell,ordo] + proba }
                   }
    }
# 2nd to (JMAX)th intervals (provided JMAX > 1)
             if ( JMAX > 1 ) {
  for (j in 2:JMAX) {
       if( Limit[j] > Limit[j-1] )
          {
           ordo =  as.integer( log( tt[sim,j]/SizeMin) / deltaSize )
           if(((ordo+1)>0) & (ordo<Nsize) ) {
           ordo = ordo + 1
           for (ell in (Limit[j-1]+1):Limit[j] ) 
             { map[ell,ordo] = map[ell,ordo] + proba }
                          }
          }
                    }
# last interval
   if( Limit[JMAX] < (Ntime-1) )
     {
      ordo =  as.integer( log( tt[sim,(JMAX+1)]/SizeMin) / deltaSize )
      if(((ordo+1)>0) & (ordo<Nsize) ) {
      ordo = ordo + 1
      for (ell in (Limit[JMAX]+2):Ntime )
         { map[ell,ordo] = map[ell,ordo] + proba }
                       }
    }
                           }
############    end "sim" loop
}
return( list( M=map,SM=log(SizeMax),sm=log(SizeMin) ) )
}
##############################################################################
#
# read data and extract allele frequencies
#           skipping loci with 0 sample size
#
# input: infile, data file
#          NSIM,  number of repetitions (simulations)
#          NBLOC, number of markers
#
# output: NbMaxAll, PkMoyTout, VarCovTout
#          and the effective number of loci: NELOC
#
.DataExtract <- function(infile,NBLOC)
{
# Read the data file
#
TOUT <- scan(infile)
#

#
# First step: derive the maximum number of  alleles "NbMaxAll" 
#             from one or several samples of the population: 
#
NbMaxAll = 1
posit = 0
#
# Analyse each sample (simu)
#

# Look the number of alleles at each makers 
for (locus in 1:NBLOC) {
posit = posit +1
nbal = TOUT[posit]
# Find the allele max
NbMaxAll = max(NbMaxAll,nbal)
# Step to the next maker 
posit = posit + nbal
}
# End
# 
# Second step: Extracting frequencies
#
# Read again
posit = 0
Genot = vector(mode="numeric",length=NbMaxAll)
PkMoy = vector(mode="numeric",length=NbMaxAll)
FreqDis = matrix(nrow=NBLOC,ncol=NbMaxAll)
VarCov = matrix(nrow = NbMaxAll,ncol=NbMaxAll)
# 
# Save the number of alleles at each markers "Nall[locus] = nbal"
#
locut = 0
	for (locus in 1:NBLOC) {
	posit = posit +1
	nbal = TOUT[posit]
	Genot = Genot * 0
   		for (al in 1:nbal) {
   		posit = posit + 1
   		Genot[al] = TOUT[posit]
   		}
# Frequencies
   		nbgenes = sum( Genot[] )
# check sample size
		if( nbgenes > 2 )  {
#
		locut = locut + 1
		if ( nbgenes < 10 ) {
 		print(c(" Marker ",locus," , low sample size = ",nbgenes), quote="FALSE")
					   }
#
 		nbgen2  = sum( Genot[]**2 )
# FreqDis(0) = (nbgen2 - nbgenes) / ( nbgenes*(nbgenes-1) )
		FreqDis[locut,1] = nbgen2 - nbgenes
# Initialization
		FreqDis[locut,2:NbMaxAll] = 0
# Browse the pairs of alleles  al1, al2 
  		if( nbal > 1 )  {
   			for (al1 in 1:(nbal-1) ) { for (al2 in (al1+1):nbal ) {
   			delta = al2-al1+1
   			FreqDis[locut,delta] =FreqDis[locut,delta] + Genot[al1]*Genot[al2]
   			}}
            	}
# Standardization
		FreqDis[locut,] = FreqDis[locut,] / ( nbgenes*(nbgenes-1) )
}
					}
# End
if ( locut < NBLOC ) FreqDis[1:locut,] -> FreqDis
# Extracting PkMoy, VarCov
for (delta in 1:NbMaxAll ) { PkMoy[delta] <- mean( FreqDis[,delta] ) }
 VarCov = var( FreqDis )
return( list(L=locut,A=NbMaxAll,P=PkMoy,V=VarCov) )
}
############################################################################## 
# function to write parameters to a file
.parDiagwrite  <- function(parafile,infile,NBLOC,JMAX,MODEL,MUTAT,NBAR,VARP1,RHOCORN,GBAR,VARP2,DMAXPLUS,Diagonale,NumberBatch,LengthBatch,SpaceBatch)
#
{
parval = vector(length=15)
#
parval=c(parafile,infile,NBLOC,JMAX,MODEL,MUTAT,NBAR,VARP1,RHOCORN,GBAR,VARP2,DMAXPLUS,Diagonale,NumberBatch,LengthBatch,SpaceBatch)
tete = vector(length=16)
tete = c("VarEff(parafile = ","infile = ","NBLOC = ","JMAX = ","MODEL = ",
"MUTAT = ","NBAR = ",
"VARP1 = ","RHOCORN  = ","GBAR = ","VARP2 = ","DMAXPLUS = ",
"Diagonale = ",
"NumberBatch = ","LengthBatch = ","SpaceBatch = ")
#
ligne = vector(length=15)
ligne[c(1:2,5)]=paste(tete[c(1:2,5)],"'",parval[c(1:2,5)],"'",",",sep="")
ligne[c(3:4,6:15)]=paste(tete[c(3:4,6:15)],parval[c(3:4,6:15)],",",sep="")
ligne[16]=paste(tete[16],parval[16],")",sep="")
ecrit= matrix(ligne,nrow=length(ligne))
write(ecrit,ncolumns=1,file=paste(parafile,".R",sep=""))
}
#
############################################################################## 
# read the mutation model
.MutMod <- function(Answer)
{
Coef = 0
rep = .extrTC(Answer)
if( rep$A == "S" ) return(Coef)
if( rep$A == "G" ) Coef = rep$C
if( rep$A == "T" ) Coef = - rep$C
return(Coef)
}
.extrTC <- function(x)
{
#extract first  non blank character
while(substring(x,1,1) == " ") {x = substring(x,2) }
Rep = substring(x,1,1)
if(Rep == 0) {Coef=0} else {Coef=as.numeric(substring(x,2))}
return(list(A=Rep,C=Coef))
}
############################################################################## 
# Adjust "scalecal" to get the chosen acceptation rate of MCMC 
# 
# Input:  the posterioir prob function Lapo
#          paraminit
#          JMAX and the lists observed and priors
#          the requested acceptation rate:  refacc
#          the firts value of the proposal scalecal
# Output:  the new scalecal value
.Refa <- function(Lapo,paraminit,JMAX,observed,priors,refacc, scalecal)
{
# refacc = 0.25
out <- metrop(Lapo,paraminit,nbatch=1000,blen=1,nspac=1,scale=scalecal,,,JMAX,observed,priors)
acc1 = out$accept
acc2 = 0.0
coef = 0.0
if( abs(acc1-refacc) > 0.02 )
{
if (acc1 < refacc){coef=0.2}else{coef=5}
scalecal = scalecal * coef
out <- metrop(Lapo, paraminit,nbatch=1000,blen=1,nspac=1,scale= scalecal,,,JMAX,observed,priors)
acc2 = out$accept
#
# Linear iteration 
while ( abs(acc2-refacc) > 0.02 )
{
# print( paste(acc1,acc2,coef) )
if ( abs(acc1-acc2) > 0.01 )
{ 
   coef = ( (refacc-acc2) / coef + acc1 - refacc )/(acc1-acc2) 
   if( coef < 0 ) coef = 0.5
   if( coef > 5 ) coef = 5.0
}
else { if (acc1 < refacc){coef=0.2}else{coef=5} }
scalecal = scalecal * coef
out <- metrop(Lapo, paraminit,nbatch=1000,blen=1,nspac=1,scale= scalecal,,,JMAX,observed,priors)
acc1 = acc2
acc2 = out$accept
}
}
return(scalecal)
}
############################################################################## 
# function .LSD(NameBATCH=NULL,TMAX=NULL,NBT=NULL,MUTAT=NULL,optionlog)
#
.LSD <- function(NameBATCH,TMAX,NBT,MUTAT,optionlog)
{
#
parafile = substring(NameBATCH,1,nchar(NameBATCH)-6)
#
BATCH = read.table(NameBATCH)
#
Lenfull = length(BATCH[1,])
JMAX    = (Lenfull-3)/2 - 1
tt      =  BATCH[,4:Lenfull]
#
# scaling according to MUTAT
# on entry, TMAX is in reduced scale G*MUTAT
TMAXplot = TMAX
if(!(MUTAT == 0) )
{
TMAXplot = TMAX /MUTAT
tt[,1:(JMAX+1)] = tt[,1:(JMAX+1)] /(4*MUTAT)
tt[(JMAX+2):(2*JMAX+1)] = tt[(JMAX+2):(2*JMAX+1)] / MUTAT
}
#
distriresu = matrix(0,nrow=(NBT+1),ncol=514)
#
#  option LOG : passage en log_10(theta ou N)
if(optionlog) 
{
tt[,1:(JMAX+1)] = log( tt[,1:(JMAX+1)] )/log(10)
distriresu = matrix(0,nrow=(NBT+1),ncol=1025)
}
#
# call ExtracTheta
extrait = .ExtracTheta(tt,TMAXplot,NBT,JMAX)
Nbresu = 0
 Resus  <- rep(0, (6*(NBT+1)) ) 
 for   (tem in 1:(NBT+1))
{
Nbresu = Nbresu+1
Resus[Nbresu] = extrait$moy [tem]
Nbresu = Nbresu+1
if(!optionlog) Resus[Nbresu] = extrait$eca[tem]
Nbresu = Nbresu+1
#
# Mode estimates by density function
liminf = 0
if(optionlog) liminf = extrait$Q[tem,1]
#  densNtem <- density(extrait$tK[,tem],from=liminf,to=extrait$Q[tem,5],adjust=0.5)
densNtem <- density(extrait$tK[,tem],from=liminf,to=extrait$Q[tem,5])
mode = densNtem$x[which.max(densNtem$y)]
#
# save the density
if(optionlog)
{
distriresu[tem,2:513]= densNtem$x[]
distriresu[tem,514:1025]=densNtem$y[]
}
else
{
distriresu[tem,2]= (densNtem$x[512]-liminf) / 511
distriresu[tem,3:514]=densNtem$y[]
}
#
Resus[Nbresu]=mode
#
Nbresu=Nbresu+1
Resus[Nbresu]=extrait$Q[tem,3]
Nbresu=Nbresu+1
Resus[Nbresu]=extrait$Q[tem,2]
Nbresu=Nbresu+1
Resus[Nbresu]=extrait$Q[tem,4]
}
#
#############################################
# Outputs
#
# Collect the results from "Resus" to give some global values 
# and to built some tables in function of "Tempsmax" (-> TMAX)
#
# Mean, mode and median of the effective size in function of the time
#
nbtemps  = length(Resus[])/6
meanmoy= vector(mode="numeric",length=nbtemps)
meanmod= vector(mode="numeric",length=nbtemps)
meanmed= vector(mode="numeric",length=nbtemps)
#
meanhar= vector(mode="numeric",length=nbtemps)
#
 {for (tem in 1:nbtemps) meanmoy[tem] = Resus[1+6*(tem-1)]}
#
 {for (tem in 1:nbtemps) meanmod[tem] = Resus[3+6*(tem-1)]}
#
 {for (tem in 1:nbtemps) meanmed[tem] = Resus[4+6*(tem-1)]}
#
 {for (tem in 1:nbtemps) meanhar[tem] = Resus[2+6*(tem-1)]}
#
# File .Nstat summary statistics of the N(T)
# 
sortieresu = matrix(ncol=7,nrow=nbtemps)
#
# adjust size scale to MUTAT
#
legabs = " Time T in the past (unit T= G*u) "
legord = " Theta(T)"
if(optionlog) legord = " log_10( Theta(T) )"
if(!(MUTAT == 0)) 
{ 
legabs = " Time G in the past (generations)"
legord = " Ne(G)"
if(optionlog) legord = " log_10( Ne(G) )"
}
#
for (tem in 1:nbtemps) 
{
sortieresu[tem,1]=TMAXplot*(tem-1)/(NBT)
for (j in 1:6) {sortieresu[tem,(j+1)]=Resus[(j+6*(tem-1))]}
}
#
if(optionlog)
{ write(t(sortieresu),ncolumns=7,file=paste(parafile,"Lstat",sep=".")) }
else
{ write(t(sortieresu),ncolumns=7,file=paste(parafile,"Nstat",sep=".")) }
#
print(" Plot of mean results ",quote="FALSE")
#
nbtm = nbtemps-1
echelle = TMAXplot/nbtm
 topordo = max(meanmoy,meanmod,meanmed,meanhar)+ 0.1
 botordo = 0
 if(optionlog) botordo = min(meanmoy,meanmod,meanmed)- 0.1
#
par(adj=0.5)
 plot(c(0,TMAXplot),c(botordo,topordo),type='n',xlab=legabs,ylab=legord)
#
par(adj=0)
if(!optionlog) title(sub=parafile,main="Mean (red), Mode (blue), Median (black) and Harmonic mean (orange)",cex.main=1,cex.sub=1)
if(optionlog) title(sub=parafile,main="Mean (red), Mode (blue) and Median (black) of log-sizes ",cex.main=1,cex.sub=1)
#
 lines(c(0:(nbtemps-1))*echelle,meanmod,col='blue')
 lines(c(0:(nbtemps-1))*echelle,meanmed,col='black')
 lines(c(0:(nbtemps-1))*echelle,meanmoy,col='red')
if(!optionlog)  lines(c(0:(nbtemps-1))*echelle,meanhar,col='orange')
#
# File .Ndist
# the densities of the posterior N(T) distributions at the chosen times
#
for (tem in 1:nbtemps) 
{
distriresu[tem,1]=TMAXplot*(tem-1)/(NBT) 
}
if(optionlog)
{
write(t(distriresu),ncolumns=1025,file=paste(parafile,"Ldist",sep="."))
# save the plot
#dev.copy2pdf(file=paste(parafile,"_LLL.pdf",sep="") )
#return()
}
else
{
write(t(distriresu),ncolumns=514,file=paste(parafile,"Ndist",sep="."))
# save the plot
#dev.copy2pdf(file=paste(parafile,"_mmm.pdf",sep="") )
#return()
}
return()
}
