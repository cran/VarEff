#############################################################################################
#  Preliminary analysis of data and building a parameter file to run VarEff
#    input : name of the data file  (infile)
#            the number of markers  (parafile)
#            a name for the current parameters (NBLOC)
#
#    output: two files: parafile.theta , three global estimates of theta
#                       parafile.R     , the command to run VarEff
#
Theta <-
function(infile=NULL,
                   parafile=NULL,
                         NBLOC=NULL, 
                         MUTAT=NULL,
                         JMAX=NULL,
                         MODEL=NULL,
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
print(" Estimation of Theta ",quote="FALSE")
print(" from microsatellites data",quote="FALSE")
print(" ",quote="FALSE")
print(" *** Data *** ",quote="FALSE")
if(!is.null(infile)){
print (paste(" Name of the data file= ",infile,sep=""),quote="FALSE")}else{
infile = readline(" Name of the data file= ")}
if(!is.null(NBLOC)){
print (paste(" Number of microsatellite markers= ",NBLOC,sep=""),quote="FALSE")}else{ 
NBLOC = as.numeric(readline(" Number of microsatellite markers= "))}
# If not given, give a name for the job
if(!is.null(parafile)){
print(paste(" Name of the job = ",parafile,sep=""),quote="FALSE")}else{
parafile = readline(" Give a name for this job = ")}
# Read the data from infile
#
##############################################################################
#
# read data and extract allele frequencies
#           skipping loci with 0 sample size
#
# input: infile, data file
#          NBLOC, number of markers
#
# output : NbMaxAll, PkMoyTout, VarCovTout
#          and the effective number of loci : NELOC
#
.DataExtract <- function(infile,NBLOC)
{
# Read the data file
#
TOUT <- scan(infile)
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
# Second step : Extracting frequencies
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
		locut = locut + 1
		if ( nbgenes < 10 ) {
 		print(c(" Marker ",locus," , low sample size = ",nbgenes), quote="FALSE")
					   }
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
#############################################
datex = .DataExtract(infile,NBLOC)
PkMoy  = datex$P
VarCov = datex$V
NbMaxAll   = datex$A
NELOC      = datex$L
print(" ************************************************************************",quote="FALSE")
print(" ",quote="FALSE")
print(" *** Priors about population sizes ***",quote="FALSE")
print(" ",quote="FALSE")
print(" Preliminary analysis of data to choose Ne priors (current and ancestral) ",quote="FALSE")
print(" thanks to Theta estimate from microsatellites data.",quote="FALSE")
print(" ",quote="FALSE")
print(" VarEff calculations make use of normalized scales for population size and times,",quote="FALSE")
print(" so that population sizes are expresssed as theta values, Theta = 4Neu,",quote="FALSE")
print(" where Ne is the census population size and u the global mutation rate.",quote="FALSE")
print(" Choosing a mutation rate allows conversion from these normalized scale",quote="FALSE")
print(" to natural values of population size Ne.",quote="FALSE")
#
if(!is.null(MUTAT)){
print (paste(" Mutation rate MUTAT= ",MUTAT,sep=""),quote="FALSE")}else{
MUTAT <- as.numeric(readline(" Give the mutation rate MUTAT = "))} 
#
# Estimation of theta_0, theta_1 et theta_2
#
thetasim <- vector(mode="numeric",length=3)
#
# Theta_1
theta1 = 2 * sum(c(1:(NbMaxAll-1))*PkMoy[2:NbMaxAll])
theta1 = theta1*(theta1+sqrt(theta1**2+1))
thetasim[2]=theta1
# Theta_2
theta2 = 2*sum((c(1:(NbMaxAll-1))**2)*PkMoy[2:NbMaxAll])
thetasim[3]= theta2
# the variance of theta_2 estimates accross markers
#squares = c(1:(NbMaxAll-1))^2
#vartheta2O = t(squares) %*% VarCov[2:NbMaxAll,2:NbMaxAll] %*% squares
# its theoretical expectation for independent markers
#vartheta2T = theta2 * ( 4 * theta2 + 1 ) / 3
#vartheta2T = sqrt(vartheta2T)
#vartheta2O = 2 * sqrt(vartheta2O)
# Theta0
theta0 = PkMoy[1]
theta0 = theta0**2
theta0 = (1/theta0)-1
theta0 = theta0/2
thetasim[1]=theta0
#
# Global values of Theta_0, Theta_1, Theta_2 
#
print(" Theta (4*Ne*u) is estimated using 3 estimators correlated to",quote="FALSE")
print(" Present, Intermediate and Ancestral population size.",quote="FALSE")
print(" Results are written to a .Theta file ",quote="FALSE")
print(" ",quote="FALSE")
print(" Mean values of estimates Theta_0, Theta_1, Theta_2",quote="FALSE")
print(thetasim,quote="FALSE")
#
print(" Imbalance indices ln(Theta_1/Theta_0) and ln(Theta_2/Theta_0)",quote="FALSE")
ib1 = log(theta1/theta0) 
ib2 = log(theta2/theta0) 
print(c( ib1 , ib2 ),quote="FALSE")
#
write(thetasim[1:3],file=paste(parafile,"Theta",sep="."))
write(c(ib1,ib2),file=paste(parafile,"Theta",sep="."),append="TRUE")
#
print(" ",quote="FALSE")
print(" The prior value of population size (Ne) and the assumed mutation rate (MUTAT)",quote="FALSE")
print(" must be in agreement with these estimated theta values."  ,quote="FALSE")
#
MinNe = as.integer(min(thetasim)/(4*MUTAT))
MaxNe = as.integer(max(thetasim)/(4*MUTAT))
print(" ",quote="FALSE")
print(paste(" Hence, choose the prior Ne in the range of ",MinNe," to ",MaxNe) ,quote="FALSE")
#
write(c(MinNe,MaxNe),file=paste(parafile,"Theta",sep="."),append="TRUE")
#print(" ",quote="FALSE")
#print(" Results of Theta, Imbalance indices, and MinNe and MaxNe are written to a .Theta file",quote="FALSE")
#print(" ",quote="FALSE")
#print(" Thanks to the MinNe and MaxNe, you can made a choice on census Ne (NBAR), ",quote="FALSE")
#print(" which must be given to execute VarEff function. ",quote="FALSE")
#print(" The prior value of census population size (NBAR) should be between MinNe and MaxNe value.",quote="FALSE")
print(" ",quote="FALSE")
#
#########################################################################
#
print(" ",quote="FALSE")
print(" ************************************************************************",quote="FALSE")
print(" ",quote="FALSE")
print(" *** Set parameters for the VarEff algorithm ***",quote="FALSE")
#
# Ask for parameters to run VarEff:
print(" ",quote="FALSE")
print(" *** Maximal distance between alleles (DMAX) *** ",quote="FALSE")
print(" ",quote="FALSE")
print(" DMAX is the maximal distance between alleles in terms of repeat units ",quote="FALSE")
print(" that you have to provide in VarEff function ",quote="FALSE")
print(" ",quote="FALSE")
print(" Look at the distribution of the f_k frequencies of pairs of alleles",quote="FALSE")
print(" at distance k found in your microsatellites markers data,",quote="FALSE")
print(" and take a large enough maximal allele distance DMAX, ",quote="FALSE")
print(" ",quote="FALSE")
print(" A suggestion is to take DMAX such that frequencies f_k",quote="FALSE")
print(" are under the red line for k > DMAX (i.e. f_k < 0.005 for k > DMAX)",quote="FALSE")
print(" (this choice is generally efficient) ",quote="FALSE")
#
Dk = length(PkMoy)
abscissa = " Distance k between alleles"
ordinates= " Frequencies f_k of pairs of alleles at distance k"
title = paste(" Data file: ", infile)
par(adj=0.5)
plot(c(-0.5,Dk+0.5),c(0,max(PkMoy)),xlab=abscissa,ylab=ordinates,main=title,type='n')
lines(-0.5+c(0:(Dk-1)),PkMoy,type='h')
esfplus = c(PkMoy,0)
lines(-0.5+c(0:Dk),esfplus,type='s')
lines(c(-0.5,Dk-0.5),c(0.005,0.005),col="red")
#
# Choose the DMAXPLUS
#
if(!is.null(DMAXPLUS)){
print (paste(" Distance between the alleles in term of repeat units= ",DMAXPLUS,sep=""),quote="FALSE")}else{
DMAXPLUS <- 1 + as.numeric(readline(" Give the DMAX = "))}
#
# check that DMAXPLUS is smaller than the largest allele distance
if( DMAXPLUS > NbMaxAll )
{
print(" ERROR !!!",quote="FALSE")
print(c(" DMAX must be smaller than ",NbMaxAll),quote="FALSE")
DMAXPLUS = 1 + as.numeric(readline(" re-enter DMAX : "))
}
# final check !
if( DMAXPLUS > NbMaxAll )
{
print(c("DMAX is reduced to ",(NbMaxAll-1)),quote="FALSE")
DMAXPLUS = NbMaxAll-1
}
#
print(" ",quote="FALSE")
print(" *** Demographic model (JMAX) *** ",quote="FALSE")
print(" Choose the number JMAX of transitions between periods with constant Ne. ",quote="FALSE")
print(" This value is an integer generally between 2 and 6. ",quote="FALSE")
if(!is.null(JMAX)){
print (paste(" Number of transitions =",JMAX,sep=""),quote="FALSE")}else{
JMAX <- as.numeric(readline(" Value of JMAX = "))}
#
print(" ",quote="FALSE")
print(" *** Mutation model (MODEL) *** ",quote="FALSE")
#
print(" ",quote="FALSE")
print(" Choose the mutation model (S, T or G) ",quote="FALSE")
print("   (S = Single Step Model, T = Two Phase Model, G = Geometric Model)",quote="FALSE")
print("   and the prior value of its C parameter,  0 < C < 1",quote="FALSE")
print("   e.g. answer ' T   0.15 ' for the Two Phase Model with a proportion 15% of two motifs steps",quote="FALSE")
if(!is.null(MODEL)){
print (paste(" Mutation model =",MODEL,sep=""),quote="FALSE")}else{
MODEL = readline(" Mutation model and C coefficient? ")}
#
print(" ",quote="FALSE")
print(" *** Prior of census population size (NBAR) *** ",quote="FALSE")
print(" Look the two last values (MinNe and MaxNe) in the file .Theta",quote="FALSE")
print(" and choice NBAR value (integer) between them.",quote="FALSE")
#
if(!is.null(NBAR)){
print (paste(" Prior value of census Ne =",NBAR,sep=""),quote="FALSE")}else{
NBAR = as.numeric(readline(" Give prior of census Ne = "))}
#
print(" ",quote="FALSE")
print(" *** Prior variance of census population size in logarithm (VARP1) *** ",quote="FALSE")
print(" More you are not sure of your prior, more increase the variance.",quote="FALSE")
print(" Ex: VARP1=3 allows for searches with 20- to 40-fold relative variations of effective size.",quote="FALSE")
if(!is.null(VARP1)){
print (paste(" Prior variance of log(Ne) =",VARP1,sep=""),quote="FALSE")}else{
VARP1 = as.numeric(readline(" Give the prior variance of log(Ne) = "))}
#
print(" ",quote="FALSE")
print(" *** Prior Prior correlation between successive population sizes (RHOCORN) *** ",quote="FALSE")
print(" Generally RHOCORN=0, meaning that there are any correlation.",quote="FALSE")
if(!is.null(RHOCORN)){
print (paste(" Prior correlation between successive population sizes =",RHOCORN,sep=""),quote="FALSE")}else{
RHOCORN = as.numeric(readline(" Give the prior correlation between successive population sizes = "))}
#
print(" ",quote="FALSE")
print(" *** Priors about time and time intervals (GBAR) (number of generations)***",quote="FALSE")
print("Number of generations since the assumed origin of the population.",quote="FALSE")
if(!is.null(GBAR)){
print (paste(" Prior time since the Origin (number of generations) =",GBAR,sep=""),quote="FALSE")}else{
GBAR =as.numeric(readline(" Give the prior time since the Origin = "))}
# 
print(" ",quote="FALSE")
print(" *** Prior variance of the time since the origin (VARP2) *** ",quote="FALSE")
print(" More you are not sure of your prior, more increase the variance.",quote="FALSE")
if(!is.null(VARP2)){
print (paste(" Prior Variance of log(time interval) =",VARP2,sep=""),quote="FALSE")}else{
VARP2 = as.numeric(readline(" Give the Prior Variance of log(time interval) = "))}
#
print(" ",quote="FALSE")
#
# MCMC
print(" *** Smoothing the Covariance matrix (Diagonale) ***",quote="FALSE")
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
print (paste(" Covariance matrix = ",Diagonale,sep=""),quote="FALSE")}else{
Diagonale= as.numeric(readline(" Enter D value = "))}
#
Diagonale = min(Diagonale,1)
if( Diagonale <= 0) 
{
print(" D must be > 0, it is reset to 0.1",quote="FALSE")
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
#  save the VarEff file 
namecom = paste(parafile,".R",sep="")
print (paste(" Name of the parameter set = ",parafile,sep=""),quote="FALSE")
print(paste(" enter the command: source('",namecom,"') to launch the VarEff function",sep=""),quote="FALSE")
.parDiagwrite(parafile,infile,NBLOC,JMAX,MODEL,MUTAT,NBAR,VARP1,RHOCORN,GBAR,VARP2,DMAXPLUS,Diagonale,NumberBatch,LengthBatch,SpaceBatch)
#############################################
return()
}