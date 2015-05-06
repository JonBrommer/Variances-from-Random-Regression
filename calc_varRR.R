calc.varRR<-function(phi,K, vvp.filename, positions=c(NA,NA)) {
	#function to calculate variances, the se in the variance and approximate 95% CI around the variance at
	#each value of a specified covariate, derived from random regression (RR) estimates
	#Fischer et al. 2004, Genet. Sel. Evol. 36:363–369
	#function input:
	#phi - vector of (standardised) covariate values. note that these need to be on the same scale 
	#as the covariate used in the Random Regression
	#K - covariance matrix of elevation and slope estimates by the RR
	#vvp.filename - name of file with covariance matrix for estimates of elevation and slope (Hessian or Average Information criterion)
	# if a matrix is passed it is assumed to be the vvp info
	#positions - vector with the positions of the uncertainty in the RR, not to be provided if vvp is matrix
 	#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
	#initialise
	if (is.matrix(vvp.filename)) {
		vvp=vvp.filename
	} else {
		vvp.1<-read.table(vvp.filename,skip=1,sep="\t",stringsAsFactor=F)
		vec.vvp<-as.vector(unlist(lapply(vvp.1,  function(m) strsplit(m," "))))
		vec.vvp<-as.numeric(vec.vvp[vec.vvp!=""])
		# unelegant brute-force finding the size of the matrix
		m.size <- function(x,n) {(x-(n*(n+1)/2))^2}
		n=round(optimize(m.size,c(0,100),x=length(vec.vvp))$min)
		m<-matrix(0,n,n)
		#below places the lements in the sub-diagonal by row
		teller=1
		for (i in 1:n) {
			for (j in 1:i) {
				m[i,j]=vec.vvp[teller]
				m[j,i]=vec.vvp[teller]
				teller=teller+1
			} #j
		} #i
		vvp=m[positions[1]:positions[2],positions[1]:positions[2]]
	} #if
	phi<-matrix(c(rep(1,length(phi)),phi),length(phi),2) #add the constant to phi
	P <- phi%*%K%*%t(phi) #variances at each value of phi in diagonal
	#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
	#below vectorises the variances in each element of K (i.e. uncertainty) and covariances between these uncertainties
	#this vectorised form is the trick needed to estimate the uncertainty following delta method
	#covariances around the estimates in K are either provided as a vvp matrix 
	#for (in diagonal): elevation (e), cov(e,s), slope (s)
	#OR, they are provide[d through the ".vvp" output of AsReml in which case the filename of this file
	#and the positions of the elements in this file are provided as "positions" (not implemented yet)
	#the vectorised form equations require 
	k2=matrix(NA,4,4)
	mapping=c(1,2,2,3)
	for (row in 1:4) {
		for (col in 1:row) {
			k2[row,col]=vvp[mapping[row],mapping[col]]
			k2[col,row]=vvp[mapping[row],mapping[col]]
		} #col
	} # row
	#back to calculating
	kron.phi=phi%x%phi
	var.var.mat=kron.phi%*%k2%*%t(kron.phi) # this is an (n x n) x (n x n) matrix
	var.var.vec=diag(var.var.mat)
	se=matrix(var.var.vec,dim(phi)[1],dim(phi)[1])
	ase=sqrt(diag(se));
	var.and.se.ci<-cbind(phi[,2],diag(P),ase,diag(P)-2*ase,diag(P)+2*ase) 
	#above cbinds covariate, variance, SE, and approximate CI per covariate value
	dimnames(var.and.se.ci)[[2]]<-c("cov.value", "variance", "se.variance","lower.ci","upper.ci")
	min.max=c(min(var.and.se.ci),max(var.and.se.ci))
	#function produces a plot to check output makes sense
	plot(var.and.se.ci[,2]~phi[,2],type="l",ylim=min.max,ylab="Variance",xlab="Standardised covariate")
	par(new=T)
	plot(var.and.se.ci[,4]~phi[,2],type="l",lty="dashed",ylim=min.max,ylab=NA,xlab=NA)
	par(new=T)
	plot(var.and.se.ci[,5]~phi[,2],type="l",lty="dashed",ylim=min.max,ylab=NA,xlab=NA)
	return(var.and.se.ci)
} #function

#example 
rm(list=ls())
K=matrix(c(5,1,1,1),2,2)
phi=seq(-1,1,by=(1/3))
vvp=matrix(c(
0.160930E-01, -0.353929E-02, 0.778355E-02,
-0.353929E-02,  0.779459E-03, -0.171652E-02  ,
0.778355E-02, -0.171652E-02,  0.378527E-02),3,3)
calc.varRR(phi,K,vvp)
# using vvp file as input
vvp.filename="C:\\Users\\joegbr\\Documents\\data\\Tutkimus\\projects\\MethodsPapers\\vvp_example.txt"
calc.varRR(phi,K, vvp.filename, positions=c(3,6))
