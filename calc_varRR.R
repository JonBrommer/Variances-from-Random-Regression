calc.varRR<-function(vvp.name, phi,K, positions=c(NA,NA)) {
	#function to calculate variances, the se in the variance and approximate 95% CI around the variance at
	#each value of a specified covariate, derived from random regression (RR) estimates
	#Fischer et al. 2004, Genet. Sel. Evol. 36:363–369
	# CURRENTLY ONLY RUNS FOR LINEAR (FORST-ORDER) POLYNOMIALS
	#function input:
	#vvp.name - can be
	#  1. name of asreml ".vvp" file with covariance matrix for estimates of elevation and slope (Average Information on variance component scale)
	#  2. a matrix with the vvp info (variances of REML-estimated variance in diagonal, covariances in off-diagonal
	#  3. an AsReml-R object of class "asreml" as produced by running the function asreml with random = ~us(pol(X,1)):SUBJECT~and rcov=~idv(units) or equivalent
	#phi - vector of (standardised) covariate values. note that these need to be on the same scale 
	#  if asreml or asreml-r was used, phi is assumed to be a scalar giving the number nodes on scale [-1,1]
	#  asreml scaling for covariate X is min(unique(X))- mean(unique(X)) and the standarised covariate [-1,1] can then be 
	#  multiplied with this scaling to obtain the data-scale covariate
	#K - covariance matrix of elevation and slope estimates by the RR. not needed when AsReml-R object is passed
	#positions - if vvp.name is a file, the positions of the uncertainty in the ".vvp" file
	#   positions=c(pos1,pos2) where [pos1,pos1] and [pos2,pos2] are the upper-left and bottm-right elements respectively
 	#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
	#initialise
	if (is.matrix(vvp.name)) { #asreml may not have been used
		vvp=vvp.name
		if (length(phi)==1) { 
			print("phi should be a vector giving the knot points")
			break()
		}
	} #is.matrix
	if (class(vvp.name)=="asreml") { #asreml used
		if (length(phi)>1) { 
			print("for asreml results, phi should be a scalar giving the number of knot points")
			break()
		}
		phi=seq(-1,1,length.out = phi)
		#construct vvp
		check=summary(vvp.name)$varcomp[,3]
		if (!is.na(sum(check))) {
			print("For this function to work, the asreml model must specify the residual error (use rcov=~idv(units) or equivalent)")
			break()	
		}
		vvp<-matrix(NA,length(check),length(check))
		r=rep(1:length(check),1:length(check))
		c=sequence(1:length(check))
		for (e in 1:length(vvp.name$ai)) {
			vvp[r[e],c[e]]<-vvp.name$ai[e]
			vvp[c[e],r[e]]<-vvp.name$ai[e]
		} #for
		#get rid of the empty row/column and name rows/cols
		nam<-row.names(summary(vvp.name)$varcomp)[!is.na(summary(vvp.name)$varcomp[,3])]
		vvp<-matrix(vvp[vvp!=0],sqrt(length(vvp[vvp!=0])),sqrt(length(vvp[vvp!=0]))
		, dimnames=list(nam,nam))
		if (length(grep("pol",nam))==0) {
			print("This function assumes Random Regression was performed using polynomials. Use pol() in asreml")
			break()
		}
		vvp<-vvp[grep("pol",nam),grep("pol",nam)]
		Kv<-summary(vvp.name)$varcomp[grep("pol",nam),1]
		K<-matrix(c(Kv[1],Kv[2],Kv[2],Kv[3]),2,2)
	} #if#asreml object
	if (is.character(vvp.name)) {
		#check of phi
		if (length(phi)>1) { 
			print("for asreml results, phi should be a scalar giving the number of knot points")
			break()
		}
		#get vvp
		vvp.1<-read.table(vvp.name,skip=1,sep="\t",stringsAsFactor=F)
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
