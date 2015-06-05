# Variances-from-Random-Regression

# Description
Calculates covariate specific variances and their SE on the basis of random regression estimates
Approach described by Fischer et al. 2004, Genet. Sel. Evol. 36:363–369
This function plots variances and approximate 95% CI based on double the SE. The function currently works only with RR of first-order polynomials

# Usage
```
calc.varRR(K,phi,vpp)
calc.varRR(K,phi,vpp,positions)
```

# Arguments 

K - the covariance matrix of random effect for elevation (intercept) and slope. These random effects are here assumed to be polynomials (a + b*x). 

phi - the covariate which *must* be scaled in the same way as in the analysis. Note that your software may internally scale this. AsReml will have the minimal value of the covariate as -1 and the mean of the set of observed covariate values as 0. Thus, if your covariate has values 0, 1, 2 (each of which may occur different times in your data); AsReml will scale these to -1, 0, 1. Again, other softare will do this differently. R, for example, does not change the scaling.

vvp - either the vvp matrix itself (Hessian) or an AsReml file and a vector indicating the position of the relevant vvp elements. The vvp matrix denotes the variances (in diagonal) and covariances associated with the estimated Random Regression covariance matrix. 

positions - AsReml will generate a ".vvp" file. AsReml counts all estimated covariances in the order as reported in the ".asr" file. For example, the Random Regression variance in intercept, covariance between intercept and slope and variance in slope are the 3rd, 4th and 5th covariances reported, respectively. In that case, the relevant part of the vvp matrix concerns these positions. You can either extract these manually from the vvp file (example 1 below), or pass the whole file to calc.varRR and specify the positions (example 2 below).

# Use
Copy the script directly into your own script.
Alternatively, save it as a file to your working directory or other directory (in which case the path tto that directory needs to be included). 
```
calc.varRR<-dget(“calc.varRR.R”)
```
or, when you have it not in your working directory e.g.
```
calc.varRR<-dget(“C:\\MyScripts\\calc.varRR.R”)
```

# example 1 
```
calc.varRR<-dget(“calc.varRR.R”)
K=matrix(c(5,1,1,1),2,2)
phi=seq(-1,1,by=(1/3))
vvp=matrix(c(
0.160930E-01, -0.353929E-02, 0.778355E-02,
-0.353929E-02,  0.779459E-03, -0.171652E-02  ,
0.778355E-02, -0.171652E-02,  0.378527E-02),3,3)
calc.varRR(phi,K,vvp)
```

# example 2: as above, but using a vvp file as input
```
## Not run:
calc.varRR<-dget(“calc.varRR.R”)
vvp.filename="C:\\vvp_example.txt"
K=matrix(c(5,1,1,1),2,2)
phi=seq(-1,1,by=(1/3))
calc.varRR(phi,K, vvp.filename, positions=c(3,5))
```
