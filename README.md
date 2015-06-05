# Variances-from-Random-Regression
Calculating covariate specific variances and their SE on the basis of random regression estimates
Approach described by Fischer et al. 2004, Genet. Sel. Evol. 36:363–369
This function plots variances and approximate 95% CI based on double the SE. The function currently works only with RR of first-order polynomials

Input is either the vvp matrix itself (Hessian) or an AsReml file and a vector indicating the position of the relevant vvp elements

# example 1 
```
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
# (not run)
vvp.filename="C:\\vvp_example.txt"
K=matrix(c(5,1,1,1),2,2)
phi=seq(-1,1,by=(1/3))
calc.varRR(phi,K, vvp.filename, positions=c(3,5))
```
