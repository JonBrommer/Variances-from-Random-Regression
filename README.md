# Variances-from-Random-Regression
calculating covariate specific variances and their SE on the basis of random regression estimates

function plots variances and approximate 95% CI based on double the SE

This currently works only with RR of first-order polynomials

input is either the vvp matrix itself or an AsReml file and a vector indicating the position of the relevant vvp elements

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
vvp.filename="C:\\vvp_example.txt"
calc.varRR(phi,K, vvp.filename, positions=c(3,5))

