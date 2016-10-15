source("Utils.R")
source("TestUtils.R")


# setwd(/home/oar/Projects/R/Rcvx)
# source("runUtils.R))


# Print ||Q-t(Q)||_oo, where Q=A'A and A a random nxn-matrix with
# entries in (-1,1). 
# This is a test if we have a hardware problem or libopenblas has 
# a problem.
#
if(TRUE){
  
   n <- 200
   nTests <- 10
   testSymmetric(n,nTests)
  
}



if(FALSE){
  
    timeEquilibration(1000)
}
# Allocates 100x100 symmetric random matrices Q with entries in (0,1).
# Then unbalances the diagonal (to become dominant) so that the matrix
# gets a large condition number.
#
# Then equilibrates the matrix and computes the condition number of the
# equilibrated matrix Qp.
# Collects the fractions cond(Qp)/cond(Q) in a vector x and reports the
# histogram of x as well as the quantiles quantile q(0.5,0.9,0.95,0.99).
#
# In this way we see how much equilibration reduces (or does not reduce)
# the condition number.
#
if(FALSE){
  
  testEquilibration(1000)
}



# Displays the factor matrices Q,R in the QR decomposition of a rank-deficient
# 5x7 matrix A with rank(A)=3.
# Q is computed both with
#                            Q=qr.Q(qrA,complete=FALSE) and
#                            Q=qr.Q(qrA,complete=TRUE)
# where qrA=qr(A).
#
if(FALSE){

    show_QR()
    show_SVD()
}




# Test if the computed nullspace is really annihilated by A.
# Here A is a random mxn matrix with entries uniform in (-1,1).
#
if(FALSE){

     m <-  100
     n <- 1000
     nTests <- 5
     testNullSpaceBasisMatrix(nTests,m,n)
}



# Test if the computed solution space list(.x0,.Q) of Ax=b satisfies
#
#     Ax0=b and
#     AQ=0 (A annihilates the column space of Q, i,e, Im(Q)\subseteq ker(A)
#
# Here A is a random mxn matrix with entries uniform in (-1,1) and b a random
# vector of length m.
#
if(TRUE){

     m <-  100
     n <- 1000
     nTests <- 5
     testSolutionSpace(nTests,m,n)
}



