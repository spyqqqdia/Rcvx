source("Constraint.R")
source("ConstraintSet.R")
source("TestConstraints.R")


# Mainly checking if affine change of variables in the constraints
# works correctly.
#
if(TRUE){
  
    n <- 20   # dimension
    d <- 10   # reduced dimension
    x0 <- 4*runif(n)-2
    F <- matrix(2*runif(n*d)-1,nrow=n,ncol=d)
    u <- 2*runif(d)-1
    basicConstraintSetTest(x0,F,u)
  
}
