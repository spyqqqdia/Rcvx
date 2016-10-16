source("globals.R")
source("StepFinders.R")
source("Solvers.R")
source("TestSolvers.R")

# setwd("C:/Projekte/R/KKT")
# source("run/runKKT.R")


# How accurate are the Newton steps
#
if(FALSE){

    testKKTSolvers(5)
}


# How much slower is QR-decomposition than block elimination with
# Cholesky factorization.
#
if(TRUE){

    n <- 500                   #  number of decision variables
    nIter <- 100
    timeKKTStepFinders(n,nIter)
}


# Do unconstrained Newton on a convex quadratic function with known
# minimum at zero.
# In principle this can be solved exactly in one Newton step (since
# the Newton step computation solves the quadratic Taylor approximation 
# exactly) but our backtracking line search algorithm uses the linear
# decrement and so may choose a shorter step (unclear to the author). 
#
if(TRUE){

    # backtracking line search params
    alpha <- 0.05
    bta <- 0.5
    maxIter <- 10
   
    n <- 1000           # dimension
    # random matrix H will be unbalanced via H -> DHD, D=diag(d)
    d_max <- 100
    delta <- (d_max-1)/(n-1)
    d <- 1.0+(0:(n-1))*delta    # linearly from 1 to d_max
   
    # terminate as soon as rate of decrease in search direction
    # is < 2*eps
    eps <- 1e-18
    testSolveNewton(n,d,maxIter,alpha,bta,eps)
}
