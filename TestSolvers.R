cat("Reading TestSolvers.R")
flush.console()



# Testing how accurately the systems for finding the Newton step are solved
# (Cholesy factorization versus QR decomposition).
#
testKKTSolvers <- function(nTests){

    msg <- "\n\n\n#----Solving _without_ equilibration----#\n"
    # switch in globals.R
    if(useEquilibration) msg <- "\n\n\n#----Solving _with_ equilibration----#\n"
    cat(msg)
  
    rho <- 1
    n <- 100
    B <- matrix(2*runif(n*n)-1,nrow=n)
    H <- crossprod(B)     # nxn
    #unbalance H
    d <- 1:100
    H <- outer(d,d)*H    # H <- diag(d)H diag(e)
    A <- matrix(runif(500),nrow=5)
   
    for(ii in 1:nTests){

        grad <- runif(n)
        r <- runif(5)
        

        cat("\n\nKKT system, no constraints:\n")
        w <- solveNewtonStep(H,grad)
        error <- sqrt(sum((H%*%w+grad)^2))
        cat("Cholesky factorization, error = ",error)
        w <- qrSolveNewtonStep(H,grad)
        error <- sqrt(sum((H%*%w+grad)^2))
        cat("\nqr.solve, error = ",error)

        cat("\n\nKKT system, no inequalities:\n")
        z <- c(-grad,r)
        w <- solveKKTStep_NoIneq(H,A,grad,r,rho)
        K <- kktMatrix(H,A)
        error <- sqrt(sum((K%*%w-z)^2))
        cat("block elimination, error = ",error)
        w <- qrSolveKKTStep_NoIneq(H,A,grad,r)
        error <- sqrt(sum((K%*%w-z)^2))
        cat("\nqr.solve, error = ",error)
    }
}




timeReport <- function(nIter,operation,opID){

    timeInfo <- system.time(for(ii in 1:nIter)operation())
    cat("\n",opID,":\n")
    print(timeInfo)
}
# Time the solution of the various systems for n decision variables to see
# how much slower the simpler algo using the QR-decomposition is.
# The number p of equality constraints is comparably small, we use p=5.
#
timeKKTStepFinders <- function(n,nIter){

    cat(
        "\n\nTiming solution of ",nIter,
        " systems of equations in dimension ",n,":\n"
    )
    rho <- 1
    B <- matrix(runif(n*n),nrow=n)
    H <- crossprod(B)
    A <- matrix(runif(5*n),nrow=5)    # 5 equality constraints

    grad <- runif(n)
    r <- runif(5)

    operation <- function(){  solveNewtonStep(H,grad); grad[2] <- runif(1) }
    opID <- "Newton with Cholesky"
    timeReport(nIter,operation,opID)
   
    operation <- function(){  qrSolveNewtonStep(H,grad); grad[2] <- runif(1) }
    opID <- "Newton step with QR"
    timeReport(nIter,operation,opID)
   
    operation <- function(){  solveKKTStep_NoIneq(H,A,grad,r,rho); grad[2] <- runif(1) }
    opID <- "KKTStep_NoIneq step with block elimination"
    timeReport(nIter,operation,opID)
   
    operation <- function(){  qrSolveKKTStep_NoIneq(H,A,grad,r); grad[2] <- runif(1) }
    opID <- "KKTStep_NoIneq step with QR"
    timeReport(nIter,operation,opID)
}

# We test the basic Newton solver (solveNewton in Solvers.R) on the function
# f(x)=||Ax||^2 with condition ||x||^2<1.
#
# As A we choose the matrix A=DMD, where M is an nxn random matrix with entries
# in (0,1) and D=diag(d) a diagonal matrix with positive entries. By choosing
# the d_i of very uneven sizes (e.g. growing like d_i=i) we can make the
# matrix A very poorly conditioned.
#
# @param alpha,bta the backtracking line search parameters.
# @param b vector with positive entries.
# @param eps: termination criterion is |u.grad|<2*eps, where u is the
# full Newton step (this is the maximal possible decrease after linearization
# in this direction).
#
testSolveNewton <-  function(n,d,maxIter,alpha,bta,eps){

    cat("\ntestSolveNewton: solver started.\n")
    M <- matrix(runif(n*n),nrow=n)
    didj <- outer(d,d)
    A <- didj*M
    hess <- 2*crossprod(A)

    f <- function(x){ v <- A%*%x; sum(v*v) }
    grad <- function(x){
         v <- A%*%x; 2*t(A)%*%v
    }
    H <- function(x) hess
    condition <- function(x) sum(x*x)<1
   
    #starting point
    x0 <- rep(0.9999/sqrt(n),n)
   
    stopifnot(condition(x0))
    x <- solveNewton(x0,f,grad,H,condition,alpha,bta,maxIter,eps)
   
    cat("Solution x should be zero and is:\n")
    print(x)
}

