cat("Reading TestUtils.R\n")
flush.console()


timeEquilibration <- function(nIter){
  
    Q <- matrix(runif(10000), ncol=100)  
    system.time(
      
        for(i in 1:nIter){
          
            equilibrate(Q)
            Q[50,50] <- runif(1)
        }
    )
}
# Test if equilibration improves the condition number of a symmetric matrix Q=A'A
# where A is a random 100x100 matrix with uniform entries in (-1,1) and dominant diagonal 
# modified so it has a bad condition number.
#
# Displays histogram of nIter trials, variable histogrammed is cond(Qp)/cond(Q),
# where Qp is the matrix Q equilibrated.
#
testEquilibration <- function(nIter){
  
    x <- rep(0,nIter)
    for(i in 1:nIter){
      
        A <- matrix(2*runif(10000)-1, ncol=100) 
        Q <- crossprod(A)
        
        #unbalance Q
        d <- (1:100)
        Q <- outer(d,d)*Q    # q <- diag(d)Q diag(d)
        
        eig <- eigen(Q,symmetric=TRUE)
        lam <- eig$values
        condQ <- max(lam)/min(lam)
        
        eq <- equilibrate(Q)
        Qp <- eq$.Qp
        eig <- eigen(Qp,symmetric=TRUE)
        lam <- eig$values
        condQp <- max(lam)/min(lam)
        
        x[i] <- condQp/condQ
    }
    q <- quantile(x,c(0.5,0.9,0.95,0.99))
    cat(paste("\nquantile(0.5,0.9,0.95,0.99):\n",q))
    hist(x,main="cond(Qp)/cond(Q)")
}




# Checking out how R::qr handles a rank-deficient mxn matrix A with m<n
# and rank(A)=p<m when computing the QR decomposition of the transpose t(A).
#
# The QR factorization t(A)=QR of the transpose is needed to compute the
# null space of A and special solutions to Ax=b.
#
show_QR  <- function(){

     msg <- "\n\nComputing QR factorization of transpose t(A)=QR where\n"
     msg <- paste(msg,"A is 5x7 with rank(A)=3:\n\n")
     cat(msg)

     A <- matrix(c(
            1,2,1,1,0,0,1,
            1,1,2,1,0,0,1,
            1,1,1,2,0,0,1,
            2,3,3,2,0,0,2,
            2,2,3,3,0,0,2
          ), nrow=5,byrow=TRUE
     )
     qrA <- qr(t(A))
     Q <- qr.Q(qrA,complete=FALSE)
     cat("\nA is 5x7 with rank 3, t(A)=QR: matrix Q (complete=FALSE):\n")
     print(Q)
     Q <- qr.Q(qrA,complete=TRUE)
     cat("\nA is 5x7 with rank 3, t(A)=QR: matrix Q (complete=TRUE):\n")
     print(Q)
     R <- qr.R(qrA)
     cat("A is 5x7 with rank 3, t(A)=QR: matrix R:\n")
     print(R)
}



# Checking out how R::svd handles a rank-deficient mxn matrix A with m<n
# and rank(A)=p<m when computing the SVD-decomposition of A=UDV'.
#
# Displays the matrices U (5x5), D (5x7) and V (7x7).
#
#
show_SVD  <- function(){

     msg <- "\n\nComputing SVD-factorization A=UDV' where\n"
     msg <- paste(msg,"A is 3x7 with full rank(A)=3:\n\n")
     cat(msg)

     A <- matrix(c(
            1,2,1,1,0,0,1,
            1,1,2,1,0,0,1,
            1,1,1,2,0,0,1
          ), nrow=3,byrow=TRUE
     )
     svdA <- svd(A,nv=7)
     U <- svdA$u
     cat("\n---U:\n"); print(U)
     d <- svdA$d
     cat("\n---d:\n"); print(d)
     V <- svdA$v
     cat("\n---V:\n"); print(V)
}




# Testing the computation of the null-space matrix of a random
# mxn matrix A with uniform entries in (-1,1).
# We must have m<n.
#
testNullSpaceBasisMatrix <- function(nTests,m,n){

    if(m>=n){
   
         msg <- "testNullSpaceBasisMatrix: we must have m<n but"
         msg <- paste(msg, "m =",m,"and n =",n)
         stop(msg)
    }

    cat("\nTesting nullSpaceMatrix.\n")
    for(ii in 1:nTests){
   
        A <- matrix(2*runif(m*n)-1,nrow=m)
        # with QR-decomposition of A
        Q <- nullSpaceBasisMatrix(A)
        X <- A%*%Q                         # must be zero
        error <- sqrt(sum(X*X)/(m*n))
        cat("test ",ii,"QR error = ",error,"\n")
        
        # with SVD-decomposition of A
        Q <- nullSpaceBasisMatrix_SVD(A)
        X <- A%*%Q                         
        error <- sqrt(sum(X*X)/(m*n))
        cat("test ",ii,"SVD error = ",error,"\n")
    }
}


# Testing the computation of the solution-space Ax=b for a random
# mxn matrix A with uniform entries in (-1,1) satisfying m<n and
# having full rank m.
#
testSolutionSpace <- function(nTests,m,n){

    if(m>=n){

         msg <- "testSolutionSpace: we must have m<n but"
         msg <- paste(msg, "m =",m,"and n =",n)
         stop(msg)
    }

    cat("\nTesting nullSpaceMatrix.\n")
    for(ii in 1:nTests){

        A <- matrix(2*runif(m*n)-1,nrow=m)
        b <- runif(m)
       
        # with QR-decomposition of A
        sol <- solutionSpace(A,b)
        x0 <- sol$.x0
        Q <- sol$.Q

        error <- sum((A%*%x0-b)^2)
        cat("\nwith QR-decomposition x0 satisfies Ax0=b, error = ",error,"\n")

        X <- A%*%Q                     # must be zero
        error <- sqrt(sum(X*X)/(m*n))
        cat("with QR-decomposition A annihilates the column space of Q, error = ",error,"\n")
        
        
        # with QR-decomposition of A
        sol <- solutionSpace_SVD(A,b)
        x0 <- sol$.x0
        Q <- sol$.Q

        error <- sum((A%*%x0-b)^2)
        cat("\nwith SVD-decomposition x0 satisfies Ax0=b, error = ",error,"\n")

        X <- A%*%Q
        error <- sqrt(sum(X*X)/(m*n))
        cat("with SVD-decomposition A annihilates the column space of Q, error = ",error,"\n")
    }
}


# Checking if we have a numerical problem (blas? processor?)
# Just how symmetric is A'A for random A?
#
testSymmetric <- function(n,nTests){
  
    msg <- "\nComputing max(abs(Q-t(Q))), where Q=A'A,\n"
    msg <- paste(msg,"for",nTests,"random matrices A with uniform entries in (-1,1):\n")
    msg <- paste(msg, "||Q-t(Q)||_oo:\n")
    cat(msg)
    for(i in 1:nTests){
        A <- matrix(2*runif(n*n)-1,nrow=n)
        Q <- crossprod(A)
        cat(max(abs(Q-t(Q))),", ")
    }
}


