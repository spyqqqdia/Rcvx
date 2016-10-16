cat("Reading Utils.R\n")
flush.console()


#clear the console
cls <- function() cat("\014")

# Equilibrates the _symmetric_ square matrix Q by multiplication on left and right by 
# the same diagonal matrix D. This preserves symmetry.
# Uses the algorithm of Ruiz with one round of ||.||_oo equilibration
# followed by 5 rounds of ||.||_2 equlibration.
#
#@return list(.d,.e,.Qp), where .d, .e is the diagonal matrix and Qp=DQE, 
# where D=diag(.d) and E=diag(.e).
#
equilibrate <- function(Q){
  
    if(class(Q)!="matrix") stop("\nequilibrate: argument must be a matrix")
    if(nrow(Q)!=ncol(Q)) stop(paste(
        "\nequilibrate: matrix Q not square, nrow(Q) = ",nrow(Q),"; ncol(Q) =",ncol(Q)
    ))
    # initialization of d,e
    n <- nrow(Q); d <- rep(1,n); e <- rep(1,n)
    # one round of ||.||_oo equilibration:
    norm_oo <- function(x){
        r <- max(abs(x)); ifelse(r>0,r,1)      # "1": do nothing with zero rows or cols.
    }
    d <- d/sqrt(apply(Q,1,norm_oo))            
    didj <- outer(d,d)                         # matrix (d_i*e_j)
    Qp <- didj*Q                               # diag(d)Qdiag(e)
    
    # 5 rounds of ||.||_2 equlibration
    norm_2  <- function(x){
        r <- sqrt(sum(x*x)/length(x)); ifelse(r>0,r,1)
    }
    for(i in 1:5){
        
        d <- d/sqrt(apply(Qp,1,norm_2))                  
        didj <- outer(d,d)
        Qp <- didj*Q   
    }
    list(.d=d,.Qp=Qp)
}



# In function callerID check that matrix A has more columns than rows.
#
checkHasMoreColsThanRows <- function(A,callerID){

    if(nrow(A)>=ncol(A)){

        msg <- paste(callerID,": we must have nrow(A) < ncol(A) but\n",sep="")
        msg <- paste(msg,"nrow(A)=",nrow(A),"and ncol(A)=",ncol(A),"\n")
        stop(msg)
    }
}


# Function computes the null space of an mxn matrix A with m<n and rank p
# using the QR-decomposition of A.
# See doc/nullspace, section 1.1, p1.
# @returns orthonormal nx(n-p) matrix F such that ker(A)=Im(F)=colspace(F).
#
nullSpaceBasisMatrix <- function(A){

    checkHasMoreColsThanRows(A,"nullSpaceBasisMatrix")

    qrA <- qr(t(A))
    Q <- qr.Q(qrA,complete=TRUE)
    p <- qrA$rank
    Q[,(p+1):n]
}

# Computes the space of solutions of Ax=b as a pair (w0,F), where w0 is a
# particular solution and the columns of the matrix F are an orthogonal basis
# for the kernel ker(A): Im(F)=colspace(F)=ker(A).
#
# Then the set of all solutions is the hyperplane
#          w0+Im(F) = w0+span(columns of F)
#
# A must be an mxn matrix with m<n and full rank(A)=m. In this case a solution
# is guaranteed to exist.
# This function uses the QR-decomposition of A.
# See doc/nullspace, section 1.1, p1.
#
# @ return list(.w0,.F).
#
solutionSpace <- function(A,b){

    checkHasMoreColsThanRows(A,"solutionSpace")
    n <- ncol(A)
    qrAt <- qr(t(A))
    Q <- qr.Q(qrAt,complete=TRUE)
    R <- qr.R(qrAt)        # incomplete, no zero rows on bottom

    # note A=R'Q' and so solve R'Q'x=b, set y=Q'x, solve R'y=b,...
    y0 <- forwardsolve(t(R),b)
    P <- Q[,1:m]
    w0 <- P%*%y0

    list(.w0=w0,.F=Q[,(m+1):n])
}



# Function computes the null space of an mxn matrix A with m<n and rank p
# using the SVD-decomposition of A.
# # See doc/nullspace, section 1.2, p3.
# @returns orthonormal nx(n-p) matrix F such that ker(A)=Im(F)=colspace(F).
#
nullSpaceBasisMatrix_SVD <- function(A){

    m <- nrow(A)
    n <- ncol(A)
    checkHasMoreColsThanRows(A,"solutionSpace_SVD")
    svdA <- svd(A,nv=n)              # A=UDV'

    Vp <- svdA$v                     # Vp=[V,v_{m+1},...,v_n]
    Vp[,(m+1):n]
}

# Computes the space of solutions of Ax=b as a pair (w0,F), where w0 is a
# particular solution and the columns of the matrix F are an orthogonal basis
# for the kernel ker(A): Im(F)=colspace(F)=ker(A).
#
# Then the set of all solutions is the hyperplane
#          w0+Im(F) = w0+span(columns of F)
#
# A must be an mxn matrix with m<n and full rank(A)=m. In this case a solution
# is guaranteed to exist.
# This function uses the SVD-decomposition of A.
# See doc/nullspace, section 1.2, p3.
#
# @ return list(.w0,.F).
#
solutionSpace_SVD <- function(A,b){

    m <- nrow(A)
    n <- ncol(A)
    checkHasMoreColsThanRows(A,"solutionSpace_SVD")
    svdA <- svd(A,nv=n)              # A=UDV'

    U <- svdA$u
    Vp <- svdA$v
    V <- Vp[,1:m]
    d <- svdA$d                 # D=diag(d)
   
    vec_c <- t(U)%*%b           
    w <- vec_c/d
    w0 <- V%*%w

    list(.w0=w0,.F=Vp[,(m+1):n])
}

# Given a solution x0 of Ax=b (underdetermined, full rank) as well as a 
# parametrization of the solutions as x=w0+Fu where F has orthonormal
# columns
#
# @return a parameter u0 such that x0 = w0+Fu0.
#
solutionSpaceParameter <- function(x0,w0,F){ t(F)%*%(x0-w0) }




