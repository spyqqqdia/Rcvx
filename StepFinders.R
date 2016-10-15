cat("Reading StepFinders.R\n")
flush.console()


# Freestanding functions to solve linear systems that occur in the 
# computation of the search direction in a variety of constrained
# and unconstrained optimization problems.
#
# Generally these problems solve the necessary KKT conditions for an extremum
# (which are nonlinear) by linearizing them and using an iterative approach
# (a generalization of Newton's root finding algorithm).
#
# The linear system which we treat here serve the purpose to compute the
# "full step" u from one iterate x to the next. The next iterate has the form 
# x+su with s\in(0,1]. The full step x -> x+u is usually not taken and the
# scalar s is computed by the backtracking line search algorithm.
#
# The step us is based on a quadratic Taylor approximation of the objective 
# function f centered at the current iterate x. The full step x -> x+u will
# hit the minimum only if f is quadratic and there are no constraints.
#
# In general one is better off by taking several shorter steps (if only to
# ensure that constraints are maintained).



# Solves the equation Hu=-grad for the Newton step u in unconstrained
# optimization. Assumes that H is positive definite.
#
# This is applied to H=hessian(f)(x), where f is the objective function.
# Positive definiteness of H then corresponds to strict convexity of f.
#
solveNewtonStep <- function(H,grad){

    n <- nrow(H)  
    # switch in globals.R
    if(useEquilibration){
    
        eql <- equilibrate(H)
        d <- eql$.d
        # transformation H -> DHD, grad <- d*grad, D=diag(d)
        DHD <- eql$.Qp
    } else {
      
        d <- rep(1,n)
        DHD <- H
    }
    v <- -d*grad

    tryCatch({

        # H=R'R, R upper triangular
        # if no error, H is pos. definite
        R <- chol(DHD)
        L <- t(R)
    },
    error = function(e) stop("solveNewton: H not positive definite.")
    )
    y <- forwardsolve(L,v)
    x <- backsolve(R,y)
    d*x
}
qrSolveNewtonStep <- function(H,grad) qr.solve(H,-grad)








# @return the KKT matrix
#                    / H, A'\
#                    \ A, 0 /
#
kktMatrix <- function(H,A){

    Zeros <- matrix(0,nrow(A),nrow(A))
    rbind(cbind(H,t(A)),cbind(A,Zeros))
}


# Checks if the dimensions of H,A,grad,r in the KKT system
#
#      Hu+A'\nu = -grad
#      Au       = r
# are compatible.
#
checkDimensionsKKT <- function(H,A,grad,r){

    n <- ncol(H)
    if(nrow(H)!=n){

         msg <- "\ncheckDimensionsKKT: H must be square but:\n"
         msg <- paste(msg,"nrow(K) =",nrow(K),"ncol(K) =",ncol(K))
         stop(msg)
    }
    if(ncol(A)!=n){

         msg <- "\ncheckDimensionsKKT: we must have ncol(A)=ncol(H) but:\n"
         msg <- paste(msg,"nrow(A) =",ncol(A),"ncol(H) =",ncol(H))
         stop(msg)
    }
    if(length(grad)!=n){

         msg <- "\ncheckDimensionsKKT: we must have length(grad)=ncol(H) but:\n"
         msg <- paste(msg,"length(grad) =",length(grad),"ncol(H) =",ncol(H))
         stop(msg)
    }
    p <- nrow(A)
    if(length(r)!=p){

         msg <- "\ncheckDimensionsKKT: we must have length(r)=nrow(A) but:\n"
         msg <- paste(msg,"length(r) =",length(r),"nrow(A) =",nrow(A))
         stop(msg)
    }
}


# Solves a linearized KKT system with only equality constraints
# (such as occurs e.g. in the barrier method):
#
#      Hu+A'\nu = -grad
#      Au       = r
#
# This system of equations has the form Bu = v, where B is the KKT matrix
#                    / H, A'\
#                    \ A, 0 /
# and v=c(-grad,r).
#
# This is used to compute the full step u (before backtracking) from one
# iteration x to the next iteration x+su when computing the solution of the
# nonlinear KKT system. In this case H=hessian(f)(x), grad = grad(f)(x) and
# r is the residual r=b-Ax.
#
# The solution uses block elimination.
# The matrix H is first balanced: H --> DHD, where D=diag(d) is computed with
# the Ruiz algorithm.
# If DHD (equivalently H) is not positive definite, it is augmented via
# H -> Hp = H+rho*A'A.
#
# @return solution w=c(x,\nu).
#
solveKKTStep_NoIneq <- function(H,A,grad,r,rho=1){

    checkDimensionsKKT(H,A,grad,r)
    n <- nrow(H)

    # switch in globals.R
    if(useEquilibration){
    
        eql <- equilibrate(H)
        d <- eql$.d
        # transformation H -> DHD, grad <- d*grad, A <- AD, D=diag(d)
        DHD <- eql$.Qp
    } else {
      
        d <- rep(1,n)
        DHD <- H
    }
    
    v <- d*grad
    p <- nrow(A)
    A <-  outer(rep(1,p),d)*A
   
    tryCatch({

        # H=R'R, R upper triangular
        # if no error, H is pos. definite
        R <- chol(DHD)
        L <- t(R)
    },
    error = function(e){

        # DHD was not positive definite
        # transformation DHD -> DHD+rA'A
        DHD <- DHD+rho*crossprod(A)
        v <- v-rho*t(A)%*%r
        R <- chol(DHD)
        L <- t(R)
    })

    rhs <- cbind(t(A),v)
    Y <- forwardsolve(L,rhs)
    X <- backsolve(t(L),Y)
    H_invAt <- X[,1:ncol(t(A))]
    AH_inv <- t(H_invAt)
    H_inv_q <- X[,1+ncol(t(A))]

    w <- -(r+AH_inv%*%v)
    S <- A%*%H_invAt
    R <- chol(S)                  # upper triangular  S=R'R
    mu <- forwardsolve(t(R),w)
    nu <- backsolve(R,mu)

    # solve Hx=LL'x=-(q+A'nu)
    y <- forwardsolve(L,-(v+t(A)%*%nu))
    x <- backsolve(t(L),y)
    c(d*x,nu)
}


# Solves a KKT system with only equality constraints (such as occurs e.g.
# in the barrier method):
#
#      Hx+A'\nu = -grad
#      Ax       = r
#
# The QR factorization on the entire KKT matrix is used.  See the previous
# function for more details.
#
# @return solution w=c(x,\nu).
#
qrSolveKKTStep_NoIneq <- function(H,A,grad,r){

    K <- kktMatrix(H,A)
    qr.solve(K,c(-grad,r))
}
