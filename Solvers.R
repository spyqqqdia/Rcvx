cat("Reading Solvers.R\n")
flush.console()


# Routines to solve the various nonlinear KKT systems by linearization and
# iteration. This uses the appropriate step finder to find the full Newton
# step and then one of the updaters to compute the next iteration until the
# iterations have converged.
#



# Routine which update one iterative solution x to one of the nonlinear
# KKT systems to the next x+s*u by starting with s=1 and backtracking
# s -> bta*s until condition(x+s*u) returns TRUE _and_ the descent condition
#          
#            f(x+su) <= f(x)+alpha*s*grad(x)'u  
# is satisfied.
#
# Stops with error if the constraint is not satisfied
# either at x or after 50 iterations.
#
# @param alpha: must be in (0,1), say alpha=0.05.
# @param bta: must be in (0,1), ideally = 0.5.
# @param fx: value f(x)
# @param grad: value of gradient grad(f)(x) at the point x.
# @return x+s*u
#
backtrack <- function(x,u,f,grad,condition,alpha,bta){

    if(!condition(x))
        stop("\nbacktrack: starting point x does not satisfy the condition.\n")
    iter <- 0
    s <-1
    fx <- f(x)
    w <- x+s*u
    d_uf <- sum(u*grad) # directional derivative of f in direction of u at x
    decreasedEnough <- f(w) <= fx+alpha*s*d_uf
    while(!(condition(w) && decreasedEnough)  && iter<50){

        w <- x+s*u
        iter <- iter+1
        s <- s*bta
        decreasedEnough <- f(w) <= fx+alpha*s*d_uf
    }
    if(iter==50)
        stop("\nbacktrack: condition not satisfied after 50 iterations.\n")
    w
}


# The solvers terminate if the norm of a residual is below some tolerance
# eps. What the residual is, depends on the algorithm. These residuals are
# functions
#                    residual=residual(u,...),
#
# where u is the search direction computed at the current iterate x.
# Other data are the gradient grad(f)(x), deviation b-Ax, etc.
#
# Since we only need the norm of the residual, this norm will be implemented
# instead of the residual itself.
# The different functions will all have the same name but be distinguished
# by their parameter signature.



# Residual for the unconstrained Newton solver: |f_u(x)|=|(u,grad(f)(x))|.
# This is the directional derivative i.e. the rate of decrease of f in the
# search direction u.
#
# @return a function residualNorm(u)
#
residualNormUC <- function(u,grad) abs(sum((u*grad)))

# Residual for the equality constrained KKT solver (no inequality constraints),
# (see boyd, 10.3.1, p531):
#                                res(x) = (f_u(x),b-Ax)

# @param grad: gradient of f at x.
# @param r: residual Ax-b
# @return a boolean function terminate=terminate(u)
#
residualNormEC <- function(u,grad,r){

    f_ux <- sum(u*grad)
    sqrt(f_ux*f_ux+crossprod(r,r)^2)
}




# Solution of the KKT condition (grad(f)(x)=0) for minimization of f without
# constraints on a region C defined by condition, when it is known that a
# minimum exists at an interior point of this region.
# This is used for the Barrier method without equality constraints.
#
# The region C={x\in R^n : condition(x)} is then irrelevant for the necessary
# first order KKT conditions which simply become grad(f)(x)=0. But the region C
# is relevant for the solver which starting at a point x0 satisfying this
# condition must maintain it for all iterates.
#
# @param x0: starting point satisfying condition(x0)
# @param f: objective function
# @param grad:  function of x, grad(x)=gradient(f)(x)
# @param H:  function of x, H(x)=hessian(f)(x)
# @param condition: function of x, condition(x)=TRUE/FALSE
# @param bta: shrink factor used in backtracking.
# @param eps: termination as soon as |(u,grad(x))|<2eps, where u is the
# next search direction (i.e. when no more decreasse in f can be found in
# the Newton search direction).
#
solveNewton <- function(x0,f,grad,H,condition,alpha,bta,maxIter,eps){

     iter <- 0
     x <- x0
     resNorm <- 1+2*eps
     while(resNorm>eps && iter<maxIter){
     
         iter <- iter+1
         fx <- f(x)
         grad_x <-  grad(x)
         u <- solveNewtonStep(H(x),grad_x)
         resNorm <- residualNormUC(u,grad_x)
         # step to next iterate
         if(resNorm>eps)
            x <- backtrack(x,u,f,grad_x,condition,alpha,bta)
     }
     if(iter==maxIter)
         cat("\nsolveNewton: iteration limit reached at iteration ",iter,"\n")
     x
}



# Solution of the KKT condition
#                                 grad_xL(x,\nu)=0
#
# for minimization of f with equality constraints of the form Ax=b but no
# inequality constraints on a region C defined by condition, when it is known
# that a minimum exists at an interior point of this region.
#
# Here L=L(x,\nu)=f(x)+\nu'(Ax-b) is the Lagrangian function of the problem.
# This is used for the Barrier method with equality constraints.
#
# The region C={x\in R^n : condition(x)} is then irrelevant for the necessary
# first order KKT conditions but is relevant for the solver which starting at a
# point x0 satisfying this condition must maintain it for all iterates.
#
# This is an infeasible start algorithm which does not assumes Ax0=b.
#
# @param x0: starting point satisfying condition(x0) but need not satisfy
# Ax0=b.
# @param f: objective function
# @param grad:  function of x, grad(x)=gradient(f)(x)
# @param f: objective function.
# @param H:  function of x, H(x)=hessian(f)(x)
# @param condition: function of x, condition(x)=TRUE/FALSE
# @param A: matrix defining equality constraints Ax=b
# @param b: right hand side of equality constraints Ax=b
# @param bta: shrink factor used in backtracking.
# @param eps: termination as soon as |(u,grad(x))|<2eps, where u is the
# next search direction (i.e. when no more decrease in f can be found in
# the Newton search direction).
#
# @return w=c(x,nu), where nu is the Lagrange multiplier at the optimal
# point x.
#
solveKKT_NoIneq <- function(x0,f,grad,H,A,b,condition,alpha,bta,maxIter,eps){

     rho <- 1         # transform H -> H+rho*A'A to get positive definiteness
     n <- length(x0)
     p <- nrow(A)
     nu0 <- rep(0,p)  # starting Lagrange multiplier
     x <- x0
     w <- c(x0,nu0)   # starting point for solution
     # residual b-Ax, remains at zero at all iterates:
     r <- rep(0,length(x0))
     
     iter <- 0
     resNorm <- 1+2*eps
     while(resNorm>eps && iter<maxIter){

         iter <- iter+1
         x <- w[1:n]
         r_p <- b-A%*%x        # primal residual
         fx <- f(x)
         grad_x <- grad(x)
         u <- solveKKTStep_NoIneq(H(x),A,grad_x,r,rho)
         resNorm <- residualNormEC(u,grad_x,r_p)
         # step to next iterate
         if(resNorm>eps)
            w <- backtrack(w,u,f,grad_x,condition,alpha,bta)
     }
     if(iter==maxIter)
         cat("\nsolveKKT_NoIneq: iteration limit reached at iteration ",iter,"\n")
     w
}
