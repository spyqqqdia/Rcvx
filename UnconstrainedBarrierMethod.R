# Class to solve convex problems by eliminating not only the
# inequality constraints (as usual in the barrier method) but also
# the equality constraints so that the problem is reduced to an
# unconstrained convex optimization.
#
# The solutions to the equality constraints Ax=b are represented as
# x = z+Fu, where z is a special solution to Ax=b and Im(F)=ker(A),
# see docs/nullspace.pdf.
#
# The problem is then reduced to a sequence of unconstrained convex
# problems in the variable u.
#
setClass("UnonstrainedBarrierMethod",
representation(
  .x0 = "numeric",                  # starting point (full dimension)
  .w0 = "numeric",                  # change of variables x = w0+Fu, w0 special solution to Ax=b
  .F  = "matrix",                   # dimension reduction matrix
  .dimReducedInequalities = "list", # inequalities after change of variables x = w0+Fu
  .tMultiplier = "numeric"          # multiplies parameter t in barrier function
),
prototype(
  .x0 = NULL,
  .w0 = NULL,
  .F  = NULL,
  .dimReducedInequalities = list(),
  .tMultiplier = 5.0
),
contains = "ConvexOptimizer"
)

# Allocates an optimizer from list of equality and inequality
# constraints.
#
# @param f  f=f(x), objective function of original variable x (not dimension reduced u).
# @param gradf gradient of f, function of x.
# @param Hf hessian of f, function of x.
# @param inequalities: SetOfInequalityConstraints.
# @param equalities: SetOfEqualityConstraints.
#
newUnconstrainedBarrierMethod <- function(ID,x0,f,gradf,Hf,inequalities,equalities){

    d <- inequalities@.dim
    d1 <- equalities@.dim
    if(d!=d1){

        msg <- "\nnewUnconstrainedBarrierMethod: dim(inequalities) != dim(equalities):\n"
        msg <- paste(msg,"dim(inequalities) =",d," and dim(equalities) =",d1)
        stop(msg)
    }
    ucbm <- new("newUnconstrainedBarrierMethod")
    ucbm@.ID <- ID
    ucbm@.dim <- d
    ucbm@.x0 <- x0
    ucbm@.f <- f
    ucbm@.gradf <- gradf
    ucbm@.Hf <- Hf
    ucbm@.inequalities <- inequalities
    ucbm@.equalities <- equalities

    A <- equalities@.A
    b <- equalities@.b

    S <- solutionSpace(A,b)
    this@.w0 <- S$.x0          # special solution to Ax=b
    w0 <- this@.w0
    this@.F <- S$.Q

    # set objective function, gradient and hessian as functions of the dimension
    # reduced variable u in the change of variables x = x0+Fu
    this@.f <- function(u) f(w0+F%*%u)
    this@.gradf <- function(u) crossprod(F,gradf(w0+F%*%u))
    this@.Hf <- function(u) crossprod(F,Hf(w0+F%*%u)%*%F)

    # inequalities after variable transform x -> u, x = w0+Fu
    ineqs <- inequalities@.constraints
    i <- 0
    for(ineq in ineqs){

        this@.dimReducedInequalities[[i]] <- affineTransformed(ineq,this@.w0,this@.F)
        i <- i+1
    }
    ucbm
}

# See ConvexOptimizer.
#
setMethod("barrierFunction",
signature(this="UnconstrainedBarrierMethod",u="numeric",t="numeric"),
definition = function(this,u,tt){

    sumIneq <- 0.0
    for(ineq in this@.dimReducedInequalities)
        sumIneq <- sumIneq - log(-valueAt(ineq,u))
    tt*this@.f(u) + sumIneq
})

# See ConvexOptimizer.
#
setMethod("gradientBarrierFunction",
signature(this="UnconstrainedBarrierMethod",u="numeric",t="numeric"),
definition = function(this,u,tt){

    sumGradIneq <- 0.0
    for(ineq in this@.dimReducedInequalities)
        sumGradIneq <- sumGradIneq - gradientAt(ineq,u)/valueAt(ineq,u)
    tt*this@.gradf(u) + sumGradIneq
})

# See ConvexOptimizer.
#
setMethod("hessianBarrierFunction",
signature(this="UnconstrainedBarrierMethod",u="numeric",t="numeric"),
definition = function(this,u,tt){

    sumHessIneq <- 0.0
    for(ineq in this@.dimReducedInequalities){

        r <- valueAt(ineq,u)
        g <- gradientAt(ineq,u)
        sumHessIneq <- sumHessIneq - hessianAt(ineq,u)/r + outer(g,g)/(r*r)
    }
    tt*this@.Hf(u) + sumHessIneq
})

# @return TRUE if x0 satisfies the inequality constraints strictly,
# FALSE otherwise.
#
setGeneric("isStrictlyFeasible",
function(this,x0) standardGeneric("isStrictlyFeasible")
)
setMethod("isStrictlyFeasible",
signature(this="UnconstrainedBarrierMethod",x0="numeric"),
definition = function(this,x0){

    # find the corresponding variable u0, see Utils.R
    u0 <- solutionSpaceParameter(x0,this@.w0,this@.F)
    res <- TRUE
    for(ineq in this@.dimReducedInequalities)
        if(valueAt(ineq,u0)>=0) res <- FALSE
    res
})


# Minimization starting from x=x0. See ConvexOptimizer.
# @param eps keep going until the duality gap at current central point
# is below eps.
# @param alpha,beta: line search parameters.
#
setMethod("optimize",
signature(
    this="UnconstrainedBarrierMethod",x0="numeric",alpha="numeric",bta="numeric",
    eps="numeric", verboseLevel="numeric"
),
definition = function(this,x0,alpha,bta,eps,verboseLevel){

    maxIter <- 500    #---FIX ME: make this a parameter in the generic function
    if(!isStrictlyFeasible(this,x0)){

        cat("\nUnconstrainedBarrierMethod::optimize: starting point not strictly feasible.")
        cat("\nStatus of constraints at starting point:\n")
        reportStatusOfConstraintsAt(this@.inequalities,x0)
        stop("Terminating.\n")
    }
    # recall duality gap at central point x(t) is m/t,  where m is the number
    # of inequalitiy constraints, [boyd], p565,566
    # recall also that the optimization is carried out in the dimension
    # reduced variable u, where x = w0+Fu
    m <- length(this@.dimReducedInequalities)
    u <- solutionSpaceParameter(x0,this@.w0,this@.F)        # x0 -> u0
    condition <- function(v) isStrictlyFeasible(this,this@w0+F%*%v)
    tt <- 1
    gg <- this@.tMultiplier
    dualGap <- eps+1
    iter <- 1
    while((dualGap>eps)&&(iter<maxIter)){
   
        f <- function(u) barrierFunction(this,u,tt)
        grad <- function(u) gradientBarrierFunction(this,u,tt)
        H <- function(u) hessianBarrierFunction(this,u,tt)
        u <- solveNewton(u,f,grad,H,condition,alpha,bta,maxIter,eps)
        dualGap <- m/tt
        tt <- tt*gg
    }
    this@w0+F%*%u
})

