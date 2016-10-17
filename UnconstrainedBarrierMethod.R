# Class to solve convex problems by eliminating not only the 
# inequality constraints (as usual in the barrier method) but also 
# the equality constraints so that the problem is reduced to an
# unconstrained convex optimization.
#
# The solutions to the equality constraints Ax=b are represented as
# x = x0+Fu, where x0 is a special solution to Ax=b and Im(F)=ker(A),
# see docs/nullspace.pdf.
#
# The problem is then reduced to a sequence of unconstrained convex 
# problems in the variable u.
#
setClass("UnonstrainedBarrierMethod",
representation(
   .x0 = "numeric",
   .F  = "matrix",
   .dimReducedInequalities = "list"
),
prototype(
   .x0 = NULL,
   .F  = NULL,
   .dimReducedInequalities = list()    # inequalities after change of variables x = x0+Fu
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
newUnconstrainedBarrierMethod <- function(ID,f,gradf,Hf,inequalities,equalities){
  
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
    ucbm@.f <- f
    ucbm@.gradf <- gradf
    ucbm@.Hf <- Hf
    ucbm@.inequalities <- inequalities
    ucbm@.equalities <- equalities
    
    A <- equalities@.A
    b <- equalities@.b
    
    S <- solutionSpace(A,b)
    this@.x0 <- S$.x0
    this@.F <- S$.Q
    
    # set objective function, gradient and hessian as functions of the dimension reduced
    # variable u in the change of variables x = x0+Fu
    this@.f <- function(u) f(x0+F%*%u)
    this@.gradf <- function(u) crossprod(F,gradf(x0+F%*%u))
    this@.Hf <- function(u) crossprod(F,Hf(x0+F%*%u)%*%F)
    
    # inequalities after variable transform x -> u, x = x0+Fu
    ineqs <- inequalities@.constraints
    i <- 0
    for(ineq in ineqs){
      
        this@.dimReducedInequalities[[i]] <- affineTransformed(ineq,this@.x0,this@.F)
        i <- i+1
    }
    ucbm
}

# See ConvexOptimizer.
#
setMethod("barrierFunction",
signature(this="UnconstrainedBarrierMethod",u="numeric",t="numeric"),
definition = function(this,u,t){
            
    sum <- 0.0
    for(ineq in this@.dimReducedInequalities)
        sum <- sum - log(-valueAt(ineq,u))
    t*this@.f(u) + sum
})

# See ConvexOptimizer.
#
setMethod("gradientBarrierFunction",
signature(this="UnconstrainedBarrierMethod",u="numeric",t="numeric"),
definition = function(this,u,t){
            
    sum <- 0.0
    for(ineq in this@.dimReducedInequalities)
        sum <- sum - gradientAt(ineq,u)/valueAt(ineq,u)
    t*this@.gradf(u) + sum
})

# See ConvexOptimizer.
#
setMethod("hessianBarrierFunction",
signature(this="UnconstrainedBarrierMethod",u="numeric",t="numeric"),
definition = function(this,u,t){
            
    sum <- 0.0
    for(ineq in this@.dimReducedInequalities){
      
        r <- valueAt(ineq,u)
        g <- gradientAt(ineq,u)
        sum <- sum - hessianAt(ineq,u)/r + outer(g,g)/(r*r)
    }
    t*this@.Hf(u) + sum
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
#
setMethod("optimize",
signature(
    this="UnconstrainedBarrierMethod",x0="numeric",alpha="numeric",bta="numeric",
    eps="numeric", verboseLevel="numeric"
),
definition = function(this,x0,alpha,bta,eps,verboseLevel){
            
    if(!isStrictlyFeasible(this,x0)){
    
        cat("\nUnconstrainedBarrierMethod::optimize: starting point not strictly feasible.") 
        cat("\nStatus of constraints at starting point:\n")
        reportStatusOfConstraintsAt(this@.inequalities,x0)
        stop("Terminating.\n")
    }
    # recall duality gap at central point x(t) is m/t,  where m is the number
    # of inequalitiy constraints, [boyd], p565,566
    
    
  
  
})





