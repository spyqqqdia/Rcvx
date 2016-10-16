# Class to solve convex problems with affine equality and convex 
# inequality constraints.
# Base class to Barrier and Primal-Dual interior point methods.
#
# The optimization will not in general be carried out in the original variable x.
# We might employ a change of variables. For example to eliminate affine equality
# constraints Ax=b we employ a change of variables x = x0 +Fu, where Im(F)=ker(A).
#
# The letter u is used to denote the variable in which the optimzation is carried 
# out.
#
setClass("ConvexOptimizer",
representation(
   .ID             = "character",
   .dim            = "numeric",
   .f              = "function",           # objective function f=f(u)
   .gradf          = "function",           # gradient of f(u)
   .Hf             = "function",           # hessian of f(u)
   .inequalities   = "SetOfInequalityConstraints",
   .equalities     = "SetOfEqualityConstraints"
),
prototype(
   .ID             = "NullOptimizer",
   .dim            = NaN,
   .f              = NULL,
   .gradf          = NULL,
   .Hf             = NULL,
   .inequalities   = NULL,
   .equalities     = NULL
))



# Minimize starting from initial point x=x0 with parameters
# alpha, bta for the backtracking line search.
#
# @param eps: termination as soon as norm of residual beow eps.
#
setGeneric("optimize",
function(this,x0,alpha,bta,eps,verboseLevel) standardGeneric("optimize")
)


# Minimize without starting point with parameters
# alpha, bta for the backtracking line search.
#
setGeneric("optimize",
function(this,alpha,bta,eps,verboseLevel) standardGeneric("optimize")
)


# Feasibility analysis (Phase I in boyd terminology).
# @return feasible starting point or certificate of infeasibility.
#
setGeneric("findStartingPoint",
function(this) standardGeneric("findStartingPoint")
)



# Abstract base class will be fully intialized from concrete 
# subclasses, therefore no constructor factory method defined here.


#-- Declaring some abstract functions which can be implemented for some
#-- concrete subclasses but do not have to be implemented for all of these.

# The objective function $f_t$ with weight t in the barrier method:
#
#            $f_t(x) = f(x)-t\sum_i\log(-g_i(x))$,
#
# where $g_i(x)\leq 0$ are the inequality constraints. In this formulation
# t -> 0 in the iterations of the optimization.
#
setGeneric("barrierFunction",
function(this,u,t) standardGeneric("barrierFunction")
)


# Gradient of barrier function.
#
setGeneric("gradientBarrierFunction",
function(this,u,t) standardGeneric("gradientBarrierFunction")
)


# Hessian (2nd derivative) of barrier function.
#
setGeneric("hessianBarrierFunction",
function(this,u,t) standardGeneric("hessianBarrierFunction")
)
