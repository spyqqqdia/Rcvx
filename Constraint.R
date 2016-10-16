cat("Reading Constraint.R\n")
flush.console()

# General constraint of an optimization defined in terms of the value 
# of a function g at a vector x (typically g(x)<=0 but can also be g(x)=0
# if g is of the affine form g(x)=Ax-b).
#
setClass("Constraint",
representation(.ID="character", .dim="numeric"),
prototype(.ID=NULL, .dim=NULL)
)
setGeneric("identity",
function(this) standardGeneric("identity")
)
# Checks if length(x)==this@.dim.
setGeneric("checkDimension",
function(this,x) standardGeneric("checkDimension")
)
setMethod("checkDimension",
signature(this="Constraint",x="numeric"),
definition = function(this,x) if(length(x)!=this@.dim){
  
    msg <- paste(this@.ID,": dimension mismatch,",sep="")
    msg <-paste(msg,"dim =",dim,", dim(x) =",length(x))
    stop(msg) 
})
# value f(x) of function defining the constraint
setGeneric("valueAt",
function(this,x) standardGeneric("valueAt")
)
setGeneric("gradientAt",
function(this,x) standardGeneric("gradientAt")
)
setGeneric("hessianAt",
function(this,x) standardGeneric("hessianAt")
)
setGeneric("isSatisfiedAt",
function(this,x) standardGeneric("isSatisfiedAt")
)
# @return the constraint after a variable transform x=x0+Fu,
# where F is a matrix. Returns the equivalent constraint on the variable u.
#
setGeneric("affineTransformed",
function(this,x0,F) standardGeneric("affineTransformed")
)





# Linear equality constraint of the form h(x)=a'x=r
# where p,x are column vectors and r a number.
#
setClass("LinearEquality",
representation(.r="numeric", .a="numeric"),
prototype(.r=NULL, .p=NULL),
contains = "Constraint"
)
# factory function (constructor).
#@return an equality constraint of the form a'x=r, 
# where a,x are column vectors and r a number.
#
newLinearEquality <- function(ID,a,r){
  
    eq <- new("LinearEquality")  
    eq@.ID <- ID
    eq@.dim <- length(a)
    eq@.a <- a
    eq@.r <- r
    eq
}
#@return left hand side p'x
setMethod("valueAt",
signature(this="LinearEquality",x="numeric"),
definition = function(this,x){
    checkDimension(this,x)
    sum(this@.a*x)  
})
#@return abs(p'x-r)<1e-14
setMethod("isSatisfiedAt",
signature(this="LinearEquality",x="numeric"),
definition = function(this,x) abs(valueAt(x)-this@.r)<1e-14
)
setMethod("gradientAt",
signature(this="LinearEquality",x="numeric"),
definition = function(this,x){
    checkDimension(this,x)  
    this@.a
})
setMethod("hessianAt",
signature(this="LinearEquality",x="numeric"),
definition = function(this,x){ 
    checkDimension(this,x)
    matrix(0,nrow=this@.dim,ncol=this@.dim)
})
# The method affineTransformed will not be defined for linear
# equalities since the purpose of the transform is to eliminate
# the equality constraints.





# Inequality constraint g(x)<=0
#
setClass("Inequality",
contains = "Constraint"
)
#@return g(x)<=0
setMethod("isSatisfiedAt",
signature(this="Inequality",x="numeric"),
definition = function(this,x) valueAt(this,x)<=0
)
#@return g(x)<-eps
setGeneric("isSatisfiedStrictlyAt",
function(this,x,eps) standardGeneric("isSatisfiedStrictlyAt")
)
setMethod("isSatisfiedStrictlyAt",
signature(this="Inequality",x="numeric",eps="numeric"),
definition = function(this,x,eps) valueAt(this,x)<(-eps)
)
#@return TRUE/FALSE according as the constraint is satisfied with equality
# at x ( |g(x)|<eps ).
setGeneric("isActiveAt",
function(this,x,eps) standardGeneric("isActiveAt")
)
setMethod("isActiveAt",
signature(this="Inequality",x="numeric",eps="numeric"),
definition = function(this,x,eps) abs(valueAt(this,x))<eps
)
# Print message wether the constraint is satisfied or not
# and wether it is active.
setGeneric("reportStatusAt",
function(this,x,eps) standardGeneric("reportStatusAt")
)
setMethod("reportStatusAt",
signature(this="Inequality",x="numeric",eps="numeric"),
definition = function(this,x,eps){
  
    msg <- paste("\nInequality",this@.ID)  
    if(!isSatisfiedAt(this,x)) 
        msg <- paste(msg,"is not satisfied at x.\n")
    else if(isActiveAt(this,x,eps))
        msg <- paste(msg,"is satisfied and active at x with tolerance",eps,"\n")
    else
       msg <- paste(msg,"is satisfied but not active at x with tolerance",eps,"\n")
    
    cat(msg)
})






# Linear inequality constraint of the form g(x)=r+p'x<=0
# with Q a symetric matrix.
# This is convex if and only if Q is positive semidefinite.
#
setClass("LinearInequality",
representation(.r="numeric", .p="numeric"),
prototype(.r=NULL, .p=NULL),
contains = "Inequality"
)
# factory method (constructor)
#@return a LinearInequality of the form r+p'x<=0
newLinearInequality <- function(ID,r,p){
  
    ineq <- new("LinearInequality")
    ineq@.ID <- ID; ineq@.dim <- length(p); 
    ineq@.r <- r;   ineq@.p <- p
    ineq
}
setMethod("valueAt",
signature(this="LinearInequality",x="numeric"),
definition = function(this,x){
    checkDimension(this,x)
    this@.r+crossprod(this@.p,x)  
})
setMethod("gradientAt",
signature(this="LinearInequality",x="numeric"),
definition = function(this,x){
    checkDimension(this,x)  
    this@.p
})
setMethod("hessianAt",
signature(this="LinearInequality",x="numeric"),
definition = function(this,x){ 
    checkDimension(this,x)
    matrix(0,nrow=this@.dim,ncol=this@.dim)
})
# @return the constraint after a variable transform x=x0+Fu,
# where F is a matrix. Returns the equivalent constraint on 
# the variable u.
#
# The inequality becomes r+p'(x0+Fu)<=0, i.e.: 
# r -> r+p'x0, p -> F'p.
#
setMethod("affineTransformed",
signature(this="LinearInequality",x0="numeric",F="matrix"),
definition = function(this,x0,F){
    checkDimension(this,x0)
    d <- this@.dim
    if(nrow(F)!=d){
      
        msg <-paste("\nConstraint",this@.ID," affineTransformed: nrow(F) must be",d)
        msg <-paste(msg,"\nbut is =",nrow(F),"\n")
        stop(msg)
    }
    ID <- paste(this@.ID,"affine transformed")
    r <- this@.r
    p <- this@.p
    s <- r+crossprod(p,x0)
    q <- t(F)%*%p
    newLinearInequality(ID,as.numeric(s),as.numeric(q))
})




# Quadratic inequality constraint of the form g(x)=r+p'x+(1/2)x'Qx<=0
# with Q a _symetric_ matrix.
# This is convex if and only if Q is positive semidefinite.
#
setClass("QuadraticInequality",
representation(.r="numeric", .p="numeric", .Q="matrix"),
prototype(.r=NULL, .p=NULL, .Q=NULL),
contains = "Inequality"
)
# factory method (constructor)
#@return a LinearInequalityConstraint of the form r+p'x<=0
newQuadraticInequality <- function(ID,r,p,Q){
  
    if(class(Q)!="matrix"){
    
        msg <- "\nnewQuadraticInequality: parameter Q must be a matrix,"
        msg <- paste(msg,"\nbut class(Q) =",class(Q))
        stop(msg)
    } 
    if(length(p)!=nrow(Q) || length(p)!=nrow(Q)){
    
        msg <- "\nnewQuadraticInequality: we must have length(p)=nrow(Q)=ncol(Q)"
        msg <- paste(msg,"\nbut length(p) = ",length(p),", nrow(Q) = ",nrow(Q),", ncol(Q) = ",ncol(Q))
        stop(msg)
    
    }
    # check symmetry
    symmErr <- max(abs(Q-t(Q)))
    if(symmErr>1e-12){
      
        msg <- "newQuadraticInequality: matrix Q not symmetric:\n"
        msg <- paste(msg,"max(abs(Q-t(Q))) =",symmErr,"\n")
        stop(symmErr)
    }
    ineq <- new("QuadraticInequality")
    ineq@.ID <- ID; ineq@.dim <- length(p); ineq@.r <- r; 
    ineq@.p <- p; ineq@.Q <- Q
    ineq
}
setMethod("valueAt",
signature(this="QuadraticInequality",x="numeric"),
definition = function(this,x){
    checkDimension(this,x)
    this@.r+sum(this@.p*x)+sum(x*(this@.Q%*%x))/2
})
setMethod("gradientAt",
signature(this="QuadraticInequality",x="numeric"),
definition = function(this,x){
    checkDimension(this,x)
    p+this@.Q%*%x
})
setMethod("hessianAt",
signature(this="QuadraticInequality",x="numeric"),
definition = function(this,x){
    checkDimension(this,x)
    this@.Q
})
# @return the constraint after a variable transform x=x0+Fu,
# where F is a matrix. Returns the equivalent constraint on 
# the variable u.
#
# The inequality becomes r+p'(x0+Fu)+0.5*(x0+Fu)'Q(x0+Fu)<=0, equivalently
# r+(p+x0)'x0+[F'p+F'Qx0]'u+0.5*u'F'QFu <= 0.
# Thus the induced transformation on the data r,p,Q is the following
# r -> r+p'x0+0.5*x0'Qx0
# p -> F'(p+Qx0)   and
# Q -> F'QF
#
setMethod("affineTransformed",
signature(this="QuadraticInequality",x0="numeric",F="matrix"),
definition = function(this,x0,F){
    checkDimension(this,x0)
    d <- this@.dim
    if(nrow(F)!=d){
      
        msg <-paste("\nConstraint",this@.ID,": affineTransformed nrow(F) must be",d)
        msg <-paste(msg,"\nbut is =",nrow(F),"\n")
        stop(msg)
    }
    r <- this@.r
    p <- this@.p
    Q <-this@.Q
    
    ID <- paste(this@.ID,"affine transformed")
    s <- r+crossprod(x0,p+0.5*Q%*%x0)
    q <- crossprod(F,p+Q%*%x0)                     # F'(p+Qx0)
    R <- crossprod(F,Q%*%F)                        # F'QF
    newQuadraticInequality(ID,as.numeric(s),as.numeric(q),R)
})


