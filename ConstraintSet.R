cat("Reading ConstraintSet.R\n")
flush.console()


# Inequality and equality constraints enter the theory and algorithms differently and so
# no common base class will be defined.
#
# Generally the equality constraints will be of the affine form Ax=b so it would seem 
# superfluous to even introduce the notion of an EqualityConstraintSet: why not just 
# parametrize this set with the matrix A and vector b?
#
# The reason we do not follow this route is the following: in practical problems each
# constraint has a definite economic meaning, we therefore want to be able to refer to
# it by name and not just as a row in a matrix.
#
# Needless to say each EqualityConstraintSet will export the corresponding matrix A and
# vector b for use in the algorithms.
#

# Set (list) of constraints of the form g_i(x)<=0. 
# Such sets can compute the gradients and Hessians of barrier functions and the like.
#
setClass("SetOfInequalityConstraints",
representation(
   .ID   = "character",
   .dim  = "numeric",
   .constraints  = "list"
),
prototype(
   .ID   = NULL,
   .dim  = NULL,
   .constraints  = list()
))
# Creates empty set of inequality constraints.
# @param ID name of set
# @param d  dimension
#
newSetOfInequalityConstraints <- function(ID,d){
  
    ineqs <- new("SetOfInequalityConstraints") 
    ineqs@.ID <- ID
    ineqs@.dim <- d
    ineqs
}
# Add a single inequality constraint ineq to this set.
setGeneric("addConstraint",
function(this,constraint) standardGeneric("addConstraint")
)
setMethod("addConstraint",
signature(this="SetOfInequalityConstraints",constraint="Inequality"),
definition = function(this,constraint){

    if(constraint@.dim!=this@.dim){
    
        msg <- paste("SetOfInequalityConstraints ",this@.ID,":\n",sep="")
        msg <- paste(msg,"adding constraint of dimension",constraint@.dim)
        msg <- paste(msg,"while excpected dimension is",this@.dim)
        stop(msg)
    }
    k <- length(this@.constraints)
    this@.constraints[[k+1]] <- constraint
cat("\nnumber of constraints: ",length(this@.constraints))
})
# Add a set of inequality constraints ineqs to this set.
setGeneric("addConstraints",
function(this,constraints) standardGeneric("addConstraints")
)
setMethod("addConstraints",
signature(this="SetOfInequalityConstraints",constraints="SetOfInequalityConstraints"),
definition = function(this,constraints){
  
    if(constraints@.dim!=this@.dim){
    
        msg <- paste("SetOfInequalityConstraints ",this@.ID,":\n",sep="")
        msg <- paste(msg,"adding constraints of dimension",constraints@.dim)
        msg <- paste(msg,"while excpected dimension is",this@.dim)
        stop(msg)
    }
    k <- length(this@.constraints)
    for(constraint in constraints@.constraints){
        this@.constraints[[k+1]] <- constraint
        k <- k+1
    }
})

# Report the status (satisfied, active) of each constraint at the point x.
#
setGeneric("reportStatusOfConstraintsAt",
function(this,x) standardGeneric("reportStatusOfConstraintsAt")
)
setMethod("reportStatusOfConstraintsAt",
signature(this="SetOfInequalityConstraints",x="numeric"),
definition = function(this,x){
    for(ineq in this@.constraints) reportStatusAt(ineq,x)
})





# Set (list) of _affine_ equality constraints of the form a'x=r, where a,x are column 
# vectors and r is a number. Such sets convert the constraints to matrix form Ax=b. 
#
setClass("SetOfEqualityConstraints",
representation(
   .ID   = "character",
   .dim  = "numeric",
   .constraints  = "list",
   .A    = "matrix",
   .b    = "numeric"
),
prototype(
   .ID   = NULL,
   .dim  = NULL,
   .constraints  = list(),
   .A = NULL,
   .b = NULL
))
# @return Set of equality constraints initialized with the constraints Ax=b.
# We can then add additional equality constraints as needed.
#
newSetOfEqualityConstraints <- function(ID,A,b){
  
    m <- nrow(A)
    n <- ncol(A)
    lb <- length(b)
    if(lb!=m){
      
        msg <- paste("newSetOfEqualityConstraints",ID,": nrow(A) must equal length(b)")
        msg <- paste(msg, "\nbut nrow(A) =",m," and length(b) =",lb,".\n")
        stop(msg)
    }
    eqs <- new("SetOfEqualityConstraints")
    eqs@.dim <- n
    eqs@.A <- A
    eqs@.b <- b
    eqs
}
setMethod("addConstraint",
signature(this="SetOfEqualityConstraints",constraint="LinearEquality"),
definition = function(this,constraint){

    if(constraint@.dim!=this@.dim){
    
        msg <- paste("SetOfEqualityConstraints ",this@.ID,":\n",sep="")
        msg <- paste(msg,"adding constraint of dimension",constraint@.dim)
        msg <- paste(msg,"while expected dimension is",this@.dim)
        stop(msg)
    }
    k <- length(this@.constraints)
    this@.constraints[[k+1]] <- constraint
    this@.b <- c(this@.b,constraint@.r)
    this@.A <- rbind(this@.A,constraint@.a)
})
setMethod("addConstraints",
signature(this="SetOfEqualityConstraints",constraints="SetOfEqualityConstraints"),
definition = function(this,constraints){
  
    if(constraints@.dim!=this@.dim){
    
        msg <- paste("SetOfEqualityConstraints ",this@.ID,":\n",sep="")
        msg <- paste(msg,"adding constraints of dimension",constraints@.dim)
        msg <- paste(msg,"while expected dimension is",this@.dim)
        stop(msg)
    }
    k <- length(this@.constraints)
    for(constraint in constraints@.constraints){
        this@.constraints[[k+1]] <- constraint
        k <- k+1
    }
    this@.b <- c(this@.b,constraints@.b)
    this@.A <- rbind(this@.A,constraints@.A)
})




    


