# Build a list of inequality constraints, then check if they are 
# satisfied at a point x=x0+Fu.
#
# Then variable transform them via x=x0+Fu and check them at a point u0
# satisfying x=x0+Fu which must yield the same violated/satisfied pattern.
#
# The constraints will be constructed so that we know in advance 
# if they are violated or satisfied at x and then the constraint ID will
# be set accordingly.
# 
# Then iterate over all constraints and do the following:
# (a) report the status of the constraint (satisfied/active)
# (b) affine transform the constraint with the variable transformation
#     x = x0+Fu
#     Then compute the value at x and of the transformed constraint at u.
#     Verify that the values are equal and report the error.
#
# @param x point where they are to be checked.
#
basicConstraintSetTest <- function(x0,F,u){
  
    n <- nrow(F)
    m <- ncol(F)
    if(length(x0)!=n || length(u)!=m){
      
        msg <- "\nWe must have length(x0)=nrow(F) and length(u)=ncol(F).\n"
        msg <- paste(msg,"but nrow(F) =",n,"ncol(F) =",m," length(x0) =",length(x0))
        msg <- paste(msg,"and length(u) = ",length(u))
        stop(msg)
    }
    x <- x0+F%*%u
    x <- as.numeric(x)
  
    ID <- "InequalityTestSet"
    n <- length(x)           # dimension
    cat("\nBuilding a list of 10 linear and 10 quadratic inequality constraints.\n")
    cat("Alternately satisfied and active with tolerance 1e-7 and violated.\n")
    # 10 linear inequalities 
    
    listOfInequalities <-list()
    k <- 0
    for(ll in 1:10){
      
        p <- 2*runif(n)-1 
        if(ll%%2==0){  # add a constraint satisfied by x
          
            ID <- paste("Satisfied LinIneq",ll)  
            r <- -crossprod(p,x)-0.9e-7
            ineq <- newLinearInequality(ID,as.numeric(r),p)
          
        } else { # add a constraint not satisfied by x
          
            ID <- paste("Violated LinIneq",ll)  
            r <- -crossprod(p,x)+1e-7
            ineq <- newLinearInequality(ID,as.numeric(r),p)
        }
        listOfInequalities[[k+1]] <- ineq 
        k <- k+1
    }
    # 10 quadratic inequalities 
    for(qq in 1:10){
      
        p <- 2*runif(n)-1 
        A <- matrix(2*runif(n*n),nrow=n)
        Q <- crossprod(A)
        if(qq%%2==0){  # add a constraint satisfied by x
          
            ID <- paste("Satisfied QuadIneq",qq)  
            r <- -crossprod(p,x)-crossprod(x,Q%*%x)-0.9e-7
            ineq <- newQuadraticInequality(ID,as.numeric(r),p,Q)
          
        } else { # add a constraint not satisfied by x
          
            ID <- paste("Violated QuadIneq",qq)  
            r <- -crossprod(p,x)-crossprod(x,Q%*%x)+1e-7
            ineq <- newQuadraticInequality(ID,as.numeric(r),p,Q)
        }
        listOfInequalities[[k+1]] <- ineq 
        k <- k+1
    }
    # now check the constraints at the point x
    cat("\nReporting status of constraints and verifying that\n")
    cat("affine transform preserves value:\n\n")
    eps <- 1e-7
    for(ineq in listOfInequalities){
      
        reportStatusAt(ineq,x,eps)
        ineqTransformed <- affineTransformed(ineq,x0,F)
        val <- valueAt(ineq,x)
        val1 <- valueAt(ineqTransformed,u)
        cat("Error in affine transform:",val-val1,"\n")
    }
}