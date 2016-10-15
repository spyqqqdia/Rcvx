cat("Reading TestFunctions.R")
flush.console()



# Even dimensional, decoupled Rosenbrook function:
#
# f(x)=sum_{i=1}^{n/2}[100(x_{2i}-x^2_{2i-1})^2+(1-x_{2i-1})^2]
#
# This function has an obvious valley along x_{2i}=x_{2i-1}^2  since the value
# at any point can be decreased by redefining the even numbered coordinates as
# x_{2i}=x_{2i-1}^2.
#
# In two dimensions the valley runs along y=x^2 and hits the minimum at
# x=1 (thus y=1).
#
# The valley is narrow because of the factor 100 (it gets narrower with
# steeper walls if this factor is increased). It is thus a good test function
# for Newton methods to see if the path to the minimum bounces around from one
# wall to the opposite wall necessitating many iterations before we converge
# at the solution.
#
# For the computation of the gradient and Hessian note that this is a sum of
# coordinate-disjoint functions f(x,y)=100*(y-x^2)^2+(x-1)^2
# with x=x_{2k-1} and y=x_{2k}.
#
# @returns list(.f,.grad,.H), where f=f(.) is the function itself,  grad its
# gradient and H its hessian (all are _functions_)
#
rosenbrook <- function(n){

    if(n%%2!=0)
        stop("\nrosenbrook: dimension must n be even but is n = ",n)
       
    m <- n/2
    evenIdx <- 2*(1:m)
    oddIdx  <- even-1
       
    # the function f itself
    f <- function(x){
   
       if(length(x)!=n){
       
          msg <- "\nrosenbrook(n): argument x must have length n ="
          msg <- paste(msg,n,"but length(x) = ",length(x))
          stop(msg)
       }

       x_even <- x[evenIdx]
       x_odd <- x[oddIdx]
       sum(100*(x_even-x_odd)^2+(x_odd-1)^2)
    }
    # the gradient of f
    grad <- function(x){
   
        if(length(x)!=n){

          msg <- "\nrosenbrook(n): argument x must have length n ="
          msg <- paste(msg,n,"but length(x) = ",length(x))
          stop(msg)
        }
        res <- rep(0,n)
        res[evenIdx] <- 200*(x[evenIdx]-x[oddIdx]^2)
        res[oddIx] <- -2*x[odd]*res[evenIdx]+2*(x[odd]-1)
       
        res
    }
    # the Hessian of f
    H <- function(x){

        if(length(x)!=n){

          msg <- "\nrosenbrook(n): argument x must have length n ="
          msg <- paste(msg,n,"but length(x) = ",length(x))
          stop(msg)
        }
        res <- matrix(0,nrow=n,ncol=n)
        i <- 1       # iterate through the odd indices < n
        while(i<n){
       
           # i=2k-1, j=2k, x:=x_{2k-1}=x[i], y:= x_{2k}=x[j]
           j <- i+1
           res[i,i] <- -400*x[j]+1200*x[i]*x[i]+2         # d^2f/dx^2
           res[i,j] <- -400*x[i]                          # d^2f/dxdy
           res[j,i] <- -400*x[i]
           res[j,j] <- 200                                # d^2f/dy^2
           i <- i+2                                       # next odd index
        }
        res
    }
    list(.f=f,.grad=grad,.H=H)
}
 
    

