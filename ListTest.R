# Warning example: R does not support reference semantics.
# If an an object is changed inside a function a local copy is made
# and operated on.
#
# You have to return that copy from the function and reassign that
# copy to the original object when using the function to change the
# original object.
#
# Note how we have to set up and use the member function add
# so the ListHolder object is changed.
#
setClass("ListHolder",
representation(
 .list   = "list"
),
prototype(
 .list   = list()
))
# The wrong way to add an element to the list.
# This function will not change this ListHolder.
#
setGeneric("addElementBadly",
function(this,i) standardGeneric("addElementBadly")
)
setMethod("addElementBadly",
signature(this="ListHolder",i="numeric"),
definition = function(this,i){
            
    k <- length(this@.list)
    this@.list[[k+1]] <- i
})
# The right way to add an element to the list, 
# return the modified copy.
# But note how it must be used in testListHolder
#
setGeneric("addElement",
function(this,i) standardGeneric("addElement")
)
setMethod("addElement",
signature(this="ListHolder",i="numeric"),
definition = function(this,i){
            
    k <- length(this@.list)
    this@.list[[k+1]] <- i
    this
})

testListHolder <- function(){
  
    cat("\nAdding elements 1,...,10 the wrong way:\n")
    lh <- new("ListHolder")  
    for(i in 1:10) add(lh,i)
    print(lh@.list)
    cat("\nAdding elements 1,...,10 the right way:\n")
    lh <- new("ListHolder")  
    for(i in 1:10) lh <- add(lh,i)
    print(lh@.list)
}

testListHolder()