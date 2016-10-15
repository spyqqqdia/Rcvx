setClass("C",
representation(
 .ID   = "character",
 .dim  = "numeric",
 .date = "Date"
),
prototype(
 .ID   = NULL,
 .dim  = NaN,
 .Date = NULL
),
contains = "Base"
)
setGeneric("memFun",
function(this,pars) standardGeneric("memFun")
)
setMethod("memFum",
signature(this="C",pars="list"),
definition = function(this,pars){
            
    print(pars)
})
