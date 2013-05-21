setMethod("RVsharing",  signature(object="pedigree"), function(object){
    RVsharing.fn( id = object$id, dad.id = ifelse(object$findex!=0,object$id[object$findex],0), mom.id = ifelse(object$mindex!=0,object$id[object$mindex],0) )
})
