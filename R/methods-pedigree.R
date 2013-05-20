## setMethod("initialize", signature(.Object="PedClass"),
##          function(.Object, ...){
##            .Object
##            callNextMethod(.Object, ... )
##          })

## setMethod( "PedClass", signature(object = "DataFrame"), function(object) {
##   ped.DF <- with( as(object, "data.frame"), {
##     DataFrame( famid = factor(famid),
##               id = factor(id),
##               fid = factor(fid),
##               mid = factor(mid),
##               sex = factor(sex),
##               dx = factor(dx) )
##   })
##   new("PedClass", ped.DF)
## })

## setMethod("trios",  signature(object="PedClass"), function(object) {
##   trio.df <- subset( as(object, "data.frame" ), !is.na(fid) & !is.na(mid) & !is.na(id) )
## #  id.char <- as.character(trio.df$id)#,
##   trio.df.2 <- data.frame( id = as.character(trio.df$id),
##                           fid = as.character(trio.df$fid),
##                           mid = as.character(trio.df$mid), stringsAsFactors=FALSE )
## #                          )
##   return(trio.df.2)
## })

setMethod("RVsharing",  signature(object="pedigree"), function(object){
    RVsharing.fn( id = object$id, dad.id = object$id[ifelse(object$findex!=0,object$findex,NA)], mom.id = object$id[ifelse(object$mindex!=0,object$mindex,NA)] )
})

## setMethod("offspring",  signature(object="PedClass"), function(object){
##     with(as(object,"data.frame"), unique(as.character(id)))
## })

## setMethod("allSubjects",  signature(object="PedClass"), function(object){
##     with(as(object,"data.frame"), unique(c(as.character(id), as.character(fid), as.character(mid))))
## })
