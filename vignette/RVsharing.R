### R code from vignette source './RVsharing.Rnw'

###################################################
### code chunk number 1: options
###################################################
  options(width=75, continue = " ")
  library("Bureau")


###################################################
### code chunk number 2: RVsharing.toy
###################################################
id <- paste0("sub", 1:8)
fa.id <- c(NA,NA,NA,"sub1","sub1", NA,"sub3","sub5")
ma.id <- c(NA,NA,NA,"sub2","sub2", NA,"sub4","sub6")
Bureau:::RVsharing( id = id, dad.id = fa.id, mom.id = ma.id )


