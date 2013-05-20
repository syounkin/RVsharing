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
test.ped <- pedigree(id = id, dadid = fa.id, momid = ma.id, sex = c(1,2,1,2,1,2,2,2))


###################################################
### code chunk number 3: showped
###################################################
test.ped


###################################################
### code chunk number 4: plotped
###################################################
plot(test.ped)


###################################################
### code chunk number 5: test
###################################################
RVsharing(test.ped)


