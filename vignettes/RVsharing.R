### R code from vignette source './RVsharing.Rnw'

###################################################
### code chunk number 1: options
###################################################
  options(width=75, continue = " ")
  library("Bureau")


###################################################
### code chunk number 2: RVsharing.toy2
###################################################
id <- paste0("sub", 1:4)
fa.id <- c(NA,NA,"sub1","sub1")
ma.id <- c(NA,NA,"sub2","sub2")
test.ped <- pedigree(id = id, dadid = fa.id, momid = ma.id, sex = c(1,2,1,2))


###################################################
### code chunk number 3: showped
###################################################
test.ped


###################################################
### code chunk number 4: plotped
###################################################
plot(test.ped)


###################################################
### code chunk number 5: kinship
###################################################
2*kinship(test.ped)


###################################################
### code chunk number 6: test
###################################################
RVsharing(test.ped)


###################################################
### code chunk number 7: RVsharing.toy
###################################################
id <- paste0("sub", 1:8)
fa.id <- c(NA,NA,NA,"sub1","sub1", NA,"sub3","sub5")
ma.id <- c(NA,NA,NA,"sub2","sub2", NA,"sub4","sub6")
test.ped <- pedigree(id = id, dadid = fa.id, momid = ma.id, sex = c(1,2,1,2,1,2,2,2))


###################################################
### code chunk number 8: showped
###################################################
test.ped


###################################################
### code chunk number 9: plotped2
###################################################
plot(test.ped)


###################################################
### code chunk number 10: kinship
###################################################
2*kinship(test.ped)


###################################################
### code chunk number 11: test
###################################################
RVsharing(test.ped)


