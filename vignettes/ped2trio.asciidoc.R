### R code from vignette source 'ped2trio.asciidoc'

###################################################
### code chunk number 1: junk
###################################################
par(mar=rep(0,4))


###################################################
### code chunk number 2: options
###################################################
  options(width=75, continue = " ")
  library("RVsharing")


###################################################
### code chunk number 3: RVsharing.toy2
###################################################
id <- paste0("sub", 1:8)
fa.id <- c(NA,NA,NA,"sub1","sub1", NA,"sub3","sub5")
ma.id <- c(NA,NA,NA,"sub2","sub2", NA,"sub4","sub6")
test.ped.first.cousins <- pedigree(id = id, dadid = fa.id, momid = ma.id, sex = c(1,2,1,2,1,2,2,2))


###################################################
### code chunk number 4: plotped2
###################################################
plot(test.ped.first.cousins)


###################################################
### code chunk number 5: firstcousins
###################################################

trio.obj <- RVsharing:::ped2trio(test.ped.first.cousins)

geno.vec <- c(1,0,0,NA,NA,0,NA,NA)
names(geno.vec) <- paste0("sub",1:8)

RVsharing:::GeneDropSim.fn( trio = trio.obj, geno.vec = geno.vec, dt.vec = c("sub7","sub8"), fd.indices = c(1,2,3,6), n = 1e4)
RVsharing(test.ped.first.cousins)


###################################################
### code chunk number 6: secondcousins
###################################################
id <- paste0("sub", 1:12)
fa.id <- c(NA,NA,NA,"sub1","sub1", NA,"sub3","sub5",NA,NA,"sub9","sub10")
ma.id <- c(NA,NA,NA,"sub2","sub2", NA,"sub4","sub6",NA,NA,"sub7","sub8")
test.ped.second.cousins <- pedigree(id = id, dadid = fa.id, momid = ma.id, sex = c(1,2,1,2,1,2,2,2,1,1,2,2))
plot(test.ped.second.cousins)

trio.obj <- RVsharing:::ped2trio(test.ped.second.cousins)

geno.vec <- c(1,0,0,NA,NA,0,NA,NA,0,0,NA,NA)
names(geno.vec) <- paste0("sub",1:12)

RVsharing:::GeneDropSim.fn( trio = trio.obj, geno.vec = geno.vec, dt.vec = c("sub11","sub12")
RVsharing(test.ped.second.cousins)



