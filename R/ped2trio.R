ped2trio <- function(ped)
{
# Converts a pedigree object into a recursive list of trio objects
# ped: a pedigree object

dv = kindepth(ped)
# md is the depth of the pedigree, i.e. the number of levels of trios
md = max(dv)

# Extracting the ids of subjects and their parents from ped object
id = as.character(ped$id)
dad.id = mom.id = rep("0",length(id))
dad.id[ped$findex>0] = ped$id[ped$findex]
mom.id[ped$mindex>0] = ped$id[ped$mindex]

trio.flag = logical(length(id))
names(trio.flag) = id

trio.list=list()
cumultrio.vec=numeric(md)

r=1
for (lev in md:2)
  {
  # Create vector of unique moms and dads at the current depth
  currentdads = unique(dad.id[dv==lev])
  currentmoms = unique(mom.id[dv==lev])
  # Create vector of unique dads at the current depth who are non founders
  dad.tokeep = dad.id[id%in%currentdads]!=0
  names(dad.tokeep) = id[id%in%currentdads]
  dad.trio = names(dad.tokeep)[dad.tokeep]
  # Loop over the non-founder dads at the current depth
  if (length(dad.trio)>0)
  {
  for (i in 1:length(dad.trio))
    {
    spousevec = unique(mom.id[dad.id==dad.trio[i]])
    for (j in 1:length(spousevec))
      {
      # Create vector of children of the couple
      offspring.vec = id[dad.id==dad.trio[i] & mom.id==spousevec[j]]
      # Creating the elements of the trio object in trio.list
      trio.obj=list(id=dad.trio[i],spouse=spousevec[j],offspring = as.list(ifelse(trio.flag[offspring.vec],paste("trio",offspring.vec,sep=""),offspring.vec)))
      trio.list[[r]]=trio.obj
      names(trio.list)[r]=paste("trio",dad.trio[i],sep="")
      r = r +1
      }
    }
  # We turn the trio.flag of the dad.trio to TRUE to indicate they now have a trio
  trio.flag[dad.trio] = TRUE
  }
  # Create vector of unique moms at the current depth who are non founders
  mom.tokeep = mom.id[id%in%currentmoms]!=0
  names(mom.tokeep) = id[id%in%currentmoms]
  mom.trio = names(mom.tokeep)[mom.tokeep]
  # Loop over the non-founder moms at the current depth
  if (length(mom.trio)>0)
  {
  for (i in 1:length(mom.trio))
    {
    spousevec = unique(dad.id[mom.id==mom.trio[i]])
    for (j in 1:length(spousevec))
      {
      # Create vector of children of the couple
      offspring.vec = id[mom.id==mom.trio[i] & dad.id==spousevec[j]]
      # Creating the elements of the trio object in trio.list
      trio.obj=list(id=mom.trio[i],spouse=spousevec[j],offspring = as.list(ifelse(trio.flag[offspring.vec],paste("trio",offspring.vec,sep=""),offspring.vec)))
      trio.list[[r]]=trio.obj
      names(trio.list)[r]=paste("trio",mom.trio[i],sep="")
      r = r +1
      }
    }
  # We turn the trio.flag of the mom.trio to TRUE to indicate they now have a trio
  trio.flag[mom.trio] = TRUE
  }
  # Record number of trios up to now
  cumultrio.vec[md-lev+1] = r - 1
  }
#  print(trio.flag)
# Depth 1
  # Create vector of unique moms and dads at the current depth
  currentdads = unique(dad.id[dv==1])
  currentmoms = unique(mom.id[dv==1])

  # All current dads are founders at depth 1. We give them all a trio
  dad.trio = currentdads
  # Loop over the non-founder dads at the current depth
  if (length(dad.trio)>0)
  {
  for (i in 1:length(dad.trio))
    {
    spousevec = unique(mom.id[dad.id==dad.trio[i]])
#    print (dad.trio[i])
#    print (spousevec)
    for (j in 1:length(spousevec))
      {
      # Create vector of children of the couple
      offspring.vec = id[dad.id==dad.trio[i] & mom.id==spousevec[j]]
      # Creating the elements of the trio object in trio.list
      trio.obj=list(id=dad.trio[i],spouse=spousevec[j],offspring = as.list(ifelse(trio.flag[offspring.vec],paste("trio",offspring.vec,sep=""),offspring.vec)))
      trio.list[[r]]=trio.obj
      names(trio.list)[r]=paste("trio",dad.trio[i],sep="")
      r = r +1
      }
    # We turn the trio.flag of the spouses to TRUE to indicate they are now members of a trio
    trio.flag[spousevec] = TRUE      
    }
#    print(trio.flag)
  }
  # Create vector of unique moms at the current depth who are not members of a trio
  mom.tokeep = !trio.flag[id%in%currentmoms]
#  print (mom.tokeep)
  mom.trio = names(mom.tokeep)[mom.tokeep]
  # Loop over the non-founder moms at the current depth
  if (length(mom.trio)>0)
  {
  for (i in 1:length(mom.trio))
    {
    spousevec = unique(dad.id[mom.id==mom.trio[i]])
    for (j in 1:length(spousevec))
      {
      # Create vector of children of the couple
      offspring.vec = id[mom.id==mom.trio[i] & dad.id==spousevec[j]]
      # Creating the elements of the trio object in trio.list
      trio.obj=list(id=mom.trio[i],spouse=spousevec[j],offspring = as.list(ifelse(trio.flag[offspring.vec],paste("trio",offspring.vec,sep=""),offspring.vec)))
      trio.list[[r]]=trio.obj
      names(trio.list)[r]=paste("trio",mom.trio[i],sep="")
      r = r +1
      }
    }
  # We turn the trio.flag of the mom.trio to TRUE to indicate they now have a trio
  trio.flag[mom.trio] = TRUE
  }
# Record the total number of trios
cumultrio.vec[md] = r - 1 
# Compute number of trios per level
ped.vec = diff(c(0,cumultrio.vec))

# from here down has been edited by sgy on Aug. 5, 2013

# Return the list of trios
#trio.list

ped2trio.list <- trio.list
trio.obj.list <- list()
k <- cumultrio.vec[md]

if( length(ped.vec) == 2 ){
  for( i in 1:ped.vec[1] ){
    trio.obj.list[[i]] <- new("Trio", id = ped2trio.list[[i]]$id, spouse = ped2trio.list[[i]]$spouse, offspring = ped2trio.list[[i]]$offspring)
  }
  trio.obj <- new("Trio", id = ped2trio.list[[k]]$id, spouse = ped2trio.list[[k]]$spouse, offspring = trio.obj.list[1:ped.vec[1]] )
  
}else if( length(ped.vec) == 3 ){
  for( i in 1:ped.vec[1] ){
    trio.obj.list[[i]] <- new("Trio", id = ped2trio.list[[i]]$id, spouse = ped2trio.list[[i]]$spouse, offspring = ped2trio.list[[i]]$offspring)
  }
  # there is a strange reversal of order here
  # without the reversal the subscript for the offspring would be i - ped.vec[2]
  for( i in ((1:ped.vec[2])+cumultrio.vec[1]) ){
    trio.obj.list[[i]] <- new("Trio", id = ped2trio.list[[i]]$id, spouse = ped2trio.list[[i]]$spouse, offspring = trio.obj.list[cumultrio.vec[1] + ped.vec[2] - i + 1 ])
  }
  trio.obj <- new("Trio", id = ped2trio.list[[k]]$id, spouse = ped2trio.list[[k]]$spouse, offspring = trio.obj.list[cumultrio.vec[1] + 1:ped.vec[2] ] )
}
return(list(trio.list = ped2trio.list, obj.list = trio.obj.list, object = trio.obj))
}
