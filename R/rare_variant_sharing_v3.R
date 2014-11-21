# By Alexandre Bureau
# 2013/06/28

RVsharing.fn = function(id, dad.id, mom.id,carriers)
{

N = length(id)
# vector of indices of final descendants
fdi = which(!(id%in%dad.id | id%in%mom.id))
nfd = length(fdi)
if (nfd < 2) stop("There are fewer than 2 descendants for which to compute a rare variant sharing probability.")
if (!is.null(carriers))
{
	missc = setdiff(carriers,id)
	if(length(missc)>0) stop(missc," not in pedigree.")
	}

# Getting the depth of each subject in the pedigree
dv = kindepth(id, dad.id, mom.id)
md = max(dv)
# Number of founders
Nf = sum(dv==0)

# Collecting the degrees between the sequenced children, the founders and the intermediate ancestors
degvec = numeric(nfd)
currentnonfounders = currentfounders = character(nfd)
active = rep(TRUE,nfd) 
# List of distance to founders of each final descendant
desfounders = list()
# List of list of final descendants and intermediate ancestors below each founder below each ancestor
foundersdegreedes = list()
# List of indicators of whether each descendant of a founder is an intermediate ancestor
iancestor.as.descendant = list()
# List of intermediate ancestors
iancestors = character(0)
# Number of founders below each intermediate ancestor
iancestors.Nf = numeric(0)
# List of vectors of degrees of relatedness of final descendants below each intermediate ancestor
ancestorsdegreedes = list()
# List of intermediate ancestors on the current pedigree degree
lev.ia = list()
ia = lia = 1

# Initializing the currentnonfounders vector with the final descendants at the deepest level
currentnonfounders[dv[fdi]==md] = id[dv==md]

# Loop from highest to lowest depth
if (md > 1)
{
for (lev in (md-1):1)
  {
  # Incrementing D for final descendants and intermediate ancestors
  degvec[dv[fdi]>lev&active] = degvec[dv[fdi]>lev&active] + 1
  # Listing ancestors at current depth of final descendants and intermediate ancestors
  # Loop over final descendants and intermediate ancestors with depth greater than current depth
  for (i in (1:nfd)[dv[fdi]>lev&active])
    {
    # The currentnonfounders are those from the previous level
    currentdad=dad.id[id == currentnonfounders[i]]
    currentmom=mom.id[id == currentnonfounders[i]]  
    # Identify non-founder among mom and dad
    currentnonfounders[i] = ifelse(is.na(dad.id[id == currentdad]),currentmom,currentdad)
    # Identify founder among mom and dad
    currentfounders[i] = ifelse(is.na(dad.id[id == currentdad]),currentdad,currentmom)
    }
  # Adding final descendent at the current level to the currentnonfounders for the next level
  currentnonfounders[dv[fdi]==lev] = id[fdi][dv[fdi]==lev]

  # Checking if there is an intermediate ancestor with more than one descendant at the current level
  tab.currentnonfounders = table(currentnonfounders[currentnonfounders>0&active])
  # If there is more than one, stop because it is not implemented
  # if (sum(tab.currentnonfounders>1) > 1) stop ("More than one intermediate ancestor at the same level with two or more descendant")
  # If there is any intermediate ancestor with more than one descendant at the current level
  if (any(tab.currentnonfounders>1))
    {
    lev.ia[[lia]] = names(tab.currentnonfounders)[tab.currentnonfounders>1]
    iancestors = c(iancestors,lev.ia[[lia]])
    for (i in which(iancestors %in% lev.ia[[lia]]))
    {
    # Adding the degrees of final descendants below the current intermediate ancestor to his list
    ancestorsdegreedes[[i]] = degvec[currentnonfounders==iancestors[i]&active]
    names(ancestorsdegreedes[[i]]) = id[fdi][currentnonfounders==iancestors[i]&active]
    # Include degrees between spouse and final descendants in list of ancestors  (assumes only one spouse)
    foundersdegreedes[[i]] = list(degvec[currentnonfounders==iancestors[i]&active])
    # Setting indicator of whether the descendant is the previous intermediate ancestor
    if (ia>1)
      {
      if (length(lev.ia[[lia]])>1) stop ("More than one intermediate ancestor at level ",lev," with two or more descendants.")
      iancestor.as.descendant[[i]] = list(id[fdi][currentnonfounders==iancestors[i]&active] == iancestors[ia-1])
      # Recording number of founders below intermediate ancestor
      }
    else iancestor.as.descendant[[i]] = list(rep(FALSE,length(foundersdegreedes[[i]][[1]])))
    # Include previous ancestors of final descendants if any
    if (any(names(desfounders) %in% id[fdi][currentnonfounders==iancestors[i]&active]))
      {
      ii = 1
      tmp = desfounders[names(desfounders) %in% id[fdi][currentnonfounders==iancestors[i]&active]]
      # Loop over final descendants 
      # check whether this is correct with more than one intermediate ancestors on the same level
      for (k in 1:length(tmp))
        {
        foundersdegreedes[[i]][(ii+1):(ii+length(tmp[[k]]))] = tmp[[k]]
        # Setting indicator of whether the descendant of all ancestors in tmp[[k]] is the previous intermediate ancestor
        iancestor.as.descendant[[i]][(ii+1):(ii+length(tmp[[k]]))] = list(ifelse (ia>1, names(tmp)[k] == iancestors[ia-1], FALSE))
        ii = ii + length(tmp[[k]])
        # Adding spouse of intermediate ancestor to list of founders of current final descendant
        
        }
      }
      iancestors.Nf[i] = ifelse(any(unlist(iancestor.as.descendant[[i]])),iancestors.Nf[ia-1],0) + length(iancestor.as.descendant[[i]])
      }
    }
  # Adding the current founder ancestral to each final descendants to his list of founders
  if (length(currentfounders[currentfounders>0][!is.na(currentfounders)])>0)
    {
    for (i in (1:nfd)[dv[fdi]>lev&active])
      {
      # If there are at least i elements in desfounders
      if (length(desfounders)>=i)
        {
        desfounders[[i]][length(desfounders[[i]])+1] = degvec[i]
        # Keeping the name of the founder
        names(desfounders[[i]])[length(desfounders[[i]])] = currentfounders[i]
        }
      else 
        {
        desfounders[[i]] = list(degvec[i])
        names(desfounders[[i]])[1] = currentfounders[i]
        }
      # Assigning the ID of the subject as name
      names(desfounders)[i] = id[fdi][i]
      }
    }
  # Finishing processing the current intermediate ancestor if there is one
  if (any(tab.currentnonfounders>1))
    {
    for (i in which(iancestors %in% lev.ia[[lia]]))
    {
    # Turning these final descendants to inactive
    active[currentnonfounders==iancestors[i]] = FALSE
    # Removing spouse(s) of intermediate ancestor from currentfounders 
    # Note: the spouse(s) have the same positions in the currentfounders vector as the 
    # intermediate ancestor in the currentnonfounders vector
    currentfounders = currentfounders[currentnonfounders != iancestors[i]]
    # Adding the intermediate ancestor to the vector of subjects with a degree
    if (any(id==iancestors[i]))
    {
    nfd = nfd + 1
    fdi[nfd] = which(id==iancestors[i])
    degvec[nfd] = 0
    active[nfd] = TRUE
    # Adding the intermediate ancestor to the vector of currentnonfounders
    currentnonfounders[nfd] = iancestors[i]
    }
    }
    # Incrementing ia
    ia = ia + length(lev.ia[[lia]])
    lia = lia + 1
    }
  }
}
# Depth 0: there should be at most 2 founders common to all subjects
# We assign one of them as a dummy "intermediate" ancestor
# Incrementing D for final descendants and intermediate ancestors
degvec[active] = degvec[active] + 1
# Listing ancestors at current depth of final descendants and intermediate ancestors
# The currentnonfounders are those from the previous level
currentdads=dad.id[id %in% currentnonfounders[active]]
currentmoms=mom.id[id %in% currentnonfounders[active]]

# The dummy intermediate ancestor has all founders below him (except himself)
iancestors.Nf[ia] = Nf - 1 
    
# If all subjects have the same dad, use him as last ancestor
if (all(currentdads==currentdads[1]))
  { 
  iancestors[ia] = currentdads[1]
  currentfounders = currentmoms
  }
# else if all subjects have the same mom, use her as last ancestor
else
  {
  if (all(currentmoms==currentmoms[1])) 
    {
    iancestors[ia] = currentmoms[1]
    currentfounders = currentdads
    }  
# else there is no common ancestor, and the probability of sharing is 0
  else return (0)
  }
     
    # Adding the degrees of final descendants below the current intermediate ancestor to his list
    ancestorsdegreedes[[ia]] = degvec[active]
    names(ancestorsdegreedes[[ia]]) = id[fdi][active]
    # Include first spouse in list of ancestors
    spousevec = unique(currentfounders)
    foundersdegreedes[[ia]]= list(degvec[currentfounders==spousevec[1]&active])
    # Setting indicator of whether the descendant is the previous intermediate ancestor
    if (ia>1)
      iancestor.as.descendant[[ia]] = list(id[fdi][currentfounders==spousevec[1]&active] %in% lev.ia[[lia-1]])
    else iancestor.as.descendant[[ia]] = list(rep(FALSE,length(foundersdegreedes[[ia]][[1]])))
    # Add additional spouses if any
    # Warning! This is going to work only if all previous intermediate descendents are under the same spouse
    if(length(spousevec)>1)
      {
      for (i in 2:length(spousevec))
        {
        foundersdegreedes[[ia]][[i]] = degvec[currentfounders==spousevec[i]&active]
        if (ia>1)
          iancestor.as.descendant[[ia]][[i]] = id[fdi][currentfounders==spousevec[i]&active] %in% lev.ia[[lia-1]]
        else iancestor.as.descendant[[ia]][[i]] = rep(FALSE,length(foundersdegreedes[[ia]][[i]]))
        }
      }  
    # Include previous ancestors of final descendants if any
    if (any(names(desfounders) %in% id[fdi][active]))
      {
      ii = length(spousevec)
      tmp = desfounders[names(desfounders) %in% id[fdi][active]]
      # Loop over final descendants 
      for (k in 1:length(tmp))
        {
        foundersdegreedes[[ia]][(ii+1):(ii+length(tmp[[k]]))] = tmp[[k]]
        # Setting indicator of whether the descendant of all ancestors in tmp[[k]] is a previous intermediate ancestor
        iancestor.as.descendant[[ia]][(ii+1):(ii+length(tmp[[k]]))] = list(ifelse (ia>1, names(tmp)[k] %in% iancestors[lev.ia[[lia-1]]], FALSE))
        ii = ii + length(tmp[[k]])
        }
      }
      # Adding the current founder couple ancestral to each final descendants to his list of founders
      # This is not required for the sharing probability computation, but is used for kinship estimation
    # print(currentfounders)
    for (i in (1:nfd)[active])
      {
      j = 1
      # If there are at least i elements in desfounders
      if (length(desfounders)>=i)
        {
        desfounders[[i]][length(desfounders[[i]])+(1:2)] = degvec[i]
        # Keeping the name of the founder
        names(desfounders[[i]])[length(desfounders[[i]])-1] = currentfounders[j]
        names(desfounders[[i]])[length(desfounders[[i]])] = iancestors[ia]
        }
      else 
        {
        desfounders[[i]] = rep(degvec[i],2)
        names(desfounders[[i]])[1] = currentfounders[j]
        names(desfounders[[i]])[2] = iancestors[ia]
        }
      # Assigning the ID of the subject as name
      names(desfounders)[i] = id[fdi][i]
      j = j+1
      }

# Computation of numerator
num = 1
for (i in 1:ia)
  num = num * 1/2^sum(ancestorsdegreedes[[i]])
# Computation for top founder or founders
# If there is only one spouse, then a couple of founders can transmit a variant to all final descendents
if (length(spousevec)==1) num = num*2
# Division by the number of founders
num = num/Nf
 
# Computation of denominator
# Probability that no variant has been transmitted
p0 = 0
# Probability that no variant has been transmitted from previous intermediate ancestor
pk = 1
for (i in 1:ia)
  {
  for (j in 1:length(foundersdegreedes[[i]]))
    p0 = p0 + prod((1-1/2^foundersdegreedes[[i]][[j]]) + ifelse(iancestor.as.descendant[[i]][[j]],(1/2^foundersdegreedes[[i]][[j]])*pk,0))
  # Updating the probability for the previous intermediate ancestor, who becomes the current intermediate ancestor
  # For now, intermediate ancestors can have only one spouse, this is why we take the indicators of the first founder attached to him
  if (i<ia) pk = prod((1-1/2^ancestorsdegreedes[[i]]) + ifelse(iancestor.as.descendant[[i]][[1]],1/2^ancestorsdegreedes[[i]]*pk,0))
  }
# At the end, add the probability from the dummy "intermediate" ancestor. He is currently the only one who can have more than one spouse
# Since only one of his spouses can be the parent of the previous intermediate ancestor, sapply returns only one non-zero term.
# The summation returns in fact the value of that single non-zero term
  # Debugging code
  # print (foundersdegreedes[[i]])
  # print (ancestorsdegreedes[[i]])
  # print (iancestor.as.descendant[[i]])
  # print (p0)
  # print (sapply(iancestor.as.descendant[[i]][1:length(spousevec)],function(lv,deg,pk) ifelse(lv, (1/2^deg) * pk,0), deg=ancestorsdegreedes[[i]],pk=pk))
# p0 = p0 + prod((1-1/2^ancestorsdegreedes[[i]]) + sum(sapply(iancestor.as.descendant[[i]][1:length(spousevec)],function(lv,deg,pk) ifelse(lv, (1/2^deg) * pk,0), deg=ancestorsdegreedes[[i]],pk=pk)))
# This remains to be tested with >1 spouse
  tmpf = as.matrix(sapply(iancestor.as.descendant[[i]][1:length(spousevec)],function(lv,deg,pk) ifelse(lv, (1/2^deg) * pk,0), deg=ancestorsdegreedes[[i]],pk=pk))  
  p0 = p0 + prod((1-1/2^ancestorsdegreedes[[i]]) + apply(tmpf,1,sum))
    # Debugging code
    # print (p0)
if (is.null(carriers))
# Sharing probability
pshare = num/(1-p0/Nf)
else
{
#  ci = which(id %in% carriers)
  noncarriers = setdiff(id[!(id%in%dad.id | id%in%mom.id)],carriers)
  fd.subsets = list(as.matrix(carriers))
  # Loop over number of non-carriers to include as "carrier" in the possible subset
  if (length(noncarriers)>1)
    for (k in 1:(length(noncarriers)-1))
    {
  	  tmp = combn(noncarriers,k)
  	  fd.subsets[[k+1]] = rbind(matrix(carriers,length(carriers),ncol(tmp)),tmp)
    }
  subsetp = numeric(length(fd.subsets))
  # Loop over possible subsets
  for (k in 1:length(fd.subsets))	
    {
    sn = nrow(fd.subsets[[k]])
    nsubs = ncol(fd.subsets[[k]])
    subsetkp = numeric(nsubs)
  	for (h in 1:nsubs)
  	  {
		fremoved = 0
		insubset = logical(sn)
  	    # Loop over intermediate ancestors
  	    for (i in 1:ia)
  	    {
		    insubset[fd.subsets[[k]][,h]%in%names(ancestorsdegreedes[[i]])] = TRUE
		    # Check if all carriers in current possible subset are descendents of current intermediate ancestor 
 	  		#if (all(fd.subsets[[k]][,h]%in%names(ancestorsdegreedes[[i]])))
 	  		if (all(insubset))
 	  		{
 	  		   if (i == 1) break
 	  		   else if(iancestors[i-1] %in% names(ancestorsdegreedes[[i]])) break
			}
		 }
		 # Compute probability of subset
		 numsub = 1
		 for (ii in 1:i)
		 {
		      numsub = numsub * 1/2^sum(ancestorsdegreedes[[ii]][c(fd.subsets[[k]][,h],iancestors)],na.rm=TRUE)
		 	  fdn.vec = setdiff(names(ancestorsdegreedes[[ii]]),fd.subsets[[k]][,h])
		 	  if (length(fdn.vec)>0)
		 	    for (fd in fdn.vec)
		 	    {
		 	  	ncf = names(desfounders[[fd]])
				founder.for.other = logical(length(ncf))
				# check if any founder of the non-carrier is also a founder for another non-carrier
				for (z in names(desfounders[names(desfounders)!=fd])) founder.for.other = pmax(founder.for.other,ncf%in%names(desfounders[[z]]))
				# Count number of founders unique to final descendant fd to remove them
				fremoved = fremoved + sum(founder.for.other==0,na.rm=TRUE)
				}
		 }
		 # Multiply by two for the spouses and divide by the number of founders
		 subsetkp[h] = numsub*2/(iancestors.Nf[i]-fremoved+1)
  	  }
  	  subsetp[k] = sum(subsetkp)
  	}
  # Add joint prob for all final descendents
  subsetp = c(subsetp,num)
  # Computation of sharing probability of observed subset
  numo = sum(subsetp*(-1)^(0:(length(subsetp)-1)))
  pshare = numo/(1-p0/Nf)   
}
new("RVsharingProb",pshare=pshare,iancestors=iancestors,desfounders=desfounders,id=as.character(id),dad.id=as.character(dad.id),mom.id=as.character(mom.id),carriers=as.character(carriers))
}

# Wrappers for pedigree object
# Returns only pshare
RVsharing.ped.pshare = function(ped)
{
RVsharing(ped)@pshare
} 
