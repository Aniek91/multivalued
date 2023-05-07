#Supplementary Materials for multirisk"

# Lisa L. Doove, Elise Dusseldorp, and Iven Van Mechelen

#In this Supplementary Material, we present code for the splitting functions
#that can be added as a set of R functions into the rpart package, and which
#implement the splitting method described in "Extension of a generic approach
#for estimating optimal treatment regimes to multivalued treatments".

library(combinat)

# The init function:
#   fix up y to deal with offsets
#   return a dummy parms list
#   numresp is the number of values produced by the eval routine's "label"
#   numy is the number of columns for y
#   summary is a function used to print one line in summary.rpart
# In general, this function would also check for bad data, see rpart.poisson
#   for instace.
# wt and y have length equal to number of nodes
mr_init <- function(y, offset, parms, wt) {   
  if (!is.null(offset)) y <- y-offset
  sfun <- function(yval, dev, wt, ylevel, digits) {
    
    paste("predicted class=", format(yval), sep='')
  }
  environment(sfun) <- .GlobalEnv
  list(y=y, parms=parms, numresp=1, numy=1, summary=sfun)
}

# The 'evaluation' function.  Called once per node.
#  Produce a label (1 or more elements long) for labeling each node,
#  and a deviance.  The latter is
#	- of length 1
#       - equal to 0 if the node is "pure" in some sense (unsplittable)
#       - does not need to be a deviance: any measure that gets larger
#            as the node is less acceptable is fine.
#       - the measure underlies cost-complexity pruning, however
mr_eval <- function(y, wt, parms){
  indiv.loss.node <- parms$indiv.loss[wt,, drop=FALSE]
  classes <- 1:ncol(parms$indiv.loss) #sort(unique(y))
  loss.per.class <- apply(indiv.loss.node, 2, sum)
  yopt <- classes[which.min(loss.per.class)]
  wme <- min(loss.per.class) 
  list(label = yopt, deviance = wme)
}

# The split function, where most of the work occurs.
#   Called once per split variable per node.
# If continuous=T
#   The actual x variable is ordered
#   y is supplied in the sort order of x, with no missings,
#   return two vectors of length (n-1):
#      goodness = goodness of the split, larger numbers are better.
#                 0 = couldn't find any worthwhile split
#        the ith value of goodness evaluates splitting obs 1:i vs (i+1):n
#      direction= -1 = send "y< cutpoint" to the left side of the tree
#                  1 = send "y< cutpoint" to the right
#         this is not a big deal, but making larger "mean y's" move towards
#         the right of the tree, as we do here, seems to make it easier to
#         read
# If continuos=F, x is a set of integers defining the groups for an
#   unordered predictor.  In this case:
#       direction = a vector of length m= "# groups".  It asserts that the               !!!
#           best split can be found by lining the groups up in this order
#           and going from left to right, so that only m-1 splits need to
#           be evaluated rather than 2^(m-1)
#       goodness = m-1 values, as before.
#
# The reason for returning a vector of goodness is that the C routine
#   enforces the "minbucket" constraint. It selects the best return value
#   that is not too close to an edge.
mr_split <- function(y, wt, x, parms, continuous)
{
  n <- length(y)
  indiv.loss.node <- parms$indiv.loss[wt,, drop=FALSE]
  classes <- unique(y)
  loss.per.class <- apply(indiv.loss.node, 2, sum)
  parent.error <- min(loss.per.class)
  
  if(continuous) {
    # continuous x variable
    goodness <- numeric(n-1)
    direction <- numeric(n-1)
    for(i in 1:(n-1)){
      left.loss <- indiv.loss.node[1:i,,drop=FALSE] 
      left.loss.per.class <- apply(left.loss, 2, sum)
      left.yopt <- classes[which.min(left.loss.per.class)]
      left.error <- min(left.loss.per.class)
      
      right.loss <-  indiv.loss.node[(i+1):n,,drop=FALSE] 
      right.loss.per.class <- apply(right.loss, 2, sum)
      right.yopt <- classes[which.min(right.loss.per.class)]
      right.error <-min(right.loss.per.class)
      
      goodness[i] <- parent.error-left.error-right.error
      direction[i] <- sign(-left.error-exp(-100))
    }
    list(goodness = goodness, direction = direction)
  } else {
    # Categorical X variable  
    ux <- sort(unique(x))
    goodness <- numeric(length(ux)-1)
    ux.seq <- seq(ux)
    
    goodness.list <- list()
    list.permn <- permn(ux.seq)
    list.permn <- list.permn[1:(length(list.permn)/2)] 
    
    for(j in 1:length(list.permn)){
      ord <- list.permn[[j]]
      n <- length(ord)      #number of categories in x
      for(i in 1:(n-1)){
        left.loss <- indiv.loss.node[x%in%ux[ord[1:i]],,drop=FALSE]    
        left.loss.per.class <- apply(left.loss, 2, sum)
        left.yopt <- classes[which.min(left.loss.per.class)]
        left.error <- min(left.loss.per.class)
        
        right.loss <-  indiv.loss.node[x%in%ux[ord[(i+1):n]],,drop=FALSE] 
        right.loss.per.class <- apply(right.loss, 2, sum)
        right.yopt <- classes[which.min(right.loss.per.class)]
        right.error <-min(right.loss.per.class)
        
        goodness[i] <- signif(parent.error-left.error-right.error) 
      }
      goodness.list[[length(goodness.list)+1]] <- goodness
    }
    
    best.permn <-  which.max(lapply(goodness.list, max))
    goodness <- goodness.list[[best.permn]]
    ord <- list.permn[[best.permn]]
    
    list(goodness= goodness,direction = ux[ord])  #direction = # categories
  } 
}

#create method multirisk
multirisk <- list(eval = mr_eval, split = mr_split, init = mr_init) 