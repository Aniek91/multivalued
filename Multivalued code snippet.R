####################################
### User defined functions rpart
### Minimizes weighted misclassification error 
### Three functions: iloss, eloss and sloss

library(combinat)
library(rpart)

# The init function:
#   fix up y to deal with offsets
#   return a dummy parms list
#   numresp is the number of values produced by the eval routine's "label"
#   numy is the number of columns for y
#   summary is a function used to print one line in summary.rpart
# In general, this function would also check for bad data, see rpart.poisson
#   for instace.
# wt and y have length equal to number of nodes
iloss <- function(y, offset, parms, wt) {   
  if (!is.null(offset)) y <- y-offset
  sfun <- function(yval, dev, wt, ylevel, digits) {
    
    paste("  predicted class=", format(yval),        #Lisa: so far, only predicted class is printed in summary of nodes. 
          #         ", expected loss=" , format(length(wt)),   #More results can be added here.
          sep='')
  }
  list(y=y, parms=0, numresp=1, numy=1, summary=sfun)
}


# The 'evaluation' function.  Called once per node.
#  Produce a label (1 or more elements long) for labeling each node,
#  and a deviance.  The latter is
#	 - of length 1
#  - equal to 0 if the node is "pure" in some sense (unsplittable)
#  - does not need to be a deviance: any measure that gets larger
#    as the node is less acceptable is fine.
#  - the measure underlies cost-complexity pruning, however
eloss <- function(y, wt, parms){
  indiv.loss.node <- indiv.loss[wt,]
  classes <- 1:ncol(indiv.loss) 
  loss.per.class <- apply(indiv.loss.node, 2, sum)
  yopt <- classes[which.min(loss.per.class)]
  wme <- min(loss.per.class)                      #weighted misclassification error
  # print(yopt)                                     
  # print(wt)                                       #print to check lenght and values of arguments
  # print(y)
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
#       direction = a vector of length m= "# groups".  It asserts that the               !!! Lisa: work around this
#           best split can be found by lining the groups up in this order
#           and going from left to right, so that only m-1 splits need to
#           be evaluated rather than 2^(m-1)
#       goodness = m-1 values, as before.
#
# The reason for returning a vector of goodness is that the C routine
#   enforces the "minbucket" constraint. It selects the best return value
#   that is not too close to an edge.
sloss <- function(y, wt, x, parms, continuous)
{
  n <- length(y)
  # print(n)
  # print(x)             # a vector with one predictor per time (ordered)
  # print(y)
  # print(wt)
  
  #Lisa: Select rows of indiv.losses that are in node
  indiv.loss.node <- indiv.loss[wt,]
  
  classes <- unique(y)                             #classes in parent node
  loss.per.class <- apply(indiv.loss.node, 2, sum) #misclassification error in parent node if node is assigned to class 1, 2, 3, ..., respectively
  parent.error <- min(loss.per.class)              #minimum misclassification error in parent node
  goodness <- vector()
  direction <- vector()
  
  if (continuous) {
    # continuous x variable
    for(i in 1:(n-1)){
      left.loss <- matrix(indiv.loss.node[1:i,], ncol=ncol(indiv.loss))       #select individuals in left child node
      left.loss.per.class <- apply(left.loss, 2, sum)                         #misclassification error in left child node if node is assigned to class 1, 2, 3, ..., respectively
      left.yopt <- classes[which.min(left.loss.per.class)]                    #class which results in minimum misclassification error
      left.error <- min(left.loss.per.class)                                  #minimum misclassification error in left child node
      
      right.loss <-  matrix(indiv.loss.node[(i+1):n,], ncol=ncol(indiv.loss)) #select individuals in right child node
      right.loss.per.class <- apply(right.loss, 2, sum)                       
      right.yopt <- classes[which.min(right.loss.per.class)]
      right.error <-min(right.loss.per.class)
      
      goodness[i] <- parent.error-left.error-right.error                      #decrease in misclassification error
      direction[i] <- sign(-left.error-exp(-100))
    }
    list(goodness = goodness, direction = direction)
  } else {
    # Categorical X variable  
    ux <- sort(unique(x))
    
    goodness.list <- list()
    list.permutations <- permn(as.numeric(ux))                                #all permutations of the elements of x
    list.permutations <- list.permutations[1:(length(list.permutations)/2)]   #permutations without reverse duplicates
    for(j in 1:length(list.permutations)){
      ord <- list.permutations[[j]]                                           #select j'th permutation
      n <- length(ord)                                                        #number of categories in x
      for(i in 1:(n-1)){                                                      #for every split in categories ... 
        left.loss <- matrix(indiv.loss.node[x%in%ux[ord[1:i]],], ncol=ncol(indiv.loss))       #select individuals that go to the left           
        left.loss.per.class <- apply(left.loss, 2, sum)                                       #same code as for non-categorical ...
        left.yopt <- classes[which.min(left.loss.per.class)]
        left.error <- min(left.loss.per.class)
        
        right.loss <-  matrix(indiv.loss.node[x%in%ux[ord[(i+1):n]],], ncol=ncol(indiv.loss))  
        right.loss.per.class <- apply(right.loss, 2, sum)
        right.yopt <- classes[which.min(right.loss.per.class)]
        right.error <-min(right.loss.per.class)
        
        goodness[i] <- signif(parent.error-left.error-right.error)            #decrease in misclassification error
      }
      goodness.list[[length(goodness.list)+1]] <- goodness                    
    }
    
    best.permn <-  which.max(lapply(goodness.list, max))                      #select best permutation
    goodness <- goodness.list[[best.permn]]                                   #(number of categories - 1) values (i.e., for every possible split of selected permutation)
    ord <- list.permutations[[best.permn]]                                    #best ordering of categories in x    
    
    list(goodness= goodness,direction = ux[ord])                              #direction = number of categories
  } 
}


#########################################
### Rpart with user defined functions ###
#########################################
alist <- list(eval = eloss, split = sloss, init = iloss)       #list of user defined functions, to be used as method in rpart
