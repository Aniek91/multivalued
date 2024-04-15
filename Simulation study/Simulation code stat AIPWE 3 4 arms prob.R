library(rpart) 

##########################
### Summary of results ###
##########################

descript.rpart <- function(x, scenario, arms){
  meanOutcome <- round(mean(unlist(x["Outcome",]),na.rm =T),2)
  stdOutcome <- round(sd(unlist(x["Outcome",]),na.rm =T),2)
  meanExpOutcome <- round(mean(unlist(x["ExpOut",]),na.rm =T),2)
  stdExpOutcome <- round(sd(unlist(x["ExpOut",]),na.rm =T),2)
  
  #Type I/II error
  treecreated <- unlist(x["tree.created",])
  p.treecreated <- round(mean(unlist(x["tree.created",]),na.rm =T),2)
  
  #Recovery
  if(arms==3){                                                                             ###fourth scenario added
    truesize <- c(2,3,1,3,3)
    truevar <- list(c("X1", "X2"), c("X1", "X2", "X3"), c("X1"),c("X1", "X2", "X3"),c("X1", "X2", "X3"))
    truesplit <- list(c(-.545, 0), c(0,0,0), c(0), c(0,0,0), c(0,0,0))}                               
  if(arms==4){
    truesize <- c(3,4,1,4,4)
    truevar <- list(c("X1", "X2", "X3"), c("X1", "X2", "X3", "X4"), c("X1"), c("X1", "X2", "X3", "X4"), c("X1", "X2", "X3", "X4"))
    truesplit <- list(c(-.745, -.545, 0), c(0,0,0,0), c(0), c(0,0,0,0), c(0,0,0,0))}
  
  p.correctsize <- mean(unlist(x["nsplits",treecreated])==truesize[scenario])
  correctsize <- round(unlist(x["nsplits",])==truesize[scenario],2)
  p.correctvar <- round(mean(sapply(c(1:ncol(x))[correctsize], function(v) (sum(x[,v]$splitvar%in%truevar[[scenario]])==truesize[scenario])&(sum(truevar[[scenario]]%in%x[,v]$splitvar)==length(truevar[[scenario]])))),2)
  correctvar <- (sapply(c(1:ncol(x)), function(it) sum(x[,it]$splitvar%in%truevar[[scenario]])==truesize[scenario]))*correctsize
  p.correctsplit <- round(mean(sapply(c(1:ncol(x))[correctvar==1], function(s) in_interval(x[,s]$splitpoint, truesplit[[scenario]], 0.05))),2)
  
  Tpred <- matrix(unlist(x["assignment",]),ncol=ncol(x))
  Topt <-  matrix(unlist(x["opt.assignment",]),ncol=ncol(x))
  correctassigned <- (Tpred-Topt)^2 #0 if correctly assigned
  p.incorrect <- round(mean(apply(correctassigned, 2, mean)),2)
  
  biasSample <- round(mean(unlist(x["biassample",]),na.rm =T),2)
  
  return(c(egopt=meanOutcome, sd=stdOutcome, Egopt=meanExpOutcome, SD=stdExpOutcome, p.quali=p.treecreated, 
           p.size=p.correctsize, p.var=p.correctvar, p.split=p.correctsplit, p.incorrentassign=p.incorrect, biasSample=biasSample))
}


in_interval <- function(splits, true, margin){
  interval <- matrix(c(true-margin, true+margin), ncol=2)
  check <-  matrix(0, ncol=length(splits), nrow=nrow(interval))
  for(h in 1:length(splits)){
    for(i in 1:nrow(interval)){
      check[i,h] <- interval[i,1] < splits[h] & splits[h] < interval[i,2]  
    }}
  correctsplitpoint <- sum(apply(check, 2, sum)>0)==ncol(check) #sum over the columns (splitpoints)
  return(correctsplitpoint)
}



#####################
### Caipwe 3 arms ###
#####################

C.aipwe.prob <- function(arms, scenario, n.arm, meandif,model,a,dataAipweE, h){  
  set.seed(h)
  cat("Iteration/seed", h, "\n")
  data.generated <- datagen3.4(arms, scenario, n.arm, meandif,model,h)
  dataAipwe <- data.generated[[1]]
  Topt <- data.generated[[2]] ###optimal assignment
  dataAipwe[,"A"] <- as.factor(dataAipwe[,"A"])
  A <- dataAipwe[,"A"]    ###
  Y <- dataAipwe[,"Y"]    ###
  N <- arms*n.arm
  
  if(arms==3){
    Y.pop <- dataAipwe[,c("Y1","Y2","Y3")]
    OutcomeOpt <- mean(sapply(1:nrow(Y.pop), function(k) Y.pop[k,Topt[k]]))
  
  #Outcome model (exactly the same as estimating one model!)
  model1 <- lm(Y ~ 1 + X1 + X2 + X3 + X4 + X5, data=dataAipwe, subset=A==1)
  model2 <- lm(Y ~ 1 + X1 + X2 + X3 + X4 + X5, data=dataAipwe, subset=A==2)
  model3 <- lm(Y ~ 1 + X1 + X2 + X3 + X4 + X5, data=dataAipwe, subset=A==3)
  u1 <- predict.glm(model1, dataAipwe)
  u2 <- predict.glm(model2, dataAipwe)
  u3 <- predict.glm(model3, dataAipwe)
  est.outcomes <- cbind(u1,u2,u3)
  
  ##Propensity model
  PRaipwe <- rep(0.5, n.arm)
  
  #estimate contrast function
  c1 <- (Y*(A==1)/PRaipwe[1]) - (((A==1) - PRaipwe[1]) / PRaipwe[1]*u1)
  c2 <- (Y*(A==2)/PRaipwe[2]) - (((A==2) - PRaipwe[2]) / PRaipwe[2]*u2)
  c3 <- (Y*(A==3)/PRaipwe[3]) - (((A==3) - PRaipwe[3]) / PRaipwe[3]*u3)
  est.aipwe <- cbind(c1,c2,c3)
  
  #response and losses
  Zaipwe <- sapply(1:nrow(dataAipwe), function(b) which.max(est.aipwe[b,]))
  dataCaipwe <- cbind(as.factor(Zaipwe), dataAipwe[,3:7])
 
  #CART
  controlrpart <- rpart.control(maxdepth=6, maxsurrogate=0, maxcompete=0)
  
  indiv.loss=t(sapply(1:nrow(dataCaipwe), function(b) max(est.aipwe[b,])-est.aipwe[b,]))
  
  
  fit <-  rpart(as.numeric(Zaipwe) ~ X1 + X2 + X3 + X4 + X5, 
                data = dataCaipwe, weights=1:N, control=controlrpart, method=multirisk,
                parms=list(indiv.loss=t(sapply(1:nrow(dataCaipwe), function(b)   
                  max(est.aipwe[b,])-est.aipwe[b,])))
                )                                                    
  
  
   options(warn=-1)
  
    if(nrow(fit$frame)>1){
      xfit <- xpred.rpart(fit,xval=20)
      xerror <- colMeans(sapply(1:ncol(xfit), function(m) sapply(1:nrow(xfit),    
                                                                 function(k) indiv.loss[k,xfit[k,m]])))
      cp1 <- fit$cptable[which.min(xerror),1]
      fit <- prune(fit, cp=cp1)
    }
    
    if(!nrow(fit$frame)>1){
      Outcome <- max(aggregate(dataAipwe$Y, list(dataAipwe$A), mean)[,2])
      Ag <- which.max(aggregate(dataAipwe$Y, list(dataAipwe$A), mean)[,2])
      assignment <- rep(Ag,N)
      ExpOut <- mean(dataAipweE[,c("Y1","Y2","Y3")][,Ag])
      nsplits <-0
      splitvar <- FALSE
      splitpoint <- FALSE
      pred.Y <- rep(Outcome,N)
      tree<- NA 
      tree.created <- F
    }else{print(fit)} 
    
    if(nrow(fit$frame)>1){     # AIPWE schatter voor expected outcome
      noderow <- fit$where
      assignment <- fit$frame[noderow,5]
      indicator <- assignment==A
      m <- sapply(1:N, function(v) est.outcomes[v,assignment[v]])
      propensity <- PRaipwe[assignment]
      Outcome <- mean(((Y*indicator)/propensity) - ((indicator-propensity)/propensity*m))
      
      nsplits <- sum(fit$frame[,"var"]!="<leaf>")
      splitvar <-  rownames(fit$splits)
      splitpoint <- fit$splits[,"index"]
      tree <- list(fit$frame, fit$where, fit$splits)
      tree.created <- T
      
      #estimate expected outcome in population under regime
      Agexp <- predict(fit, dataAipweE)
      Y.true <- dataAipweE[,c("Y1","Y2","Y3")]
      ExpOut <- mean(sapply(1:length(Agexp), function(s) Y.true[s, Agexp[s]]))
    }
  }
  
  if(arms==4){
    Y.pop <- dataAipwe[,c("Y1","Y2","Y3","Y4")]
    OutcomeOpt <- mean(sapply(1:nrow(Y.pop), function(k) Y.pop[k,Topt[k]]))
    
    #Outcome model (exactly the same as estimating one model!)
    model1 <- lm(Y ~ 1 + X1 + X2 + X3 + X4 + X5, data=dataAipwe, subset=A==1)
    model2 <- lm(Y ~ 1 + X1 + X2 + X3 + X4 + X5, data=dataAipwe, subset=A==2)
    model3 <- lm(Y ~ 1 + X1 + X2 + X3 + X4 + X5, data=dataAipwe, subset=A==3)
    model4 <- lm(Y ~ 1 + X1 + X2 + X3 + X4 + X5, data=dataAipwe, subset=A==4)
    u1 <- predict.glm(model1, dataAipwe)
    u2 <- predict.glm(model2, dataAipwe)
    u3 <- predict.glm(model3, dataAipwe)
    u4 <- predict.glm(model3, dataAipwe)
    est.outcomes <- cbind(u1,u2,u3,u4)
    
    ##Propensity model
    PRaipwe <- table(A)/N
    
    #estimate contrast function
    c1 <- (Y*(A==1)/PRaipwe[1]) - (((A==1) - PRaipwe[1]) / PRaipwe[1]*u1)
    c2 <- (Y*(A==2)/PRaipwe[2]) - (((A==2) - PRaipwe[2]) / PRaipwe[2]*u2)
    c3 <- (Y*(A==3)/PRaipwe[3]) - (((A==3) - PRaipwe[3]) / PRaipwe[3]*u3)
    c4 <- (Y*(A==4)/PRaipwe[4]) - (((A==4) - PRaipwe[4]) / PRaipwe[4]*u4)
    est.aipwe <- cbind(c1,c2,c3,c4)
    
    #response and losses
    Zaipwe <- sapply(1:nrow(dataAipwe), function(b) which.max(est.aipwe[b,]))
    dataCaipwe <- cbind(as.factor(Zaipwe), dataAipwe[,3:7])
    
    #CART
    controlrpart <- rpart.control(maxdepth=6, maxsurrogate=0, maxcompete=0)
    
    indiv.loss=t(sapply(1:nrow(dataCaipwe), function(b) max(est.aipwe[b,])-est.aipwe[b,]))
    
    
    fit <-  rpart(as.numeric(Zaipwe) ~ X1 + X2 + X3 + X4 + X5, 
                  data = dataCaipwe, weights=1:N, control=controlrpart, method=multirisk,
                  parms=list(indiv.loss=t(sapply(1:nrow(dataCaipwe), function(b)   
                    max(est.aipwe[b,])-est.aipwe[b,])))
    )                                                    
    
    
    options(warn=-1)
    
      if(nrow(fit$frame)>1){
        xfit <- xpred.rpart(fit,xval=20)
        xerror <- colMeans(sapply(1:ncol(xfit), function(m) sapply(1:nrow(xfit),    
                                                                   function(k) indiv.loss[k,xfit[k,m]])))
        cp1 <- fit$cptable[which.min(xerror),1]
        fit <- prune(fit, cp=cp1)
      }
      
      if(!nrow(fit$frame)>1){
        Outcome <- max(aggregate(dataAipwe$Y, list(dataAipwe$A), mean)[,2])
        Ag <- which.max(aggregate(dataAipwe$Y, list(dataAipwe$A), mean)[,2])
        assignment <- rep(Ag,N)
        ExpOut <- mean(dataAipweE[,c("Y1","Y2","Y3","Y4")][,Ag])
        nsplits <-0
        splitvar <- FALSE
        splitpoint <- FALSE
        pred.Y <- rep(Outcome,N)
        tree<- NA 
        tree.created <- F
      }else{print(fit)} 
      
      if(nrow(fit$frame)>1){     # AIPWE schatter voor expected outcome
        noderow <- fit$where
        assignment <- fit$frame[noderow,5]
        indicator <- assignment==A
        m <- sapply(1:N, function(v) est.outcomes[v,assignment[v]])
        propensity <- PRaipwe[assignment]
        Outcome <- mean(((Y*indicator)/propensity) - ((indicator-propensity)/propensity*m))
        
        nsplits <- sum(fit$frame[,"var"]!="<leaf>")
        splitvar <-  rownames(fit$splits)
        splitpoint <- fit$splits[,"index"]
        tree <- list(fit$frame, fit$where, fit$splits)
        tree.created <- T
        
        #estimate expected outcome in population under regime
        Agexp <- predict(fit, dataAipweE)
        Y.true <- dataAipweE[,c("Y1","Y2","Y3","Y4")]
        ExpOut <- mean(sapply(1:length(Agexp), function(s) Y.true[s, Agexp[s]]))
      }
  }
  
  biassample <- Outcome - OutcomeOpt 
  cat("The expected sample outcome is", Outcome, "\n")
  cat("The expected population outcome is", ExpOut, "\n")
  
  return(list(Outcome=Outcome,ExpOut=ExpOut, tree.created=tree.created, nsplits=nsplits, splitvar=splitvar,splitpoint=splitpoint,
              assignment=assignment, opt.assignment=Topt, tree=tree, data.seed=dataAipwe, biassample=biassample)) 
}
