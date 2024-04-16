
datagen3.4 <- function(arms, scenario, n.arm, meandif,model, h){
  set.seed(h)
  n <- arms*n.arm
  covar <- as.data.frame(cbind(round(matrix(rnorm(n*5),ncol=5),2),as.factor(sample(letters[1:4], n, replace = T))))
  covar[,6] <- as.factor(covar[,6])
  error <- rnorm(n,mean=0,sd=1)
  
  if(arms==3){
    gopt <- vector()
  A <- sample(1:3,n, replace=T)
  if(scenario==1){
    for(i in 1:n){
      if (covar[i,1]<= -.545) gopt[i] <- 1
      if ((covar[i,1]> -.545)*(covar[i,2]<=0)) gopt[i] <- 2
      if ((covar[i,1]> -.545)*(covar[i,2]> 0)) gopt[i] <- 3
    }
  }
  if(scenario==2){
    for(i in 1:n){
      if ((covar[i,1]>covar[i,2])*(covar[i,1]>covar[i,3])) {gopt[i] <- 1} else {
        if ((covar[i,1]<=covar[i,2])*(covar[i,1]<=covar[i,3])) {gopt[i] <- 2} else {
          gopt[i] <- 3 }}
    }
  }
  if(scenario==3){gopt <- rep(1,n)}
  
 
  if(model==1 & scenario!=3){ Y <- 1+(.25*covar[,1])+(.25*covar[,2])-(.25*covar[,5])-(meandif*as.numeric(A!=gopt)) + error}
  if(model==1 & scenario!=3){ Y1 <- 1+(.25*covar[,1])+(.25*covar[,2])-(.25*covar[,5])-(meandif*as.numeric(1!=gopt))}
  if(model==1 & scenario!=3){ Y2 <- 1+(.25*covar[,1])+(.25*covar[,2])-(.25*covar[,5])-(meandif*as.numeric(2!=gopt))}
  if(model==1 & scenario!=3){ Y3 <- 1+(.25*covar[,1])+(.25*covar[,2])-(.25*covar[,5])-(meandif*as.numeric(3!=gopt))}
  if(model==1 & scenario==3){ Y <- 1+(.25*covar[,1])+(.25*covar[,2])-(.25*covar[,5])-(meandif*as.numeric(A!=gopt)) + error}
  if(model==1 & scenario==3){ Y1 <- 1+(.25*covar[,1])+(.25*covar[,2])-(.25*covar[,5])}
  if(model==1 & scenario==3){ Y2 <- 1+(.25*covar[,1])+(.25*covar[,2])-(.25*covar[,5])-meandif}
  if(model==1 & scenario==3){ Y3 <- Y2}
  
  dataG <-  cbind(round(Y,2),as.factor(A),covar,Y1,Y2,Y3)
  colnames(dataG) <- c("Y","A","X1","X2","X3","X4","X5","X6","Y1","Y2","Y3")
  }
  
  if(arms==4){
    gopt <- vector()
    A <- sample(1:4,n, replace=T)
    if(scenario==1){
      for(i in 1:n){
        if (covar[i,1]<= -.745) gopt[i] <- 1
        if ((covar[i,1]> -.745)*(covar[i,2]<= -.545)) gopt[i] <- 2
        if ((covar[i,1]> -.745)*(covar[i,2] > -.545)*(covar[i,3] <= 0)) gopt[i] <- 3
        if ((covar[i,1]> -.745)*(covar[i,2] > -.545)*(covar[i,3] > 0)) gopt[i] <- 4
      }
    }
    if(scenario==2){
      for(i in 1:n){
        if ((covar[i,1]> covar[i,2])*(covar[i,3]> covar[i,4])) {gopt[i] <- 1} 
        if ((covar[i,1]<=covar[i,2])*(covar[i,3]<=covar[i,4])) {gopt[i] <- 2} 
        if ((covar[i,1]> covar[i,2])*(covar[i,3]<=covar[i,4])) {gopt[i] <- 3} 
        if ((covar[i,1]<=covar[i,2])*(covar[i,3]> covar[i,4])) {gopt[i] <- 4} 
      }
    }
    if(scenario==3){gopt <- rep(1,n)}
    
    if(scenario==4){
      for(i in 1:n){
        if (((covar[i,1]^2)> covar[i,2])*((covar[i,3]^2)> covar[i,4])) {gopt[i] <- 1} 
        if (((covar[i,1]^2)<=covar[i,2])*((covar[i,3]^2)<=covar[i,4])) {gopt[i] <- 2} 
        if (((covar[i,1]^2)> covar[i,2])*((covar[i,3]^2)<=covar[i,4])) {gopt[i] <- 3} 
        if (((covar[i,1]^2)<=covar[i,2])*((covar[i,3]^2)> covar[i,4])) {gopt[i] <- 4} 
      }
    }
    
    if(scenario==5){
      for(i in 1:n){
        if (((covar[i,1]^2)> covar[i,2])*(covar[i,3]> covar[i,4])) {gopt[i] <- 1} 
        if (((covar[i,1]^2)<=covar[i,2])*(covar[i,3]<=covar[i,4])) {gopt[i] <- 2} 
        if (((covar[i,1]^2)> covar[i,2])*(covar[i,3]<=covar[i,4])) {gopt[i] <- 3} 
        if (((covar[i,1]^2)<=covar[i,2])*(covar[i,3]> covar[i,4])) {gopt[i] <- 4} 
      }
    }
    
    if(model==1 & scenario!=3){ Y <- 1+(.25*covar[,1])+(.25*covar[,2])-(.25*covar[,5])-(meandif*as.numeric(A!=gopt)) + error}
    if(model==1 & scenario!=3){ Y1 <- 1+(.25*covar[,1])+(.25*covar[,2])-(.25*covar[,5])-(meandif*as.numeric(1!=gopt))}
    if(model==1 & scenario!=3){ Y2 <- 1+(.25*covar[,1])+(.25*covar[,2])-(.25*covar[,5])-(meandif*as.numeric(2!=gopt))}
    if(model==1 & scenario!=3){ Y3 <- 1+(.25*covar[,1])+(.25*covar[,2])-(.25*covar[,5])-(meandif*as.numeric(3!=gopt))}
    if(model==1 & scenario!=3){ Y4 <- 1+(.25*covar[,1])+(.25*covar[,2])-(.25*covar[,5])-(meandif*as.numeric(4!=gopt))}
    if(model==1 & scenario==3){ Y <- 1+(.25*covar[,1])+(.25*covar[,2])-(.25*covar[,5])-(meandif*as.numeric(A!=gopt)) + error}
    if(model==1 & scenario==3){ Y1 <- 1+(.25*covar[,1])+(.25*covar[,2])-(.25*covar[,5])}
    if(model==1 & scenario==3){ Y2 <- 1+(.25*covar[,1])+(.25*covar[,2])-(.25*covar[,5])-meandif}
    if(model==1 & scenario==3){ Y3 <- Y2}
    if(model==1 & scenario==3){ Y4 <- Y2}
    
    dataG <-  cbind(round(Y,2),as.factor(A),covar,Y1,Y2,Y3,Y4)
    colnames(dataG) <- c("Y","A","X1","X2","X3","X4","X5","X6","Y1","Y2","Y3","Y4")
  }
  
  return(list(as.data.frame(dataG), gopt))
}

