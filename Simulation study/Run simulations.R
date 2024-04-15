source('E:/KU Leuven/kul/KUL - multivalued/Paper may 2023/code/Multivalued code snippet.R')
source('E:/KU Leuven/kul/KUL - multivalued/Paper may 2023/code/New simulations/Simulation code stat AIPWE 3 4 arms.R')
source('E:/KU Leuven/kul/KUL - multivalued/Paper may 2023/code/New simulations/Simulation code stat AIPWE 3 4 arms prob.R')
source('C:/Users/aniek/OneDrive/Desktop/Simulations multivalued/Simulation code RG 3 4 arms.R')
source('E:/KU Leuven/kul/KUL - multivalued/Paper may 2023/code/New simulations/Simulation code RF 3 4 arms.R')
source('C:/Users/aniek/OneDrive/Desktop/Simulations multivalued/Data generation RCT 3 4 arms.R')
##


##########################################################################
### Call data generation, regime estimation and summary of the results ###
##########################################################################
simulation.stat.AIPWE.emp <- function(arms, scenario, n.arm, meandif, model, iterations, a){
  dataAipweE <- datagen3.4(arms, scenario, 10^5/arms, meandif,model,h=1)[[1]]
  result <- sapply(1:iterations, function(h) C.aipwe(arms, scenario, n.arm, meandif,model,a,dataAipweE,h))
  save(result, file=paste("C.AIPWE arms=",arms,"_s=", scenario,"_n.arm=", n.arm,"_c=", meandif,"_m=",model,"_I=", iterations,"_cp=.01_",".RData",sep=""))
  write.csv2(descript.rpart(x=result, scenario=scenario, arms=arms), file=paste("C.AIPWE arms=",arms,"_s=", scenario,"_n.arm=", n.arm,"_c=", meandif,"_m=",model,"_I=", iterations,"_cp=.01,",".csv",sep=""), row.names=FALSE)
}  

simulation.stat.AIPWE.prob <- function(arms, scenario, n.arm, meandif, model, iterations, a){
  dataAipweE <- datagen3.4(arms, scenario, 10^5/arms, meandif,model,h=1)[[1]]
  result <- sapply(1:iterations, function(h) C.aipwe.prob(arms, scenario, n.arm, meandif,model,a,dataAipweE,h))
  save(result, file=paste("C.AIPWE.prob arms=",arms,"_s=", scenario,"_n.arm=", n.arm,"_c=", meandif,"_m=",model,"_I=", iterations,"_cp=.01",".RData",sep=""))
  write.csv2(descript.rpart(x=result, scenario=scenario, arms=arms), file=paste("C.AIPWE arms=",arms,"_s=", scenario,"_n.arm=", n.arm,"_c=", meandif,"_m=",model,"_I=", iterations,"_cp=.01",".csv",sep=""), row.names=FALSE)
}  

simulationRG <- function(arms, scenario, n.arm, meandif, model, iterations){
  dataAipweE <- datagen3.4(arms, scenario, 10^5/arms, meandif,model,h=1)[[1]]
  result <- sapply(1:iterations, function(h) RG(arms, scenario, n.arm, meandif,model,dataAipweE,h))
  save(result, file=paste("RG arms=",arms,"_s=", scenario,"_n.arm=", n.arm,"_c=", meandif,"_m=",model,"_I=", iterations,".RData",sep=""))
  write.csv2(descript.RG(x=result, scenario=scenario), file=paste("RG arms=",arms,"_s=", scenario,"_n.arm=", n.arm,"_c=", meandif,"_m=",model,"_I=", iterations,".csv",sep=""), row.names=FALSE)
}  

simulation.RF <- function(arms, scenario, n.arm, meandif, model, iterations, a){
  dataAipweE <- datagen3.4(arms, scenario, 10^5/arms, meandif,model,h=1)[[1]]
  result <- sapply(1:iterations, function(h) RF(arms, scenario, n.arm, meandif,model,a,dataAipweE,h))
  save(result, file=paste("RF arms=",arms,"_s=", scenario,"_n.arm=", n.arm,"_c=", meandif,"_m=",model,"_I=", iterations,"_cp=.01",".RData",sep=""))
  write.csv2(descript.rpart(x=result, scenario=scenario, arms=arms), file=paste("C.AIPWE arms=",arms,"_s=", scenario,"_n.arm=", n.arm,"_c=", meandif,"_m=",model,"_I=", iterations,"_cp=.01,",".csv",sep=""), row.names=FALSE)
}  




#######################  
### Run simulations ###
####################### 
setwd("E:/KU Leuven/kul/KUL - multivalued/Paper may 2023/code/New simulations")

for(arms in c(3,4)){
  for(scenario in c(1,2,3)){
    for(n.arm in c(50,100,200)){ 
      for(meandif in c(1,.5)){
          simulation.stat.AIPWE.emp(arms, scenario, n.arm, meandif, model=1, iterations=500, a=20)
        }
      }
    }
  }    



for(arms in c(3)){
  for(scenario in c(2)){
    for(n.arm in c(100)){ 
      for(meandif in c(1)){
          simulation.stat.AIPWE.prob(arms, scenario, n.arm, meandif, model=1, iterations=500, a=20)
      }
    }
  }    
}

for(arms in c(3,4)){
  for(scenario in c(1,2,3)){
    for(n.arm in c(50,100,200)){ 
      for(meandif in c(1,.5)){
        if(!file.exists(paste("RG arms=",arms,"_s=", scenario,"_n.arm=", n.arm,"_c=", meandif,"_m=1","_I=500",".RData",sep=""))){
          simulationRG(arms, scenario, n.arm, meandif,model=1,iterations=500) 
        }
      }
    }
  }    
}

for(arms in c(3,4)){
  for(scenario in c(1,2,3)){
    for(n.arm in c(50,100,200)){ 
      for(meandif in c(1,.5)){
        simulation.RF(arms, scenario, n.arm, meandif, model=1, iterations=500, a=20)
      }
    }
  }    
}


