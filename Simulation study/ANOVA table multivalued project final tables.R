library(irr)

loaddataRG <- function(arms, n.arm, scenario, meandif){
  load(file=paste("RG arms=",arms,"_s=", scenario,"_n.arm=", n.arm,"_c=", meandif,"_m=1_I=500.RData",sep=""))
  return(result)
}
#loaddataRG(3,50,1,1)




loaddataCaipwe <- function(arms, n.arm, scenario, meandif,crit){
  load(file=paste("C.AIPWE arms=",arms,"_s=", scenario,"_n.arm=", n.arm,"_c=", meandif,"_m=1_I=500_cp=.01_",".RData",sep=""))
  return(result)
}
#loaddataCaipwe(3,50,1,1,"error")


loaddataCaipwe.prob <- function(arms, n.arm, scenario, meandif,crit){
  load(file=paste("C.AIPWE.prob arms=",arms,"_s=", scenario,"_n.arm=", n.arm,"_c=", meandif,"_m=1_I=500_cp=.01",".RData",sep=""))
  return(result)
}

loaddataRF <- function(arms, n.arm, scenario, meandif,crit){
  load(file=paste("RF arms=",arms,"_s=", scenario,"_n.arm=", n.arm,"_c=", meandif,"_m=1_I=500_cp=.01",".RData",sep=""))
  return(result)
}

crit<- "error"

makedata<- function(){
  AnovaData<- matrix(c(0), 18000, 29) #2*3*3*2=36 *500= 18000 3*3*2=18*500=9000
  colnames(AnovaData)<- c("dataset","arms", "n.arms", "scenario", "c", "1egopt", "1Egopt", "1type1", "1type2", "1assign", "1gain", 
                                                          "2egopt", "2Egopt", "2type1", "2type2", "2assign", "2gain",
                                                          "3egopt", "3Egopt", "3type1", "3type2", "3assign", "3gain",
                                                          "4egopt", "4Egopt", "4type1", "4type2", "4assign", "4gain"
                          )
  #1=rg, 2=Caipwe, 3=Caipwe prob, 4=RF   
  
  i<- 1
  setwd("C:/Users/aniek/OneDrive/Desktop/Simulations multivalued/run 2")  
  for(arms in c(3,4)){                                                                  ### 4
     for(n.arm in c(50, 100, 200)){ 
      for(scenario in c(1,2,3)){
        for(meandif in c(1, 0.5)){
          
          x <- loaddataRG(arms, n.arm, scenario, meandif)
          xx <- loaddataCaipwe(arms, n.arm, scenario, meandif, crit)
          xxx <- loaddataCaipwe.prob(arms, n.arm, scenario, meandif, crit)
          xxxx <- loaddataRF(arms, n.arm, scenario, meandif, crit)
          
                  AnovaData[seq(from=i,to=i+499),1]<- seq(from=i, to=i+499)         #Replication number
                  AnovaData[seq(from=i,to=i+499),2]<- rep(arms, 500)          
                  AnovaData[seq(from=i,to=i+499),3]<- rep(n.arm,500)    
                  AnovaData[seq(from=i,to=i+499),4]<- rep(scenario, 500)          
                  AnovaData[seq(from=i,to=i+499),5]<- rep(meandif, 500)          
                  
                  
                  AnovaData[seq(from=i,to=i+499),"1egopt"] <- as.numeric(x[1,])     #e^gopt
                  AnovaData[seq(from=i,to=i+499),"2egopt"] <- as.numeric(xx[1,])    
                  AnovaData[seq(from=i,to=i+499),"3egopt"] <- as.numeric(xxx[1,])    
                  AnovaData[seq(from=i,to=i+499),"4egopt"] <- as.numeric(xxxx[1,])     
                  
                  AnovaData[seq(from=i,to=i+499),"1Egopt"] <- as.numeric(x[2,])     #E^gopt: 
                  AnovaData[seq(from=i,to=i+499),"2Egopt"] <- as.numeric(xx[2,])           
                  AnovaData[seq(from=i,to=i+499),"3Egopt"] <- as.numeric(xxx[2,])     
                  AnovaData[seq(from=i,to=i+499),"4Egopt"] <- as.numeric(xxxx[2,])       
                
                  if (scenario==3){                                                 #Type 1 error
                    AnovaData[seq(from=i,to=i+499),"1type1"]<- as.numeric(x[3,])       
                    AnovaData[seq(from=i,to=i+499),"2type1"]<- as.numeric(xx[3,]) 
                    AnovaData[seq(from=i,to=i+499),"3type1"]<- as.numeric(xxx[3,])       
                    AnovaData[seq(from=i,to=i+499),"4type1"]<- as.numeric(xxxx[3,]) 
                  }
                  if (!scenario==3){                                           
                    AnovaData[seq(from=i,to=i+499),"1type1"]<- NA
                    AnovaData[seq(from=i,to=i+499),"2type1"]<- NA
                    AnovaData[seq(from=i,to=i+499),"3type1"]<- NA
                    AnovaData[seq(from=i,to=i+499),"4type1"]<- NA
                  }
                  
                  if(!scenario==3) {                                                  #Type 2 error
                    AnovaData[seq(from=i,to=i+499),"1type2"] <- 1-as.numeric(x[3,])
                    AnovaData[seq(from=i,to=i+499),"2type2"] <- 1-as.numeric(xx[3,])
                    AnovaData[seq(from=i,to=i+499),"3type2"] <- 1-as.numeric(xxx[3,])
                    AnovaData[seq(from=i,to=i+499),"4type2"] <- 1-as.numeric(xxxx[3,])
                  }
                  if(scenario==3) {                              
                    AnovaData[seq(from=i,to=i+499),"1type2"] <- NA
                    AnovaData[seq(from=i,to=i+499),"2type2"] <- NA
                    AnovaData[seq(from=i,to=i+499),"3type2"] <- NA
                    AnovaData[seq(from=i,to=i+499),"4type2"] <- NA
                  }
                 
                  # Winst behaald berekenen
                  if(arms==3){Emargbest <- matrix(c(0.35, 0.68, 0.34, 0.67, 1,1,0.61, 0.81,0.45, .72), nrow=2)}
                  if(arms==4){Emargbest <- matrix(c(0.27, 0.64, 0.25, 0.63, 1,1, 0.52, 0.76, 0.36, 0.68), nrow=2)}
                  rownames(Emargbest) <- c("1", "0.5")
                  colnames(Emargbest) <- c("1", "2", "3","4","5")
                  Egopt<- rep(1.00,500)
                  Eg1<- as.numeric(x[2,])
                  Eg2<- as.numeric(xx[2,])
                  Eg3<- as.numeric(xxx[2,])
                  Eg4<- as.numeric(xxxx[2,])
                  
                  if(scenario!=3){
                  AnovaData[seq(from=i,to=i+499),"1gain"]<- (Eg1-Emargbest[paste(meandif),paste(scenario)])/(Egopt-Emargbest[paste(meandif),paste(scenario)])
                  AnovaData[seq(from=i,to=i+499),"2gain"]<- (Eg2-Emargbest[paste(meandif),paste(scenario)])/(Egopt-Emargbest[paste(meandif),paste(scenario)])
                  AnovaData[seq(from=i,to=i+499),"3gain"]<- (Eg3-Emargbest[paste(meandif),paste(scenario)])/(Egopt-Emargbest[paste(meandif),paste(scenario)])
                  AnovaData[seq(from=i,to=i+499),"4gain"]<- (Eg4-Emargbest[paste(meandif),paste(scenario)])/(Egopt-Emargbest[paste(meandif),paste(scenario)])
                  }
                  
                  if(scenario==3){
                  AnovaData[seq(from=i,to=i+499),"1gain"]<- NA
                  AnovaData[seq(from=i,to=i+499),"2gain"]<- NA
                  AnovaData[seq(from=i,to=i+499),"3gain"]<- NA
                  AnovaData[seq(from=i,to=i+499),"4gain"]<- NA
                  }
          
                
                  #incorrect assignment
                  if(scenario!=3){
                  Topt <-  matrix(unlist(x["opt.assignment",]),ncol=500)
                  
                  Tpred <- matrix(unlist(x["assignment",]),ncol=500)
                  print(nrow(Tpred))
                  correctassigned <- (Tpred==Topt) #TRUE if correctly assigned
                  AnovaData[seq(from=i,to=i+499),"1assign"] <- 1-apply(correctassigned, 2, mean)  #misclassification rate
                  
                  Tpred <- matrix(unlist(xx["assignment",]),ncol=500)
                  print(nrow(Tpred))
                   correctassigned <- (Tpred==Topt)
                  AnovaData[seq(from=i,to=i+499),"2assign"] <- 1-apply(correctassigned, 2, mean)  
                  
                  Tpred <- matrix(unlist(xxx["assignment",]),ncol=500)
                  print(nrow(Tpred))
                  correctassigned <- (Tpred==Topt)
                  AnovaData[seq(from=i,to=i+499),"3assign"] <- 1-apply(correctassigned, 2, mean)  
                  
                  Tpred <- matrix(unlist(xxxx["assignment",]),ncol=500)
                  print(nrow(Tpred))
                  correctassigned <- (Tpred==Topt)
                  AnovaData[seq(from=i,to=i+499),"4assign"] <- 1-apply(correctassigned, 2, mean)  

                  }
                  if(scenario==3){
                  Tpred <- matrix(unlist(x["assignment",]),ncol=500)
                  correctassigned <- (Tpred==1) #TRUE if correctly assigned
                  AnovaData[seq(from=i,to=i+499),"1assign"] <- 1-apply(correctassigned, 2, mean) 
                  
                  
                  Tpred <- matrix(unlist(xx["assignment",]),ncol=500)
                  correctassigned <- (Tpred==1) #TRUE if correctly assigned
                  AnovaData[seq(from=i,to=i+499),"2assign"] <- 1-apply(correctassigned, 2, mean)  
                  
                  Tpred <- matrix(unlist(xxx["assignment",]),ncol=500)
                  correctassigned <- (Tpred==1) #TRUE if correctly assigned
                  AnovaData[seq(from=i,to=i+499),"3assign"] <- 1-apply(correctassigned, 2, mean)  
                  
                  Tpred <- matrix(unlist(xxxx["assignment",]),ncol=500)
                  correctassigned <- (Tpred==1) #TRUE if correctly assigned
                  AnovaData[seq(from=i,to=i+499),"4assign"] <- 1-apply(correctassigned, 2, mean)  
                  
                  }
          
          i=i+500

                }
             
            }    
          
        }
  }
    
  save(AnovaData, file="AnovaData.Rdata")
}

makedata()
load("AnovaData.Rdata")

library(stats)
library(reshape2)
summary(AnovaData)
head(AnovaData)
colnames(AnovaData)
str(AnovaData)


### Type II error
dataType2 <- aggregate(AnovaData[,c(9, 15, 21, 27)] ~ arms + n.arms + scenario + c, data=AnovaData, FUN=mean)

write.table(dataType2, file="DataType2.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

### Type I error
dataType1 <- aggregate(AnovaData[,c(8, 14, 20, 26)] ~ arms + n.arms + c, data=AnovaData, FUN=mean)
write.table(dataType1, file="DataType1.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

## Analysis
#NPG
library(rstatix)

NPGdata <-  data.frame(AnovaData[,c(1:5,17,23,29)])
NPGdata$dataset <- factor(NPGdata$dataset)

data_long <- gather(NPGdata, method, value, X2gain:X4gain, factor_key=TRUE)
data_long$n.arms factor(data_long$n.arms)
data_long$arms <-  factor(data_long$arms)
data_long$scenario <-  factor(data_long$scenario)
data_long$c  <-  factor(data_long$c)

res.aov <- anova_test(
  data  <-  data_long, dv = value, wid = dataset,
  within = method, between = c(arms, n.arms, scenario,c), detailed=F,
  effect.size="pes"
)

anova_res  <-  get_anova_table(res.aov)
write.csv(anova_res, "NPG_anova.csv")

data_long <-  data_long[!is.na(data_long$value),]
sum((data_long$value-mean(data_long$value, na.rm=T))^2)

### table 
data_long <-  data_long[data_long$scenario!=3,]
tapply(data_long$value, data_long$arms, mean)
tapply(data_long$value, data_long$n.arms, mean)
tapply(data_long$value, data_long$c, mean)
tapply(data_long$value, data_long$scenario, mean)
tapply(data_long$value, data_long$method, mean)



# assignment without benchmark
assignmentdata  <-  data.frame(AnovaData[,c(1:5,16,22,28)])
assignmentdata$dataset <- factor(assignmentdata$dataset)

data_long <- gather(assignmentdata, method, value, X2assign:X4assign, factor_key=TRUE)
data_long$n.arms <- factor(data_long$n.arms)
data_long$arms <- factor(data_long$arms)
data_long$scenario <- factor(data_long$scenario)
data_long$value  <-  1-data_long$value
data_long$c  <-  factor(data_long$c)

res.aov <- anova_test(
  data = data_long, dv = value, wid = dataset,
  within = method, between = c(arms, n.arms, scenario,c), detailed=TRUE,
  effect.size="pes"
)

anova_res  <-  get_anova_table(res.aov)
write.csv(anova_res, "assignanova.csv")

data_long <-  data_long[!is.na(data_long$value),]
sum((data_long$value-mean(data_long$value, na.rm=T))^2)


tapply(data_long$value, data_long$arms, mean)
tapply(data_long$value, data_long$n.arms, mean)
tapply(data_long$value, data_long$c, mean)
tapply(data_long$value, data_long$scenario, mean)
tapply(data_long$value, data_long$method, mean)



## Type II without benchmark
dataType2$condition  <-  c(1:24)
data_long <- gather(dataType2, method, value, `2type2`:`4type2`, factor_key=TRUE)
data_long$n.arms <- factor(data_long$n.arms)
data_long$arms <- factor(data_long$arms)
data_long$c <-  factor(data_long$c)



res.aov <- anova_test(
  data = data_long, 
  formula = value ~ 1+ method+arms+n.arms+ c+ scenario+
    method:scenario+ method:arms+method:n.arms+ method:c+ 
    arms:n.arms+ arms:c+ arms:scenario+ n.arms:c+ n.arms:scenario+ c:scenario+
    method:arms:n.arms+ method:arms:c+ method:n.arms:c+ method: arms:scenario+ method: n.arms:scenario+
    method:c:scenario + arms:n.arms:scenario+ arms:c:scenario + n.arms:c:scenario+ 
    n.arms:arms:c +
    n.arms:arms:c:method+ n.arms:arms:c:scenario+ arms:c:scenario:method+
    n.arms:c:scenario:method+arms:n.arms:scenario:method,
  detailed=TRUE,
  effect.size="pes"
)


anova_res  <-  get_anova_table(res.aov)
write.csv(anova_res, "type2.csv")

sum((data_long$value-mean(data_long$value))^2)


tapply(data_long$value, data_long$arms, mean)
tapply(data_long$value, data_long$n.arms, mean)
tapply(data_long$value, data_long$c, mean)
tapply(data_long$value, data_long$scenario, mean)
tapply(data_long$value, data_long$method, mean)


### type 1
data_long <- gather(dataType1, method, value, `2type1`:`4type1`, factor_key=TRUE)
data_long$n.arms <- factor(data_long$n.arms)
data_long$arms <- factor(data_long$arms)
data_long$c <-  factor(data_long$c)
data_long$condition <- c(1:12)


res.aov <- anova_test(
  data = data_long, 
  formula = value ~method+arms+n.arms+ c+  method:arms+method:n.arms+ method:c+ arms:n.arms+ arms:c+n.arms:c+ method:arms:n.arms+ method:arms:c+ method:n.arms:c+ n.arms:arms:c,
  detailed=TRUE,
  effect.size="pes"
)

anova_res  <-  get_anova_table(res.aov)
write.csv(anova_res, "type1.csv")


tapply(data_long$value, data_long$arms, mean)
tapply(data_long$value, data_long$n.arms, mean)
tapply(data_long$value, data_long$c, mean)
tapply(data_long$value, data_long$scenario, mean)
tapply(data_long$value, data_long$method, mean)

data_long$method <- as.character(data_long$method)
data_long$method[data_long$method =="2type1"] <-  "AIPWE Empirical"
data_long$method[data_long$method=="3type1"] <-  "AIPWE Probabilities"
data_long$method[data_long$method=="4type1"] <- "Random Forest"

library(ggplot2)
type1_plot <-   ggplot(data=data_long, aes(x=n.arms, y=value, group=method)) +
  stat_summary(fun.y=mean, aes(group=method, linetype=method), geom="line") +
  stat_summary(fun.y=mean, aes(group=method, shape=method), geom="point") + 
  ylab("Type I error")+
  xlab("Sample size per arm")+
  coord_cartesian(ylim=c(0,1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

type1_plot
