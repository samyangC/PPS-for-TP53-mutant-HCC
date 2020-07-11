library(survival)
load("testdata.Rdata")
dataEXP <- data[,3:ncol(data)]
dataEXP <- scale(dataEXP)
data <- cbind(data[,1:2],dataEXP)

coxR=data.frame()
coxf<-function(x){
fmla1 <- as.formula(Surv(OStime,OS)~data[,x])
mycox <- coxph(fmla1,data=data)
}

for(a in colnames(data[,3:ncol(data)])){
mycox=coxf(a)
coxResult = summary(mycox)
coxR=rbind(coxR,cbind(id=a,HR=coxResult$coefficients[,"exp(coef)"],
P=coxResult$coefficients[,"Pr(>|z|)"]))
}

sigGene <- coxR[as.numeric(as.character(coxR$P))<0.05,]
data1 <- cbind(data[,1:2],data[,as.character(sigGene$id)])
save(data1,file="testdata2.Rdata")

# choose replace = FALSE for sampling without replacement (SRSWOR)
library(survival)
load("testdata2.Rdata")

patients=rownames(data1)
outTab=data.frame()

for(gene in colnames(data1[,3:ncol(data1)])){
Mboot <- replicate(1000, expr = {
    indices <- sample(patients, size =nrow(data1)*0.7, replace = TRUE)
    data<-data1[indices,]
    fmla1 <- as.formula(Surv(data[,"OStime"],data[,"OS"])~data[,gene])
	mycox <- coxph(fmla1,data=data)
	coxResult = summary(mycox)
	P=coxResult$coefficients[,"Pr(>|z|)"]
})
times=length(Mboot[which(Mboot<0.05)])
outTab=rbind(outTab,cbind(gene=gene,times=times))
}

bootGene <- outTab[as.numeric(as.character(outTab$times))>900,]
data2 <- cbind(data1[,1:2],data1[,as.character(bootGene$gene)])
save(data2,file="testdata3.Rdata")

library("randomForestSRC")
library("Hmisc")
load("testdata3.Rdata")

# conatructing RSF model
res.rsf <- rfsrc(Surv(OStime, OS) ~ ., data2, nodesize = 20, proximity=T, tree.err = T, 
                        forest = T, ntree = 1000, splitrule = "logrank", importance = TRUE)

# MD-based virable selection
# method="vh" is used for problems where the number of variables is substantially larger than the sample size 
# (e.g., p/n is greater than 10). It is always prefered to use method="md", but to find more variables, 
# or when computations are high, variable hunting may be preferred.
 res.trc   <-c()
 res.trcoob<-c()
 res.testc <-c()
 topvars   <-vector(mode="character",length=25) 
    for (j in 1:1000) { 
     print(paste("trying for",j,"times"))
     vars<-var.select(object=res.rsf,
                     cause =1,
                     method = "md", 
                     conservative = c("high"), 
                     ntree = 1000,
                     nodesize = 20, splitrule = "logrank", nsplit = 3, xvar.wt = NULL,
                     refit = T, fast = T,
                     na.action = c("na.impute"), 
                     always.use = NULL, nrep = 10,
                     prefit =  list(action = T, ntree = 1000,
                                    nodesize = 20, nsplit = 3),
                     verbose = TRUE)
    
    # calculating C-index
    trc<-rcorr.cens(-vars$rfsrc.refit.obj$predicted, 
                    Surv(data2$OStime, data2$OS))["C Index"]
    trcoob<-rcorr.cens(-vars$rfsrc.refit.obj$predicted.oob, 
                    Surv(data2$OStime, data2$OS))["C Index"]
    
    # collecting C-index in every iteration
    res.trc    <-rbind(res.trc, trc)
    res.trcoob <-rbind(res.trcoob, trcoob)
    topvars    <-rbind(topvars,vars$topvars)
    }

# result intergrating
 result<-data.frame(res.trc,res.trcoob,row.names = 1:nrow(res.trc))
 colnames(result)<-c("res.trc.cindex","res.trcoob.cindex")
 topvars<-as.matrix(topvars)[-1,] 
 rownames(topvars)<-c(1:dim(topvars)[1])

# finding best model
 bestresult<-result[result$res.trcoob.cindex==max(result$res.trcoob.cindex),] 
 bestvars<-unique(topvars[rownames(bestresult),]) 

data3 <- cbind(data2[,1:2],data2[,as.character(bestvars)])
save(data3,file="testdata4.Rdata")
