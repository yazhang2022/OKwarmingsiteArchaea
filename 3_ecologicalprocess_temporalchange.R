#ecological processes temporal change
########################################

library("vegan")
library("ieggr")
library("iCAMP")
code.wd="/SourceCodes"

setwd("/data")
treat.file="treatment info and env data.txt"
treat=lazyopen(treat.file)
head(treat)

datarch=read.table("ecological processes by year.txt", sep="\t", header= T)
head(datarch)
datarch$sto=1-datarch$Heterogeneous.Selection-datarch$Homogeneous.Selection
prefixi="iCAMP_archaea16S.ProcessImportance_"

alpha=list()
for(k in 3:length(datarch)){
  alpha[[k]]=datarch[,c(1,2,k)]
  write.csv(as.matrix(alpha[[k]]),file = paste0(prefixi,colnames(datarch)[k],".csv"))
}

#Bacteria files (optional) #for comparison
datbact=read.table("bacteria ecological processes by year.txt", sep="\t", header= T)
head(datbact)
datbact$sto=1-datbact$Heterogeneous.Selection-datbact$Homogeneous.Selection
prefixi="iCAMP_bacteria16S.ProcessImportance_"

delta=list()
for(k in 3:length(datbact)){
  delta[[k]]=datbact[,c(1,2,k)]
  write.csv(as.matrix(delta[[k]]),file = paste0(prefixi,colnames(datbact)[k],".csv"))
}

beta.file=list()
beta.file[[1]]="iCAMP_archaea16S.ProcessImportance_HeS.csv"
beta.file[[2]]="iCAMP_archaea16S.ProcessImportance_HoS.csv"
beta.file[[3]]="iCAMP_archaea16S.ProcessImportance_DL.csv"
beta.file[[4]]="iCAMP_archaea16S.ProcessImportance_HD.csv"
beta.file[[5]]="iCAMP_archaea16S.ProcessImportance_DR.csv"
beta.file[[6]]="iCAMP_archaea16S.ProcessImportance_Sto.csv"
names(beta.file)=c("HeS","HoS","DL","HD","DR","Sto")

beta=list()
for(i in 1:length(beta.file))
{
  betai=lazyopen(beta.file[[i]])
  betai=col3.dist(betai,to.dist=FALSE)
  sampc=match.name(name.check = rownames(treat),both.list = list(betai=betai))
  betai=sampc$betai
  sampc=match.name(name.check = rownames(treat),both.list = list(betai=betai))
  print(betai[1:5,1:5])
  beta[[i]]=betai
}
names(beta)=names(beta.file)

betaBAC.file=list()
betaBAC.file[[1]]="iCAMP_bacteria16S.ProcessImportance_HeS.csv"
betaBAC.file[[2]]="iCAMP_bacteria16S.ProcessImportance_HoS.csv"
betaBAC.file[[3]]="iCAMP_bacteria16S.ProcessImportance_DL.csv"
betaBAC.file[[4]]="iCAMP_bacteria16S.ProcessImportance_HD.csv"
betaBAC.file[[5]]="iCAMP_bacteria16S.ProcessImportance_DR.csv"
betaBAC.file[[6]]="iCAMP_bacteria16S.ProcessImportance_Sto.csv"
names(betaBAC.file)=c("HeS","HoS","DL","HD","DR","Sto")

betaBAC=list()
for(i in 1:length(betaBAC.file))
{
  betaBACi=lazyopen(betaBAC.file[[i]])
  betaBACi=col3.dist(betaBACi,to.dist=FALSE)
  sampc=match.name(name.check = rownames(treat),both.list = list(betaBACi=betaBACi))
  betaBACi=sampc$betaBACi
  sampc=match.name(name.check = rownames(treat),both.list = list(betaBACi=betaBACi))
  print(betaBACi[1:5,1:5])
  betaBAC[[i]]=betaBACi
}
names(betaBAC)=names(betaBAC.file)

#########################
# 1 # temporal changes
library("lme4")
library("car")
source(paste0(code.wd,"/rsquaredglmm.r"))
source(paste0(code.wd,"/dvt.r"))

prefixm="16SArchaea_sto"

dv=list()
for(i in 1:length(beta))
{
  tryCatch({
  message("-----Now i=",i,". ",date())
  prefixi=paste0(prefixm,".",names(beta)[i])
  betai=beta[[i]]
  betaBACi=betaBAC[[i]] #bacterial optional
  sampc=match.name(both.list = list(betai,betaBACi))
  dv[[i]]=dvt(betai=betai,betaBi=betaBACi,treat=treat,boot=TRUE,rand=1000,prefixi)
  if(i==1)
  {
    write.csv(dv[[i]]$summary,file = paste0(prefixi,".Summary.csv"),append = FALSE)
  }else{
    write.csv(dv[[i]]$summary,file = paste0(prefixi,".Summary.csv"),append = TRUE)
  }
  save(dv,file = paste0(prefixi,".detail.rda"))
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

