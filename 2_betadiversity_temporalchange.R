#temporal changes

########################################
wdm="/data"
prefixm="W16SArchaea"
code.wd="/SourceCodes"

com.file="resample_otu_table.txt"
treat.file="treatment info and env data.txt"
beta.file=list()
beta.file[[1]]="betadiversity_Sorensen.txt"
beta.file[[2]]="betadiversity_unweightedUnifrac.txt"
names(beta.file)=c("Sorensen","Unifrac.uw")

#Bacteria files (optional) #for comparison
betaBAC.file=list()
betaBAC.file[[1]]="bacteria_betadiversity_Sorensen.txt"
betaBAC.file[[2]]="bacteria_betadiversity_unweightedUnifrac.txt"
names(betaBAC.file)=c("Sorensen","Unifrac.uw")

##########
library("vegan")
library("ieggr")
code.wd=iwd(code.wd)
wdm=iwd(wdm)
treat=lazyopen(treat.file)

head(treat)
warming=as.factor(treat$Warm)
block=as.factor(treat$block)
year=as.factor(treat$year)

######################
wdc<-function(wdi,...)
{
  wdci=paste0(wdm,"/",wdi)
  if(!dir.exists(wdci)){dir.create(wdci)}
  iwd(wdci)
}

idc=lazyopen(idc.file)

######################
beta=list()
for(i in 1:length(beta.file))
{
  betai=lazyopen(beta.file[[i]])
  print(dim(betai))
  print(betai[1:5,1:5])
  sampc=match.name(name.check = rownames(treat),both.list = list(betai=betai))
  betai=sampc$betai
  sampc=match.name(name.check = rownames(treat),both.list = list(betai=betai))
  print(betai[1:5,1:5])
  beta[[i]]=betai
}

names(beta)=names(beta.file)

#Bacteria files (optional) #for comparison 
beta.BAC=list()
for(i in 1:length(betaBAC.file))
{
  betaBACi=lazyopen(betaBAC.file[[i]])
  print(dim(betaBACi))
  print(betaBACi[1:5,1:5])
  sampc=match.name(name.check = rownames(treat),both.list = list(betaBACi=betaBACi))
  betaBACi=sampc$betaBACi
  sampc=match.name(name.check = rownames(treat),both.list = list(betaBACi=betaBACi))
  print(betaBACi[1:5,1:5])
  beta.BAC[[i]]=betaBACi
}

names(beta.BAC)=names(betaBAC.file)

#############################
# 1 # dissimilarity test
ds.wd=wdc("Dissimilarity")

for(i in 1:length(beta))
{
  message("Now i=",i,". ",date())
  prefixi=paste0(prefixm,".",names(beta)[i])
  b.dis=as.dist(beta[[i]])
  
  # 1.1 # adonis
  ad.blk=adonis(b.dis~warming+year*block,permutations = 999)
  ad.nblk=adonis(b.dis~warming+year,permutations = 999)
  
  sink(file=paste0(prefixi,".Adonis.txt"))
  print("---- With Block ----")
  print(ad.blk)
  print("---- Without Block ----")
  print(ad.nblk)
  sink()
  
  # 1.2 # ANOSIM
  ano.blk=anosim(b.dis, grouping = warming,permutations = 999,
                 strata = as.factor(paste(block,year,sep = ".")))
  ano.nblk=anosim(b.dis, grouping = warming,permutations = 999,
                  strata = as.factor(year))
  sink(file=paste0(prefixi,".ANOSIM.txt"))
  print("---- With Block ----")
  print(ano.blk)
  print("---- Without Block ----")
  print(ano.nblk)
  sink()
  
  # 1.3 # MRPP
  mrp.blk=mrpp(dat = b.dis, grouping = warming,permutations = 999,
               strata = as.factor(paste(block,year,sep = ".")))
  mrp.nblk=mrpp(dat = b.dis, grouping = warming,permutations = 999,
                strata = as.factor(year))
  mrp.blk
  mrp.nblk
  sink(file=paste0(prefixi,".MRPP.txt"))
  print("---- With Block ----")
  print(mrp.blk)
  print("---- Without Block ----")
  print(mrp.nblk)
  sink()
}

#########################
# 2 # temporal changes
dv.wd=wdc("Convergent")
library("lme4")
library("car")
source(paste0(code.wd,"/rsquaredglmm.r"))
source(paste0(code.wd,"/dvt.r"))

dv=list()
for(i in 1:length(beta))
{
  message("-----Now i=",i,". ",date())
  prefixi=paste0(prefixm,".",names(beta)[i])
  betai=beta[[i]]
  betaBACi=beta.BAC[[i]]
  sampc=match.name(both.list = list(betai,betaBACi))
  dv[[i]]=dvt(betai=betai,betaBi=betaBACi,treat=treat,boot=TRUE,rand=1000,prefixi=paste0("Aarchaea.Bbacteria.",names(beta)[i]))
  if(i==1)
  {
    write.csv(dv[[i]]$summary,file = paste0(prefixi,".Convergent.Summary.csv"),append = FALSE)
  }else{
    write.csv(dv[[i]]$summary,file = paste0(prefixi,".Convergent.Summary.csv"),append = TRUE)
  }
  save(dv,file = paste0(prefixi,".Convergent.detail.rda"))
}



