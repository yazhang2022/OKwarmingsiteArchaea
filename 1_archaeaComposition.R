####operating on OTU tables
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")

library("phyloseq")
library("ggplot2")      # graphics
library("dplyr")        # filter and reformat data frames
library("ieggr")

#########################################
# 2 # set working files
#########################################
setwd("/data")
otu.file="resample_otu_table.txt"
treat.file="treatment info and env data.txt"
class.file="classifier_archaea.txt"

comm=t(lazyopen(otu.file))
dim(comm)
comm[1:5,1:5]
treat=lazyopen(treat.file)
head(treat)
treat=treat[,4, drop=FALSE]
head(treat)

comm<-comm[match(row.names(treat),row.names(comm)),]
comm[1:5,1:5]
dim(comm)
clas=lazyopen(class.file)
clas[1:5,1:5]

spcd=match.name(cn.list=list(comm=comm),rn.list=list(clas=clas))
comm1=spcd$comm
clas1=spcd$clas

dim(comm1)
dim(clas1)

clas2=as.matrix(clas1)
clas2[1:5,1:5]

#########################################
# 3 # main scripts
#########################################
set.seed(134)
OTU = otu_table(comm1, taxa_are_rows = FALSE)
sampled = sample_data(treat, errorIfNULL = FALSE)
TAX = tax_table(clas2)
physeq=phyloseq(OTU,TAX,sampled)
sample_names(physeq)
rank_names(physeq)
sample_variables(physeq)

colr=c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue",
       "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")

p1=plot_bar(physeq, x="treatment", fill = "Phylum") 
d=p1$data
head(d)
write.csv(d,"composition.txt")

#statistical tests
d=read.table("composition.txt", header = T, sep = "\t")
library(dplyr)
phylsum=d %>% group_by(treatment, Sample, Phylum) %>% summarise(Abundance = sum(Abundance))
head(phylsum)

phylsum$Phylum=as.factor(phylsum$Phylum)
gy.level<-levels(phylsum$Phylum)

mresult<-matrix(data=NA, nrow=length(gy.level), ncol=3)
colnames(mresult)=c("phylum", "W", "p-value")
for(i in (1:length(gy.level))){
  phylsum2=phylsum[which(phylsum$Phylum == gy.level[i]), ]
  mresult[i,1]=gy.level[i]
  mresult[i,2]=wilcox.test(Abundance ~ treatment, data=phylsum2, exact = FALSE, alternative = "greater")$statistic
  mresult[i,3]=wilcox.test(Abundance ~ treatment, data=phylsum2, exact = FALSE, alternative = "greater")$p.value
}

write.csv(mresult, "stat.phylum.oneside.wdecrease.abund.csv")

#################
#################
##################Draw figures

##########abundance
##Phylum level
unique(d$OTU)
a1=ggplot(d, aes(x = treatment, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values =colr)+
  labs(x = "Treatment", y = "Cumulative abundance (phylum level)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
a1



