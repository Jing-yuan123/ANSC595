rm(list = ls())

source("https://bioconductor.org/biocLite.R")
biocLite("phyloseq")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")

library(ape)
library(ggplot2)
library(gplots)
library(lme4)
library(phangorn)
library(phyloseq)
library(plotly)
library(tidyr)
library(dplyr)
library(vegan)
library(VennDiagram)
library(readxl)

OTU = read.table("gene.opti_mcc.0.03.subsample.shared", header=TRUE, sep="\t")
tax = read.table("gene.taxonomy", header=TRUE, sep="\t")
meta = read.csv("metadata.csv")
data = read.table("gene.opti_mcc.groups.ave-std.summary", header=TRUE, sep="\t")

row.names(OTU) = OTU$Group
OTU.clean = OTU[,-which(names(OTU) %in% c("label", "numOtus", "Group"))]

row.names(tax) = tax$OTU
tax.clean = tax[row.names(tax) %in% colnames(OTU.clean),]
tax.clean$Taxonomy <- gsub(pattern = "[(]\\d*[)]", replacement = "", x=tax.clean$Taxonomy)
tax.clean = separate(tax.clean, Taxonomy, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"), sep=";")
tax.clean = tax.clean[,-which(names(tax.clean) %in% c("Size", "Strain", "OTU"))]


data.clean=data[-229:-456,]
#data.clean2=subset(data,method=="ave")
row.names(data.clean) = data.clean$group
data.clean = data.clean[,-which(names(data.clean) %in% c("group"))]

row.names(meta) = meta$Run
meta.clean = meta[,-which(names(meta) %in% c("Run"))]

OTU.clean = OTU.clean[order(row.names(OTU.clean)),]

meta = merge(meta.clean, data.clean, by.x =0, by.y = 0)
row.names(meta)=meta$Row.names
meta= meta[,-which(names(meta) %in% c("Row.names"))]

#metaseq.average <- meta %>% group_by(Calf.BG1) %>%summarise(overall.average=mean(nseqs))
#metaseq.average

OTU.new = merge(meta.clean, OTU, by.x =0, by.y = 0)
OTU.average <- OTU.new %>% group_by(Calf.BG1) %>%summarise(overall.average=mean(numOtus))
#OTU.average
OTU.new.clean =OTU.new[,-c(1:6)]
OTU1.clean =OTU.new.clean[,-c(2:11)]

str(meta)
meta$Calf.BG1= factor(meta$Calf.BG1, c("1","2","3","4","5","6"))
levels(meta$Calf.BG1) <- c("BG1", "BG2", "BG3","BG4", "BG5", "BG6")

OTU1.clean$Calf.BG1= factor(OTU1.clean$Calf.BG1, c("1","2","3","4","5","6"))
levels(OTU1.clean$Calf.BG1) <- c("BG1", "BG2", "BG3","BG4", "BG5", "BG6")

set.seed(8765)

graphics.off()
hist(meta$shannon, main="Shannon diversity", xlab="", breaks=10)
hist(meta$invsimpson, main="Invsimpson diversity", xlab="", breaks=10)
hist(meta$chao, main="Chao richness", xlab="", breaks=15)
hist(meta$ace, main="ACE richness", xlab="", breaks=15)

shapiro.test(meta$shannon)
shapiro.test(meta$invsimpson)
shapiro.test(meta$chao)
shapiro.test(meta$ace)

kruskal.test(shannon ~ Calf.BG1, data=meta)
pairwise.wilcox.test(meta$shannon, meta$Calf.BG1, p.adjust.method="fdr")
par(mfrow = c(1, 1))
#boxplot(shannon ~ Calf.BG1, data=meta, ylab="Shannon")
shannon_pd <- ggplot(meta, aes(Calf.BG1, shannon)) + 
  geom_boxplot(aes(color = Calf.BG1)) + 
  #ylim(c(0.5,1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Shannon", x = "") 
ggsave("C:/Users/yuanj/Desktop/microbiology/shannon.png", height = 3, width = 3)

kruskal.test(invsimpson ~ Calf.BG1, data=meta)
pairwise.wilcox.test(meta$invsimpson, meta$Calf.BG1, p.adjust.method="fdr")
par(mfrow = c(1, 1))
#boxplot(invsimpson ~ Calf.BG1, data=meta, ylab="Invsimpson")
invsimpson_pd <- ggplot(meta, aes(Calf.BG1, invsimpson)) + 
  geom_boxplot(aes(color = Calf.BG1)) + 
  #ylim(c(0.5,1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="invsimpson", x = "") 
ggsave("C:/Users/yuanj/Desktop/microbiology/invsimpson.png", height = 3, width = 3)

kruskal.test(chao ~ Calf.BG1, data=meta)
pairwise.wilcox.test(meta$chao, meta$Calf.BG1, p.adjust.method="fdr")
par(mfrow = c(1, 1))
#boxplot(chao ~ Calf.BG1, data=meta, ylab="chao")
chao_pd <- ggplot(meta, aes(Calf.BG1, chao)) + 
  geom_boxplot(aes(color = Calf.BG1)) + 
  #ylim(c(0.5,1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="chao", x = "") 
ggsave("C:/Users/yuanj/Desktop/microbiology/chao.png", height = 3, width = 3)

kruskal.test(ace ~ Calf.BG1, data=meta)
pairwise.wilcox.test(meta$ace, meta$Calf.BG1, p.adjust.method="fdr")
par(mfrow = c(1, 1))
#boxplot(ace ~ Calf.BG1, data=meta, ylab="ace richness")
ace_pd <- ggplot(meta, aes(Calf.BG1, ace)) + 
  geom_boxplot(aes(color = Calf.BG1)) + 
  #ylim(c(0.5,1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="ace", x = "") 
ggsave("C:/Users/yuanj/Desktop/microbiology/ace.png", height = 3, width = 3)




#------------------------------------------------------------------------------------------------------------------------------------------------------------
BC.nmds = metaMDS(OTU.clean, distance="bray", k=2, trymax=1000)
autotransform = TRUE
graphics.off()
plot(BC.nmds, type="n", main="Bray-Curtis")
points(BC.nmds, display="sites", pch=20, col=c("blue", "green", "red","yellow","black","pink")[meta$Calf.BG1])
legend(-3, 2, legend=c("BG1","BG2","BG3","BG4","BG5","BG6"), col=c("blue", "green", "red","yellow","black","pink"), pch=20)

BC.nmds$stress
nmds <-as.data.frame(BC.nmds$points)
metanmds <- merge(meta, nmds, by.x = 0, by.y = 0)
row.names(metanmds) <- metanmds[,1]
metanmds <- metanmds[,-1]
str(metanmds)

NMDS.mean <- aggregate(metanmds[,31:32], list(group=metanmds$Calf.BG1), mean)
colnames(NMDS.mean) <- c('design', 'groupX', 'groupY')
metanmds <- merge(metanmds, NMDS.mean, by.x = "Calf.BG1", by.y="design")
str(metanmds)

ggplot(metanmds, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(color=Calf.BG1)) +
  labs(x='NMDS 1', y= 'NMDS 2', caption = paste('Ordination stress: ', round(BC.nmds$stress, digits = 3))) +
  stat_ellipse(aes(color=Calf.BG1), level = 0.95) +
  theme(legend.title = element_blank()) 

ggsave("C:/Users/yuanj/Desktop/microbiology/nmds_ellipses_all.png", height = 3, width = 4)


#***************************************************************************************************************************************
J.nmds = metaMDS(OTU.clean, distance="jaccard", k=2, trymax=1000)
autotransform = TRUE
plot(J.nmds, type="n", main="Jaccard")
points(J.nmds, display="sites", pch=20, col=c("blue", "green", "red","yellow","black","pink")[meta$Calf.BG1])
legend(-2.5, 1.5, legend=c("BG1","BG2","BG3","BG4","BG5","BG6"), col=c("blue", "green", "red","yellow","black","pink"), pch=20)

J.nmds$stress
Jnmds <-as.data.frame(J.nmds$points)
Jmetanmds <- merge(meta, Jnmds, by.x = 0, by.y = 0)
row.names(Jmetanmds) <- Jmetanmds[,1]
Jmetanmds <- Jmetanmds[,-1]
str(Jmetanmds)

JNMDS.mean <- aggregate(Jmetanmds[,31:32], list(group=Jmetanmds$Calf.BG1), mean)
colnames(JNMDS.mean) <- c('design', 'JgroupX', 'JgroupY')
Jmetanmds <- merge(Jmetanmds, JNMDS.mean, by.x = "Calf.BG1", by.y="design")
str(Jmetanmds)

ggplot(Jmetanmds, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(color=Calf.BG1)) +
  labs(x='NMDS 1', y= 'NMDS 2', caption = paste('Ordination stress: ', round(J.nmds$stress, digits = 3))) +
  stat_ellipse(aes(color=Calf.BG1), level = 0.95) +
  theme(legend.title = element_blank()) 

ggsave("C:/Users/yuanj/Desktop/microbiology/Jnmds_ellipses_all.png", height = 3, width = 4)

#plot standard error (se) ellipses    these code did not work(don't use)
#plot(BC.nmds, type="n", main="Bray-Curtis")
#legend(-3, 2, legend=c("BG1","BG2","BG3","BG4","BG5","BG6"), col=c("blue", "green", "red","yellow","black","pink"), pch=20)
#ordiellipse(BC.nmds, groups=meta$Calf.BG1, display="sites", kind="se", conf=0.99, label=FALSE, col="blue", draw="polygon", alpha=200, show.groups = c("BG1"), border=FALSE)
#ordiellipse(BC.nmds, groups=meta$Calf.BG1, display="sites", kind="se", conf=0.99, label=FALSE, col="green", draw="polygon", alpha=200, show.groups = c("BG2"), border=FALSE)
#ordiellipse(BC.nmds, groups=meta$Calf.BG1, display="sites", kind="se", conf=0.99, label=FALSE, col="red", draw="polygon", alpha=200, show.groups = c("BG3"), border=FALSE)
#ordiellipse(BC.nmds, groups=meta$Calf.BG1, display="sites", kind="se", conf=0.99, label=FALSE, col="yellow", draw="polygon", alpha=200, show.groups = c("BG4"), border=FALSE)
#ordiellipse(BC.nmds, groups=meta$Calf.BG1, display="sites", kind="se", conf=0.99, label=FALSE, col="black", draw="polygon", alpha=200, show.groups = c("BG5"), border=FALSE)
#ordiellipse(BC.nmds, groups=meta$Calf.BG1, display="sites", kind="se", conf=0.99, label=FALSE, col="pink", draw="polygon", alpha=200, show.groups = c("BG6"), border=FALSE)

#3D plot
BC.nmds.3D = metaMDS(OTU.clean, distance="bray", k=3, trymax=1000)
BCxyz = scores(BC.nmds.3D, display="sites")
BCxyz
my_colors_3D = c("blue", "green", "red","yellow","black","pink")
plot_ly(x=BCxyz[,1], y=BCxyz[,2], z=BCxyz[,3], type="scatter3d", mode="markers", color=meta$Calf.BG1, colors=my_colors_3D)
par(mfrow=c(1,2))
plot(BCxyz[,1], BCxyz[,2], main="Bray-Curtis 1:2", pch=20, col=my_colors[meta$Calf.BG1])
legend(-5.4, 3, legend=c("BG1","BG2","BG3","BG4","BG5","BG6"), col=my_colors, pch=20)
plot(BCxyz[,1], BCxyz[,3], main="Bray-Curtis 1:3", pch=20, col=my_colors[meta$Calf.BG1])
plot(BCxyz[,2], BCxyz[,3], main="Bray-Curtis 2:3", pch=20, col=c("blue", "green", "red","yellow","black","pink")[meta$Calf.BG1])

#PERMANOVA
BC.dist=vegdist(OTU.clean, distance="bray")
adonis(BC.dist ~ Calf.BG1, data = meta, permutations = 1000)

J.dist=vegdist(OTU.clean, distance="jaccard")
adonis(J.dist ~ Calf.BG1, data = meta, permutations = 1000)

#*******************************************************************************************************************************************************
pairwise.adonis <- function(x, factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni')
{
  library(vegan)
  
  co = combn(unique(as.character(factors)),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  
  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    
    ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
  print("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
  return(pairw.res)
  
}
pairwise.adonis(OTU1.clean[,-1], OTU1.clean$Calf.BG1, sim.method = 'bray', p.adjust.m ='BH')

#----------------------------------------------------------------------------------------------------------------------------------------------------------
simper(OTU.clean, meta$Calf.BG1, permutations=100)
kruskal.test(OTU.clean$Otu00004 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00001 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00003 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00007 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00011 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00010 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00008 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00009 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00006 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00005 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00012 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00002 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00014 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00019 ~ meta$Calf.BG1)

kruskal.test(OTU.clean$Otu00023 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00025 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00030 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00022 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00013 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00021 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00024 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00026 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00016 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00015 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00020 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00074 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00029 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00027 ~ meta$Calf.BG1)

kruskal.test(OTU.clean$Otu00034 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00017 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00031 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00043 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00036 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00038 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00055 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00035 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00046 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00053 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00033 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00062 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00048 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00131 ~ meta$Calf.BG1)

kruskal.test(OTU.clean$Otu00104 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00051 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00032 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00052 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00047 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00064 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00082 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00061 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00057 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00069 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00058 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00054 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00049 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00018 ~ meta$Calf.BG1)

kruskal.test(OTU.clean$Otu00076 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00070 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00092 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00065 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00040 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00081 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00050 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00090 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00071 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00129 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00060 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00094 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00072 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00039 ~ meta$Calf.BG1)

kruskal.test(OTU.clean$Otu00028 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00044 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00045 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00200 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00067 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00037 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00086 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00154 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00041 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00102 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00059 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00078 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00108 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00087 ~ meta$Calf.BG1)

kruskal.test(OTU.clean$Otu00098 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00148 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00101 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00119 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00118 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00105 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00113 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00085 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00095 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00063 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00143 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00091 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00164 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00066 ~ meta$Calf.BG1)

kruskal.test(OTU.clean$Otu00123 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00042 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00173 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00281 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00193 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00075 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00117 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00169 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00068 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00077 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00120 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00107 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00114 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00083 ~ meta$Calf.BG1)

kruskal.test(OTU.clean$Otu00141 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00080 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00106 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00084 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00116 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00229 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00111 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00132 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00137 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00056 ~ meta$Calf.BG1)
kruskal.test(OTU.clean$Otu00088 ~ meta$Calf.BG1)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
OTU.UF = otu_table(as.matrix(OTU.clean), taxa_are_rows=FALSE)
tax.UF = tax_table(as.matrix(tax.clean))
meta.UF = sample_data(meta)

physeq = phyloseq(OTU.UF, tax.UF, meta.UF)
#install.packages("randomForest")
library(randomForest)
ntaxa(physeq)

prunescale = 0.0001
seq_depth = sample_sums(physeq)[1]
tax.mean <- taxa_sums(physeq)/nsamples(physeq)
rf.samples = c("BG1", "BG6")
physeq.prune <- prune_taxa(tax.mean > prunescale*seq_depth, physeq)
physeq.prune2 <- prune_samples(sample_data(physeq.prune)$Calf.BG1 == rf.samples[1] | sample_data(physeq.prune)$Calf.BG1 == rf.samples[2], physeq.prune)

physeq.prune2

# Make a dataframe of training data with OTUs as column and samples as rows
predictors <- (otu_table(physeq.prune2))
dim(predictors)

# Make one column for our outcome/response variable 
response <- as.factor(sample_data(physeq.prune2)$Calf.BG1)

# Combine them into 1 data frame
rf.data <- data.frame(response, predictors)

set.seed(2)

physeq.classify <- randomForest(response~., data = rf.data, ntree = 500)
print(physeq.classify)

names(physeq.classify)

# Make a data frame with predictor names and their importance
imp <- importance(physeq.classify)
imp <- data.frame(predictors = rownames(imp), imp)

# Order the predictor levels by importance
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
library(stringr)
imp.sort$predictors <- str_remove(imp.sort$predictors, "[X]")

# Select the top 10 predictors
imp.20 <- imp.sort[1:20, ]
imp.20$predictors <-  factor(imp.20$predictors, levels = imp.20$predictors)

# ggplot
ggplot(imp.20, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle(paste("Most important OTUs for classifying between", rf.samples[1], "and ", rf.samples[2]))

# What are those OTUs?
otunames <- imp.20$predictors
r <- rownames(tax_table(physeq)) %in% otunames
#rfkey=as.data.frame(tax_table(physeq)[r, ])
rf.imp.asv <- as.data.frame(tax_table(physeq)[r, ])
imp.20_tax <- merge(imp.20, rf.imp.asv, by.x = "predictors", by.y = 0)

my_colors_rf <- c(
  '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black"
)

# New ggplot
ggplot(imp.20_tax, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", aes(fill = Family)) +
  scale_fill_manual(values = my_colors_rf) +
  coord_flip() +
  ggtitle(paste("Most important ASVs between", rf.samples[1], "and", rf.samples[2]))
ggsave(paste0("RandomForest_", rf.samples[1], "_", rf.samples[2], "_names.png"), height = 4, width = 7)
write.table(imp.20_tax, "rf_top20_taxa.tsv", sep="\t", row.names=FALSE, quote = FALSE)


#-------------------------------------------------------------------------------------------------------------------------------------------------------------
plot_bar(physeq, fill="Phylum")  #don't use prune2
plot_bar(physeq, x="Calf.BG1", fill="Phylum") 
plot_bar(physeq, x="Calf.BG1", fill="Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

top5P.names = sort(tapply(taxa_sums(physeq), tax_table(physeq)[, "Phylum"], sum), TRUE)[1:5]
top5P = subset_taxa(physeq, Phylum %in% names(top5P.names))
plot_bar(top5P, x="Calf.BG1", fill="Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

#plot_bar(top5P, x="Calf.BG1", fill="Phylum", facet_grid = ~Phylum) + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

#There are many more options within `ggplot2` to alter this figure. [This document](https://www.rstudio.com/wp-content/uploads/2016/11/ggplot2-cheatsheet-2.1.pdf) has many helpful tips.

#These lines are slightly modified from the last version. These lines will use the otu table from the object `physeq` and the taxonomy table called tax.clean that we created in the very beginning of the script. The advantage of this taxonomy table is that it includes ASVs that are unclassified at any level. So if you have an ASV that you know is a "Bacilli" at the class level, but is unclassified after that, in your genus bar plot this ASV will be called "c__Bacilli". The "c" is for Class level. So with that designation you will know that it is unclassified at the genus level, but it was classified at the Class level as Bacilli. If you use the phyloseq taxonomy table, these unclassified tax are not reported at all.

physeq_otu_table <- data.frame(otu_table(physeq))
OTU.physeq = otu_table(as.matrix(physeq_otu_table), taxa_are_rows=FALSE)
tax.physeq = tax_table(as.matrix(tax.clean))    
meta.physeq = sample_data(meta)
physeq_bar_plot = phyloseq(OTU.physeq, tax.physeq, meta.physeq)


# Set colors for plotting
my_colors <- c(
  '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black"
)

my_level <- c("Phylum", "Family", "Genus")
rm(taxa.summary)

ml ="Phylum"
for(ml in my_level){
  print(ml)
  
taxa.summary <- physeq %>%
  tax_glom(taxrank = ml, NArm = FALSE) %>%  # agglomerate at `ml` level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt()  %>%                               # Melt to long format
  group_by(Calf.BG1, get(ml)) %>%
  summarise(Abundance.average=mean(Abundance)) 
names(taxa.summary)[2] <- ml

physeq.taxa.average <- taxa.summary %>% 
  group_by(get(ml)) %>%
  summarise(overall.average=mean(Abundance.average))
names(physeq.taxa.average)[1] <- ml

# merging the phyla means with the metadata #
physeq_meta <- merge(taxa.summary, physeq.taxa.average)

abund_filter <- 0.01
physeq_meta_filtered <- filter(physeq_meta, overall.average>abund_filter)
str(physeq_meta_filtered)

physeq_meta_filtered$Calf.BG1 = factor(physeq_meta_filtered$Calf.BG1, c("BG1", "BG2", "BG3", "BG4", "BG5","BG6"))

# Plot 
ggplot(physeq_meta_filtered, aes(x = Calf.BG1, y = Abundance.average, fill = get(ml))) + 
  #facet_grid(.~subject) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  # Remove x axis title
  #theme(axis.title.x = element_blank()) + 
  ylim(c(0,1)) +
  guides(fill = guide_legend(reverse = F, keywidth = .5, keyheight = .5, ncol = 1)) +
  theme(legend.text=element_text(size=8)) +
  #theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  ylab("Relative Abundance") +
  ggtitle(paste0(ml, " (>", abund_filter * 100,"%) Composition of microbiome samples")) 
ggsave(paste0(ml, "BarPlot_AllSamples.png"), height = 5)
}

#--------------------------------------------------------------------------------------------------------------------------------------------------------------

#heat map
#Sort the OTUs by abundance and pick the top 20
top20OTU.names = names(sort(taxa_sums(physeq), TRUE)[1:20])
#Cut down the physeq data to only the top 10 Phyla
top20OTU = prune_taxa(top20OTU.names, physeq)

top20OTU

plot_heatmap(top20OTU)
plot_heatmap(top20OTU, sample.label="Calf.BG1", sample.order="Calf.BG1")
ggsave(paste0( "top20_heatmap01.png"), height = 5, width=20)
plot_heatmap(top20OTU, sample.label="Calf.BG1", sample.order="Calf.BG1", taxa.label="Genus")+facet_grid(~Calf.BG1)
ggsave(paste0( "top20_heatmap02.png"), height = 5, width=20)
plot_heatmap(top20OTU, sample.label="Calf.BG1", sample.order="Calf.BG1", taxa.label="Genus", taxa.order="Phylum")
ggsave(paste0( "top20_heatmap03.png"), height = 5, width=20)
plot_heatmap(top20OTU, sample.label="Calf.BG1", sample.order="Calf.BG1", taxa.label="Genus", taxa.order="Phylum", low="white", high="purple", na.value="grey")
ggsave(paste0( "top20_heatmap04.png"), height = 5, width=20)
plot_heatmap(top20OTU, "NMDS", "bray", title="Bray-Curtis")

heatmap.2(as.matrix(BC.dist))
#Rainbow colors
rc <- rainbow(nrow(as.matrix(BC.dist)), start=0, end=0.9)
heatmap.2(as.matrix(BC.dist), col=rc)
ggsave(paste0( "top20_heatmapBCR.png"), height = 5, width=20)
