#-------------------------importing packages------------------------------

library(FactoMineR)
library(factoextra)
library(ggpubr)
library(ggplot2)

#---------------------importing and filtering data-------------------------

info <- read.table('info.tsv', sep='\t',header=TRUE)
groups<-split(info,info$group)
n_R<-dim(groups$R)[1]
n_NR<-dim(groups$NR)[1]
A <- read.table('A_homoCounts.tsv',sep='\t',row.names=1,header=TRUE)
B <- read.table('B_homoCounts.tsv',sep='\t',row.names=1,header=TRUE)
C <- read.table('C_homoCounts.tsv',sep='\t',row.names=1,header=TRUE)
all_pa<- read.table('count_a.tsv',sep='\t',row.names=1,header=TRUE)
sup_pa <- read.table('count_st.tsv',sep='\t',row.names=1,header=TRUE)
info <- info[,1:2]
hed1<-read.table("HED.tsv",sep='\t',header=TRUE)
hed2<-read.table("HED_cat.tsv",sep='\t',header=TRUE)
data <- merge(hed2, info, by.y=1,by.x=1, all.x = TRUE, all.y = TRUE)
data <- na.omit(data)
rownames(data) <- data$patient
data$patient <- NULL

#--------------------------------MCA--------------------------------------

#MCA analysis
res.mca <- MCA(data, quali.sup=4, graph=FALSE)
#results 
ind <- res.mca$ind
var <- res.mca$var
group <- res.mca$quali.sup

sink(file="MCA.txt", append = FALSE, type = "output")
print("individuals")
ind
print("variables")
var
print("groups")
group
sink(file = NULL)

#graphs

png("contribution.png",res=300, width=2000, height=2000)
fviz_screeplot(res.mca, addlabels = TRUE, ylim = c(0, 45))
dev.off()

png("variables.png",res=300, width=2000, height=2000)
fviz_mca_var(res.mca, repel=TRUE, col.var = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), ggtheme=theme_minimal())
dev.off()

png("groups.png",res=300, width=2000, height=2000)
fviz_mca_ind(res.mca,  habillage = "group",palette = c("#00AFBB", "#E7B800"),addEllipses = TRUE, repel= TRUE, ggtheme = theme_minimal())
dev.off()

png("biplot.png",res=300, width=2000, height=2000)
fviz_mca_biplot(res.mca, label="var", habillage = "group",palette = c("#00AFBB", "#E7B800"), repel = TRUE, graph.type = c("ggplot","classic"), addEllipses = TRUE)
dev.off()

#----------------------------fisher test-----------------------------------

sink(file ='homo_fisher.txt', append = FALSE, type = "output")
print("Gene A")
fisher.test(A)
print("Gene B")
fisher.test(B)
print("Gene C")
fisher.test(C)
sink(file = NULL)

#--------------------------------wilkoxon test------------------------------
row_sub=apply(as.matrix(hed1), 1, function(row) all(row !=0 ))
hed=hed1[row_sub,]
hed_stat <- function(gene,info,hed) { 
hed<-as.data.frame(hed)
df<-cbind(hed[[gene]],info$group)
df=as.data.frame(df)
names(df)[1]<- "distance"
names(df)[2] <- "group"
NR<- df[df$group == '1', ]$distance 
R<- df[df$group == '2', ]$distance
sink(file ='HED_wilcox.txt', append = TRUE, type = "output")
print(paste("Gene:",gene,sep=' '))
print(wilcox.test(distance ~ group, data = df, exact = FALSE))
sink(file = NULL)
}

hed_stat("A",info,hed)
hed_stat("B",info,hed)
hed_stat("C",info,hed)

#-----------------------------fisher on freq-------------------------------
pvalues<-c()
f <- function(y) {
d <- c(y[1],n_R-y[1],y[2],n_NR-y[2])
t=matrix(d, nrow = 2, ncol = 2)
m=matrix(unlist(t),2)
ft=fisher.test(m)
pvalues<-c(pvalues, ft$p.value)
return (pvalues)
}



sink(file ='freq_fisher.txt', append = TRUE, type = "output")
print('Alleles')
print(apply(all_pa,1, f))
print('Supertypes')
print(apply(sup_pa,1, f))
sink(file = NULL)
quit()
