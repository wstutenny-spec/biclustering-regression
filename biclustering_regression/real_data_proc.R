library(dplyr)
library(GEOquery)
#preprocessing gene data
#get phenotype data from csv file
gene_data_pheno=read.csv("/Users/wangshutu/Documents/GitHub/biclustering_regression/GSE3330_B6BTBRlivergenotype60mice.csv")
gene_data_pheno=gene_data_pheno[c(1:60),c(1:26)]#only keep phenotype
names(gene_data_pheno)=sub("X.","",names(gene_data_pheno))#remove special character


#get expression data from GSE3330
expr=getGEO("GSE3330")
df=exprs(expr[[1]])#get expression data
dim(df)
gene_data_exprs=t(df)

gene_data=as.data.frame(cbind(gene_data_pheno,gene_data_exprs))
save(gene_data,file="gene_data.Rdata")










#preprocessing afqt data
asvab=read.csv("/Users/wangshutu/Documents/GitHub/biclustering_regression/ASVAB.csv")
head(asvab);dim(asvab)
#rename variable
names(asvab)=c("ID","SEX","Birth_month","Birth_year","Survey_type","PIAT_97","Race","PIAT_99",
               "PIAT_98","PIAT_00","PIAT_01",
               "GS_pos","AR_pos","WK_pos","PC_pos","NO_pos","CS_pos","AI_pos","SI_pos","MK_pos","MC_pos","EI_pos","AO_pos",
               "GS_neg","AR_neg","WK_neg","PC_neg","NO_neg","CS_neg","AI_neg","SI_neg","MK_neg","MC_neg","EI_neg","AO_neg",
               "ACT_comp","ACT_eng","ACT_math","ACT_read","SAT_verbal","SAT_math","PIAT_02")
summary(asvab)

#PIAT
asvab$PIAT_97[asvab$PIAT_97%in%c(-1,-2,-3,-4,-5)]=NA
asvab$PIAT_98[asvab$PIAT_98%in%c(-1,-2,-3,-4,-5)]=NA
asvab$PIAT_99[asvab$PIAT_99%in%c(-1,-2,-3,-4,-5)]=NA
asvab$PIAT_00[asvab$PIAT_00%in%c(-1,-2,-3,-4,-5)]=NA
asvab$PIAT_01[asvab$PIAT_01%in%c(-1,-2,-3,-4,-5)]=NA
asvab$PIAT_02[asvab$PIAT_02%in%c(-1,-2,-3,-4,-5)]=NA

#ACT
asvab$ACT_comp[asvab$ACT_comp%in%c(-1,-2,-3,-4,-5)]=NA
asvab$ACT_eng[asvab$ACT_eng%in%c(-1,-2,-3,-4,-5)]=NA
asvab$ACT_math[asvab$ACT_math%in%c(-1,-2,-3,-4,-5)]=NA
asvab$ACT_read[asvab$ACT_read%in%c(-1,-2,-3,-4,-5)]=NA

#SAT
asvab$SAT_math[asvab$SAT_math%in%c(-1,-2,-3,-4,-5)]=NA
asvab$SAT_verbal[asvab$SAT_verbal%in%c(-1,-2,-3,-4,-5)]=NA

#ASVAB
for(i in 1:24){
  asvab[,(i+11)][asvab[,(i+11)]%in%c(-1,-2,-3,-5)]=NA
  asvab[,(i+11)][asvab[,(i+11)]==-4]=0
}

asvab=asvab%>%mutate(GS=GS_pos+GS_neg,AR=AR_pos+AR_neg,WK=WK_pos+WK_neg,PC=PC_pos+PC_neg,
                     NO=NO_pos+NO_neg,CS=CS_pos+CS_neg,AI=AI_pos+AI_neg,SI=SI_pos+SI_neg,
                     MK=MK_pos+MK_neg,MC=MC_pos+MC_neg,EI=EI_pos+EI_neg,AO=AO_pos+AO_neg)%>%
  select(GS:AO,ACT_comp:SAT_math)

#lots of missings in PIAT, just use ACT and SAT score as response
#maybe remove ACT comp as it is the average of other parts of ACT
summary(asvab)
dim(na.omit(asvab))


















