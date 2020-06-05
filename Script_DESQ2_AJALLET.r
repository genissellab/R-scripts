# R version 3.6.2 (2019-12-12) -- "Dark and Stormy Night" // This file is written in the statistical language R
# Script By Arthur Jallet , finale version 4/June/2020
###########
# Outline: 
###########

	# 1. DESEq2 Model 2: line 27 to 322

	# 2. DESEq2 Model 3 (data MGGP01): line 323 to 1653
	
		# 2.a  Model fitting, creation of output data tables containing adj. p-values and fold-changes, and extraction of relevant gene sets: line 329 to 681
		# 2.b  Repeatability of evolution under fluctuations: line 682 to 701
		# 2.c  GO terms enrichment analyses: line 702 to 734
		# 2.d  Coexpression analysis (WGCNA clustering): line 735 to 888
		# 2.e  Analysis of expression plasticity in evolved lineages: line 889 to 942
		# 2.f  Normalized density of DE genes along the genome for each evolved lineages: line 943 to 1653

	# 3. DESeq2 Model 3 (data MGGP44): line 1654 to 2984
	
		# 3.a  Model fitting, creation of output data tables containing adj. p-values and fold-changes, and extraction of relevant gene sets: line 1660 to 2012
		# 3.b  Repeatability of evolution under fluctuations: line 2013 to 2032
		# 3.c  GO terms enrichment analyses: line 2033 to 2068
		# 3.d  Coexpression analysis (WGCNA clustering): line 2069 to 2220
		# 3.e  Analysis of expression plasticity in evolved lineages: line 2221 to 2274
		# 3.f  Normalized density of DE genes along the genome for each evolved lineages: line 2275 to 2984

#######################		
#######################
# 1. DESEq2 Model 2 (all 32 evolved samples)
#######################
#######################

# library(Hmisc)
source("C:/Users/Arthur/Desktop/Enrichment.r") # Home-made function for GO enrichment analysis (uses the "goseq" R package)

library(DESeq2)
Mat<-read.table("C:/Users/Arthur/Desktop/TabCount_11_juin_vrai.txt")
M=matrix(nrow=10972,ncol=40)
rownames(M)<-Mat[,1]
colnames(M)<-colnames(Mat[,2:41])
for (j in 1:40) {for (i in 1:10972){M[i,j]<-as.numeric(Mat[i,j+1])}}
dim(M)
r=unique(c(which(M[,1]>400000),which(M[,2]>400000),which(M[,3]>400000),which(M[,4]>400000),which(M[,5]>400000),which(M[,6]>400000),which(M[,7]>400000),which(M[,8]>400000),which(M[,9]>400000),which(M[,10]>400000),which(M[,11]>400000),which(M[,12]>400000),which(M[,13]>400000),which(M[,14]>400000),which(M[,15]>400000),which(M[,16]>400000),which(M[,17]>400000),which(M[,18]>400000),which(M[,19]>400000),which(M[,20]>400000),which(M[,21]>400000),which(M[,22]>400000),which(M[,23]>400000),which(M[,24]>400000),which(M[,25]>400000),which(M[,26]>400000),which(M[,27]>400000),which(M[,28]>400000),which(M[,29]>400000),which(M[,30]>400000), which(M[,31]>400000),which(M[,32]>400000),which(M[,33]>400000),which(M[,34]>400000),which(M[,35]>400000),which(M[,36]>400000),which(M[,37]>400000),which(M[,38]>400000),which(M[,39]>400000),which(M[,40]>400000)))
Juin<-read.table("C:/Users/Arthur/Desktop/data_11_juin.txt")
Outliers_chr=vector()
for (i in 1:38){Outliers_chr[i]=as.character(Juin[which(Juin[,10]==Mat[r,1][i])[1],1])}
Outliers_start=vector()
for (i in 1:38){Outliers_start[i]=as.character(Juin[which(Juin[,10]==Mat[r,1][i])[1],4])}
Outliers_end=vector()
for (i in 1:38){Outliers_end[i]=as.character(Juin[which(Juin[,10]==Mat[r,1][i])[1],5])}
Outliers_name=vector()
for (i in 1:38){Outliers_name[i]=as.character(Juin[which(Juin[,10]==Mat[r,1][i])[1],10])}
Outliers=cbind(Outliers_chr,Outliers_start,Outliers_end,Outliers_name)
M<-M[-c(which(is.element(Mat[,1],Outliers[(which(Outliers[,1]=="chr_7")),4][grep("67",Outliers[(which(Outliers[,1]=="chr_7")),4])])==TRUE)),]
dim(M)
M44<-M[,c(5:20,25:40)] # The four samples corresponding to ancestor MGGP01 (cols. 1 to 4) and the four samples corresponding to ancestor MGG44 (cols 21 to 24) are excluded
coldata<-read.table("C:/Users/Arthur/Desktop/coldata.txt")
coldata44<-coldata[c(5:20,25:40),1:3] # The four samples corresponding to ancestor MGGP01 and the four samples corresponding to ancestor MGG44 are excluded
dim(M44)
# [1] 10950 genes x 32 samples
# rownames(M44): contains the 10950 gene names (annotation from Grandaubert et al., 2011)

coldata44

          # Genotype Regime Tassay
# A1123171        #1   St23    T17
# A1123172        #1   St23    T17
# A1123234        #1   St23    T23
# A1123235        #1   St23    T23
# A1217171        #1   St17    T17
# A1217172        #1   St17    T17
# A1217234        #1   St17    T23
# A1217235        #1   St17    T23
# A12F171         #1  Fluct    T17
# A12F172         #1  Fluct    T17
# A12F234         #1  Fluct    T23
# A12F235         #1  Fluct    T23
# A13F171         #1  Fluct    T17
# A13F172         #1  Fluct    T17
# A13F234         #1  Fluct    T23
# A13F235         #1  Fluct    T23
# A44323171      #44   St23    T17
# A44323172      #44   St23    T17
# A44323234      #44   St23    T23
# A44323235      #44   St23    T23
# A44117171      #44   St17    T17
# A44117172      #44   St17    T17
# A44117234      #44   St17    T23
# A44117235      #44   St17    T23
# A441F171       #44  Fluct    T17
# A441F172       #44  Fluct    T17
# A441F234       #44  Fluct    T23
# A441F235       #44  Fluct    T23
# A443F171       #44  Fluct    T17
# A443F172       #44  Fluct    T17
# A443F234       #44  Fluct    T23
# A443F235       #44  Fluct    T23

des44 <- DESeqDataSetFromMatrix(countData=M44[,], colData=coldata44[,], design= ~ Genotype + Regime + Tassay + Genotype:Regime)
des44 <- estimateSizeFactors(des44)
des44 <- estimateDispersions(des44)

rld=rlog(des44)
library(ggplot2)
source("C:/Users/Arthur/Desktop/Plotty_pc_ANOVA.r") # Home-made function that performs a Principal Component Analysis
Plotty(rld, intgroup =c("Genotype","Tassay","Regime"),10950) # Log-transformed counts (10950 genes x 32 samples) are given as inputs
# The output is a vector called "d", that contains the computed PC scores along the first four Principal Components for the 32 samples
colnames(d)
# [1] "PC1" "PC2" "pca.x...3." "pca.x...4." "group" "Genotype" "Tassay" "Regime" "name" 

# Analysis of variance (ANOVA) on the 1st Principal Component, with the three qualitative covariates (Genotype, Tassay, Regime) included
PC1<-d[,c(1,6,7,8)]
mod0<-lm(PC1 ~ 1, data=PC1)
mod<-lm(PC1 ~ Genotype + Tassay + Regime, data=PC1)
anova(mod)
# Response: PC1
          # Df Sum Sq Mean Sq   F value  Pr(>F)    
# Genotype   1  35184   35184 2128.9145 < 2e-16 ***
# Tassay     1     64      64    3.8934 0.05879 .  
# Regime     2     93      47    2.8197 0.07725 .  
# Residuals 27    446      17                      
# ---
# See Table 2

# ANOVA on the 2nd Principal Component, with the three qualitative covariates (Genotype, Tassay, Regime) included
PC2<-d[,c(2,6,7,8)]
mod0<-lm(PC2 ~ 1, data=PC2)
mod<-lm(PC2 ~ Genotype + Tassay + Regime, data=PC2)
anova(mod)
# Response: PC2
          # Df Sum Sq Mean Sq F value   Pr(>F)    
# Genotype   1   59.8    59.8  0.4755   0.4964    
# Tassay     1 6664.5  6664.5 53.0074 7.86e-08 ***
# Regime     2    1.1     0.6  0.0044   0.9956    
# Residuals 27 3394.7   125.7                     
# ---
# See Table 2

# ANOVA on the 3rd Principal Component, with the three qualitative covariates (Genotype, Tassay, Regime) included
PC3<-d[,c(3,6,7,8)]
mod0<-lm(pca.x...3. ~ 1, data=PC3)
mod<-lm(pca.x...3. ~ Genotype + Tassay + Regime, data=PC3)
anova(mod)
# Response: pca.x...3.
          # Df  Sum Sq Mean Sq F value   Pr(>F)   
# Genotype   1   11.33   11.33  0.1221 0.729444   
# Tassay     1  975.56  975.56 10.5149 0.003144 **
# Regime     2 1474.70  737.35  7.9473 0.001932 **
# Residuals 27 2505.04   92.78   
# See Table 2

# ANOVA on the 4th Principal Component, with the three qualitative covariates (Genotype, Tassay, Regime) included
PC4<-d[,c(4,6,7,8)]
mod0<-lm(pca.x...4. ~ 1, data=PC4)
mod<-lm(pca.x...4. ~ Genotype + Tassay + Regime, data=PC4)
anova(mod)
# Response: pca.x...4.
          # Df  Sum Sq Mean Sq F value    Pr(>F)    
# Genotype   1   24.43   24.43  0.4977    0.4865    
# Tassay     1   58.64   58.64  1.1948    0.2840    
# Regime     2 2676.06 1338.03 27.2644 3.317e-07 ***
# Residuals 27 1325.05   49.08    
# See Table 2

Cont_all=matrix(nrow=10950,ncol=14) 
Cont_all_raw=matrix(nrow=10950,ncol=14)
Cont_all_l2fc=matrix(nrow=10950,ncol=14)
Cont_all[,1]<-rownames(M44)
Cont_all_raw[,1]<-rownames(M44)
Cont_all_l2fc[,1]<-rownames(M44)

# Figure 2
# "Genotype" effect (16 samples of MGGP01 background and 16 samples of MGGP44 background) 
# assesed by an LRT test comparing the full model versus the model with no Genotype effect included
des44 <- nbinomLRT(des44,reduced= ~ Regime + Tassay)
res44<-results(des44,name="Genotype_.44_vs_.1")
Cont_all[,2]=res44[,6]
Cont_all_raw[,2]=res44[,5]
Cont_all_l2fc[,2]=res44[,2]
length(which(as.numeric(Cont_all[,2])<0.05))
# [1] 7022 genes displaying a significant effect of the Genotype/Background (see the Venn diagram of Figure 2)
G_effect<-Cont_all[which(as.numeric(Cont_all[,2])<0.05),1]
# colnames(Cont_all)[2]<-c("Genotypic_effect")
# colnames(Cont_all_l2fc)[2]<-c("Genotypic_effect")
# colnames(Cont_all_raw)[2]<-c("Genotypic_effect")

# Figure 2
# "Temperature of assay" effect (16 samples assayed at 17° and 16 samples assayed at 23°C) 
# assesed by an LRT test comparing the full model versus the model with no Tassay effect included
des44 <- nbinomLRT(estimateDispersions(estimateSizeFactors(DESeqDataSetFromMatrix(countData=M44[,], colData=coldata44[,], design= ~ Genotype + Regime + Tassay + Genotype:Regime))),reduced= ~ Genotype + Regime + Genotype:Regime)
res44<-results(des44,name="Tassay_T23_vs_T17")
Cont_all[,3]=res44[,6]
Cont_all_raw[,3]=res44[,5]
Cont_all_l2fc[,3]=res44[,2]
length(which(as.numeric(Cont_all[,3])<0.05))
# 5637 genes displaying a significant effect of the Temperature of assay (see the Venn diagram of Figure 2)
T_effect<-Cont_all[which(as.numeric(Cont_all[,3])<0.05),1]
# colnames(Cont_all)[3]<-c("Temperature_of_assay_effect")
# colnames(Cont_all_l2fc)[3]<-c("Temperature_of_assay_effect")
# colnames(Cont_all_raw)[3]<-c("Temperature_of_assay_effect")

# Figure 2
# "Regime" effect (8 samples from the Stable_17°C regime, 8 samples from the Stable_23°C regime and 16 samples from the fluctuating regime) 
# assesed by an LRT test comparing the full model versus the model with no Regime effect included
des44 <- nbinomLRT(estimateDispersions(estimateSizeFactors(DESeqDataSetFromMatrix(countData=M44[,], colData=coldata44[,], design= ~ Genotype + Regime + Tassay + Genotype:Regime))),reduced= ~ Genotype + Tassay)
res44<-results(des44)
Cont_all[,4]=res44[,6]
Cont_all_raw[,4]=res44[,5]
Cont_all_l2fc[,4]=res44[,2]
length(which(as.numeric(Cont_all[,4])<0.05))
# [1] 4146 genes displaying a significant effect of the regime of selection (see the Venn diagram of Figure 2)
R_effect<-Cont_all[which(as.numeric(Cont_all[,4])<0.05),1]
# colnames(Cont_all)[4]<-c("Regime_effect")
# colnames(Cont_all_l2fc)[4]<-c("Regime_effect")
# colnames(Cont_all_raw)[4]<-c("Regime_effect")

# Drawing of the Figure 2
# install.packages('VennDiagram')
library(VennDiagram)
venn.diagram(
    x = list(G_effect , T_effect , R_effect),
    category.names = c("Genetic background effect" , "Temperature of assay effect \n " , "Selection regimes effect"),
    filename = "C:/users/ajallet/Desktop/#14_venn_diagramm.tiff",
    output = TRUE ,
    imagetype="tiff" ,
    height = 7000 , 
    width = 7000 , 
    resolution = 800,
    compression = "none",
    lwd = 2,
    lty = 'blank',
    fill = c('coral', 'purple', 'yellow2'),
    cex = 2,
    fontface = "bold",
    fontfamily = "sans",
    cat.cex = 1.4,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 180),
    cat.dist = c(0.055, 0.055, 0.035),
    cat.fontfamily = "sans",
    rotation = 1
)

# Testing of the interaction teme Genotype:Regime with an LRT test comparing the full model versus the model with no interaction included
des44 <- nbinomLRT(estimateDispersions(estimateSizeFactors(DESeqDataSetFromMatrix(countData=M44[,], colData=coldata44[,], design= ~ Genotype + Regime + Tassay + Genotype:Regime))),reduced= ~ Genotype + Regime + Tassay)
res44<-results(des44)
Cont_all[,8]=res44[,6]
Cont_all_raw[,8]=res44[,5]
Cont_all_l2fc[,8]=res44[,2]
length(which(as.numeric(Cont_all[,8])<0.05))
# [1] 1728 genes displaying a significant Genotype:Regime term

# The six single main contrasts comparing pairs of regimen within each genetic background (3 contratst for MGGP01; 3 contrasts for MGGP44)
des44 <- DESeqDataSetFromMatrix(countData=M44[,], colData=coldata44[,], design= ~ Genotype + Regime + Tassay + Genotype:Regime)
des44 <- estimateSizeFactors(des44)
des44 <- estimateDispersions(des44)
des44 <- nbinomWaldTest(des44)
resultsNames(des44)
# [1] "Intercept"              "Genotype_.44_vs_.1"     "Regime_St17_vs_Fluct"  
# [4] "Regime_St23_vs_Fluct"   "Tassay_T23_vs_T17"      "Genotype.44.RegimeSt17"
# [7] "Genotype.44.RegimeSt23"
 # ==> There will be seven values to provide in combination to calculate the considered contrast (of the Wald test).
 # R has taken "MGGP01" as the reference of the Genotypic effect, "T17" as the reference of the Tassay effect and "Fluct" as the reference of the Regime effect

# Effect St17 vs Fluct (only MGGP01 data considered)
res44<-results(des44,contrast=c(0,0,1,0,0,0,0))
Cont_all[,9]=res44[,6]
Cont_all_raw[,9]=res44[,5]
Cont_all_l2fc[,9]=res44[,2]

# Effect St17 vs Fluct (only MGGP44 data considered)
res44<-results(des44,contrast=c(0,0,1,0,0,1,0))
Cont_all[,10]=res44[,6]
Cont_all_raw[,10]=res44[,5]
Cont_all_l2fc[,10]=res44[,2]

pch=rep(176,length(which(as.numeric(Cont_all[,10])<0.05 | as.numeric(Cont_all[,9])<0.05)))
plot(as.numeric(Cont_all_l2fc[which(as.numeric(Cont_all[,10])<0.05 | as.numeric(Cont_all[,9])<0.05),10])[1:10950],as.numeric(Cont_all_l2fc[which(as.numeric(Cont_all[,10])<0.05 | as.numeric(Cont_all[,9])<0.05),9])[1:10950],ylim=c(-c(8),8),xlim=c(-c(8),20),xlab="Log2 Fold-Change (MGGP44_St17 / MGGP44_Fluct)",ylab=" Log2 Fold-Change (MGGP01_St17 / MGGP01_Fluct)",las=1,pch=pch,col="grey50")
abline(h=0)
abline(v=0)
abline(a=0,b=1,col="grey50")
for (i in which(is.element(Cont_all[which(as.numeric(Cont_all[,10])<0.05 | as.numeric(Cont_all[,9])<0.05),1],Cont_all[which(as.numeric(Cont_all[,10])<0.05 & as.numeric(Cont_all[,9])<0.05),1])==TRUE)) {points(as.numeric(Cont_all_l2fc[which(as.numeric(Cont_all[,10])<0.05 | as.numeric(Cont_all[,9])<0.05),10])[i],as.numeric(Cont_all_l2fc[which(as.numeric(Cont_all[,10])<0.05 | as.numeric(Cont_all[,9])<0.05),9])[i],col="red",pch=20)}

# Effect St23 vs Fluct (only MGGP01 data considered)
res44<-results(des44,contrast=c(0,0,0,1,0,0,0))
Cont_all[,11]=res44[,6]
Cont_all_raw[,11]=res44[,5]
Cont_all_l2fc[,11]=res44[,2]

# Effect St23 vs Fluct (only MGGP44 data considered)
res44<-results(des44,contrast=c(0,0,0,1,0,0,1))
Cont_all[,12]=res44[,6]
Cont_all_raw[,12]=res44[,5]
Cont_all_l2fc[,12]=res44[,2]

pch=rep(176,length(which(as.numeric(Cont_all[,12])<0.05 | as.numeric(Cont_all[,11])<0.05)))
plot(as.numeric(Cont_all_l2fc[which(as.numeric(Cont_all[,12])<0.05 | as.numeric(Cont_all[,11])<0.05),12])[1:10950],as.numeric(Cont_all_l2fc[which(as.numeric(Cont_all[,12])<0.05 | as.numeric(Cont_all[,11])<0.05),11])[1:10950],ylim=c(-c(8),8),xlim=c(-c(8),20),xlab="Log2 Fold-Change (MGGP44_St23 / MGGP44_Fluct)",ylab="Log2 Fold-Change (MGGP01_St23 / MGGP01_Fluct)",las=1,pch=pch,col="grey50")
abline(h=0)
abline(v=0)
abline(a=0,b=1,col="grey60")
for (i in which(is.element(Cont_all[which(as.numeric(Cont_all[,12])<0.05 | as.numeric(Cont_all[,11])<0.05),1],Cont_all[which(as.numeric(Cont_all[,12])<0.05 & as.numeric(Cont_all[,11])<0.05),1])==TRUE)) {points(as.numeric(Cont_all_l2fc[which(as.numeric(Cont_all[,12])<0.05 | as.numeric(Cont_all[,11])<0.05),12])[i],as.numeric(Cont_all_l2fc[which(as.numeric(Cont_all[,12])<0.05 | as.numeric(Cont_all[,11])<0.05),11])[i],col="red",pch=20)}

# Effect St17 vs St23 (only MGGP01 data considered)
res44<-results(des44,contrast=c(0,0,1,-1,0,0,0))
Cont_all[,13]=res44[,6]
Cont_all_raw[,13]=res44[,5]
Cont_all_l2fc[,13]=res44[,2]

# Effect St17 vs St23 (only MGGP44 data considered)
res44<-results(des44,contrast=c(0,0,1,-1,0,1,-1))
Cont_all[,14]=res44[,6]
Cont_all_raw[,14]=res44[,5]
Cont_all_l2fc[,14]=res44[,2]

pch=rep(176,length(which(as.numeric(Cont_all[,14])<0.05 | as.numeric(Cont_all[,13])<0.05)))
plot(as.numeric(Cont_all_l2fc[which(as.numeric(Cont_all[,14])<0.05 | as.numeric(Cont_all[,13])<0.05),14])[1:10950],as.numeric(Cont_all_l2fc[which(as.numeric(Cont_all[,14])<0.05 | as.numeric(Cont_all[,13])<0.05),13])[1:10950],ylim=c(-c(8),8),xlim=c(-c(8),20),xlab="Log2 Fold-Change (MGGP44_St17 / MGGP44_St23)",ylab="Log2 Fold-Change (MGGP01_St17 / MGGP01_St23)",las=1,pch=pch,col="grey50")
abline(h=0)
abline(v=0)
abline(a=0,b=1,col="grey60")
for (i in which(is.element(Cont_all[which(as.numeric(Cont_all[,14])<0.05 | as.numeric(Cont_all[,13])<0.05),1],Cont_all[which(as.numeric(Cont_all[,14])<0.05 & as.numeric(Cont_all[,13])<0.05),1])==TRUE)) {points(as.numeric(Cont_all_l2fc[which(as.numeric(Cont_all[,14])<0.05 | as.numeric(Cont_all[,13])<0.05),14])[i],as.numeric(Cont_all_l2fc[which(as.numeric(Cont_all[,14])<0.05 | as.numeric(Cont_all[,13])<0.05),13])[i],col="red",pch=20)}

#######################
#######################
# 2. DESEq2 Model 3, data MGGP01
#######################
#######################

#################
# 2.a Model fitting, creation of output data tables containing adj. p-values and fold-changes, and extraction of relevant gene sets 
#################

# library(Hmisc)
source("C:/Users/Arthur/Desktop/Enrichment.r") # Home-made function for GO enrichment analysis (uses the "goseq" R package)

rm(Cont_all)
rm(Cont_all_l2fc)
rm(Cont_all_raw)

library(DESeq2)
Mat<-read.table("C:/Users/Arthur/Desktop/TabCount_11_juin_vrai.txt")
M=matrix(nrow=10972,ncol=40)
rownames(M)<-Mat[,1]
colnames(M)<-colnames(Mat[,2:41])
for (j in 1:40) {for (i in 1:10972){M[i,j]<-as.numeric(Mat[i,j+1])}}
dim(M)
r=unique(c(which(M[,1]>400000),which(M[,2]>400000),which(M[,3]>400000),which(M[,4]>400000),which(M[,5]>400000),which(M[,6]>400000),which(M[,7]>400000),which(M[,8]>400000),which(M[,9]>400000),which(M[,10]>400000),which(M[,11]>400000),which(M[,12]>400000),which(M[,13]>400000),which(M[,14]>400000),which(M[,15]>400000),which(M[,16]>400000),which(M[,17]>400000),which(M[,18]>400000),which(M[,19]>400000),which(M[,20]>400000),which(M[,21]>400000),which(M[,22]>400000),which(M[,23]>400000),which(M[,24]>400000),which(M[,25]>400000),which(M[,26]>400000),which(M[,27]>400000),which(M[,28]>400000),which(M[,29]>400000),which(M[,30]>400000), which(M[,31]>400000),which(M[,32]>400000),which(M[,33]>400000),which(M[,34]>400000),which(M[,35]>400000),which(M[,36]>400000),which(M[,37]>400000),which(M[,38]>400000),which(M[,39]>400000),which(M[,40]>400000)))
Juin<-read.table("C:/Users/Arthur/Desktop/data_11_juin.txt")
Outliers_chr=vector()
for (i in 1:38){Outliers_chr[i]=as.character(Juin[which(Juin[,10]==Mat[r,1][i])[1],1])}
Outliers_start=vector()
for (i in 1:38){Outliers_start[i]=as.character(Juin[which(Juin[,10]==Mat[r,1][i])[1],4])}
Outliers_end=vector()
for (i in 1:38){Outliers_end[i]=as.character(Juin[which(Juin[,10]==Mat[r,1][i])[1],5])}
Outliers_name=vector()
for (i in 1:38){Outliers_name[i]=as.character(Juin[which(Juin[,10]==Mat[r,1][i])[1],10])}
Outliers=cbind(Outliers_chr,Outliers_start,Outliers_end,Outliers_name)
M<-M[-c(which(is.element(Mat[,1],Outliers[(which(Outliers[,1]=="chr_7")),4][grep("67",Outliers[(which(Outliers[,1]=="chr_7")),4])])==TRUE)),]
dim(M)
M44<-M[,1:20] # Because the first 20 columns correspond to MGGP01 data (cols. 1:4 = Anc; cols. 5:8: St23; cols. 9:12= St17; cols. 13:20: both lineages fluct.)
coldata<-read.table("C:/Users/Arthur/Desktop/coldata.txt")
levels(coldata$Regime)<-c("Anc","Fluct","St17","St23","FluctBis")
coldata$Regime[c(17:20,37:40)]="FluctBis" # columns 13 to 16 = first lineage evolved under fluct.; cols 17:20: second lineage evolved under fluct.
coldata44<-coldata[1:20,2:3] # Because the first 20 columns correspond to MGGP01 data
levels(coldata44$Regime)
coldata44$Regime
dim(coldata44)
dim(M44)
# [1] 10950 genes x 20 MGGP01 samples (16 evolved samples + 4 samples corresponding to the data of the ancestor)
# rownames(M44): contains the 10950 gene names (annotation from Grandaubert et al., 2011)

coldata44
          # Regime Tassay
# A01171        Anc    T17
# A01172        Anc    T17
# A01234        Anc    T23
# A01235        Anc    T23
# A1123171     St23    T17
# A1123172     St23    T17
# A1123234     St23    T23
# A1123235     St23    T23
# A1217171     St17    T17
# A1217172     St17    T17
# A1217234     St17    T23
# A1217235     St17    T23
# A12F171     Fluct    T17
# A12F172     Fluct    T17
# A12F234     Fluct    T23
# A12F235     Fluct    T23
# A13F171  FluctBis    T17
# A13F172  FluctBis    T17
# A13F234  FluctBis    T23
# A13F235  FluctBis    T23

# LRT test comparing the full model on MGGP01 data versus the model with no Regime effect included
des44 <- DESeqDataSetFromMatrix(countData=M44[,1:20], colData=coldata44[1:20,], design= ~ Regime + Tassay + Regime:Tassay)
des44 <- estimateSizeFactors(des44)
des44 <- estimateDispersions(des44)
des44 <- nbinomLRT(des44,reduced= ~Tassay )
res44<-results(des44)
Cont_all_01=matrix(nrow=10950,ncol=11)
Cont_all_01_raw=matrix(nrow=10950,ncol=11)
Cont_all_01_l2fc=matrix(nrow=10950,ncol=11)
Cont_all_01[,1]<-rownames(M44)
R_effect_01<-Cont_all_01[which(as.numeric(res44[,6])<0.05),1]
length(R_effect_01)
# [1] 3420

Cont_all_01=matrix(nrow=10950,ncol=11) # Will contain FDR-corrected p-values
Cont_all_01_raw=matrix(nrow=10950,ncol=11) # Will contain raw p-values
Cont_all_01_l2fc=matrix(nrow=10950,ncol=11) # Will contain Log2 Fold-Changes

# After gene names, the first ten columns of Cont_all_01 evaluates the following contrasts between lineages
colnames(Cont_all_01)<-c("Gene_id","St23_v_St17","Fluct_v_FluctBis","Fluct_v_St17","FluctBis_v_St17","Fluct_v_St23","FluctBis_v_St23","Anc_v_St17","Anc_v_St23","Anc_v_Fluct","Anc_v_FluctBis")

Cont_all_01[,1]<-rownames(M44)
Cont_all_01_raw[,1]<-rownames(M44)
Cont_all_01_l2fc[,1]<-rownames(M44)

for (i in 1:36){
Cont_all_01=cbind(Cont_all_01,rep(NA,10950))
Cont_all_01_raw=cbind(Cont_all_01_raw,rep(NA,10950))
Cont_all_01_l2fc=cbind(Cont_all_01_l2fc,rep(NA,10950))}

# The 10 main contrasts ([,2] to [,11]) comparing lineages (based on four samples for each, 2 obtained at Tassay 17°C & 2 at Tassay 23°C) with Wald tests

des44 <- DESeqDataSetFromMatrix(countData=M44[,1:20], colData=coldata44[1:20,], design= ~ Regime + Tassay + Regime:Tassay)
des44 <- estimateSizeFactors(des44)
des44 <- estimateDispersions(des44)
des44 <- nbinomWaldTest(des44)

resultsNames(des44)
 # [1] "Intercept"                "Regime_Fluct_vs_Anc"     
 # [3] "Regime_St17_vs_Anc"       "Regime_St23_vs_Anc"      
 # [5] "Regime_FluctBis_vs_Anc"   "Tassay_T23_vs_T17"       
 # [7] "RegimeFluct.TassayT23"    "RegimeSt17.TassayT23"    
 # [9] "RegimeSt23.TassayT23"     "RegimeFluctBis.TassayT23"
 # ==> There will be ten values to provide in combination to calculate the considered contrast (of the Wald test).
 # R has taken "Anc" as the reference of the Regime effect, "T17" as the reference of the Tassay effect

# Main effect St23 - St17
res44<-results(des44,contrast=c(0,0,-1,1,0,0,0,-0.5,0.5,0))
Cont_all_01[,2]=res44[,6]
Cont_all_01_raw[,2]=res44[,5]
Cont_all_01_l2fc[,2]=res44[,2]

# Main effect Fluct - FluctBis
res44<-results(des44,contrast=c(0,1,0,0,-1,0,0.5,0,0,-0.5))
Cont_all_01[,3]=res44[,6]
Cont_all_01_raw[,3]=res44[,5]
Cont_all_01_l2fc[,3]=res44[,2]

# Main effect Fluct - St17
res44<-results(des44,contrast=c(0,1,-1,0,0,0,0.5,-0.5,0,0))
Cont_all_01[,4]=res44[,6]
Cont_all_01_raw[,4]=res44[,5]
Cont_all_01_l2fc[,4]=res44[,2]

# Main effect FluctBis - St17
res44<-results(des44,contrast=c(0,0,-1,0,1,0,0,-0.5,0,0.5))
Cont_all_01[,5]=res44[,6]
Cont_all_01_raw[,5]=res44[,5]
Cont_all_01_l2fc[,5]=res44[,2]

# Main effect Fluct - St23
res44<-results(des44,contrast=c(0,1,0,-1,0,0,0.5,0,-0.5,0))
Cont_all_01[,6]=res44[,6]
Cont_all_01_raw[,6]=res44[,5]
Cont_all_01_l2fc[,6]=res44[,2]

# Main effect FluctBis - St23
res44<-results(des44,contrast=c(0,0,0,-1,1,0,0,0,-0.5,0.5))
Cont_all_01[,7]=res44[,6]
Cont_all_01_raw[,7]=res44[,5]
Cont_all_01_l2fc[,7]=res44[,2]

# Main effect Anc - St17
res44<-results(des44,contrast=c(0,0,-1,0,0,0,0,-0.5,0,0))
Cont_all_01[,8]=res44[,6]
Cont_all_01_raw[,8]=res44[,5]
Cont_all_01_l2fc[,8]=res44[,2]
length(which(as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),8])<0.05 & abs(as.numeric(Cont_all_01_l2fc[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),8]))>1))
# [1] 244 genes significantly DE (and with a FC >2) in Stable 17 lineage vs the MGGP01 ancestor

DE_S17vAnc_01<-Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),1][which(as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),8])<0.05 & abs(as.numeric(Cont_all_01_l2fc[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),8]))>1)]
median(2^(abs(as.numeric(Cont_all_01_l2fc[which(is.element(Cont_all_01[,1],DE_S17vAnc_01)==TRUE),8]))))
# [1] Median FC for these 244 genes: 3.077134 (see Table 3)

# Main effect Anc - St23
res44<-results(des44,contrast=c(0,0,0,-1,0,0,0,0,-0.5,0))
Cont_all_01[,9]=res44[,6]
Cont_all_01_raw[,9]=res44[,5]
Cont_all_01_l2fc[,9]=res44[,2]
length(which(as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),9])<0.05 & abs(as.numeric(Cont_all_01_l2fc[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),9]))>1))
# [1] 436 genes significantly DE (and with a FC >2) in Stable 23 lineage vs the MGGP01 ancestor

DE_S23vAnc_01<-Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),1][which(as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),9])<0.05 & abs(as.numeric(Cont_all_01_l2fc[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),9]))>1)]
median(2^(abs(as.numeric(Cont_all_01_l2fc[which(is.element(Cont_all_01[,1],DE_S23vAnc_01)==TRUE),9]))))
# [1] Median FC for these 436 genes: 3.084504 (see Table 3)

# Main effect Anc - Fluct
res44<-results(des44,contrast=c(0,-1,0,0,0,0,-0.5,0,0,0))
Cont_all_01[,10]=res44[,6]
Cont_all_01_raw[,10]=res44[,5]
Cont_all_01_l2fc[,10]=res44[,2]
length(which(as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),10])<0.05 & abs(as.numeric(Cont_all_01_l2fc[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),10]))>1))
# [1] 666 genes significantly DE (and with a FC >2) in Fluct_Rep1 lineage vs the MGGP01 ancestor

DE_F1vAnc_01<-Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),1][which(as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),10])<0.05 & abs(as.numeric(Cont_all_01_l2fc[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),10]))>1)]
median(2^(abs(as.numeric(Cont_all_01_l2fc[which(is.element(Cont_all_01[,1],DE_F1vAnc_01)==TRUE),10]))))
# [1] Median FC for these 666 genes: 2.935035 (see Table 3)

# Main effect Anc - FluctBis
res44<-results(des44,contrast=c(0,0,0,0,-1,0,0,0,0,-0.5))
Cont_all_01[,11]=res44[,6]
Cont_all_01_raw[,11]=res44[,5]
Cont_all_01_l2fc[,11]=res44[,2]
length(which(as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),11])<0.05 & abs(as.numeric(Cont_all_01_l2fc[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),11]))>1))
# [1] 809 genes significantly DE (and with a FC >2) in Fluct_Rep2 lineage vs the MGGP01 ancestor 

DE_F2vAnc_01<-Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),1][which(as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),11])<0.05 & abs(as.numeric(Cont_all_01_l2fc[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),11]))>1)]
median(2^(abs(as.numeric(Cont_all_01_l2fc[which(is.element(Cont_all_01[,1],DE_F2vAnc_01)==TRUE),11]))))
# [1] Median FC for these 809 genes: 2.714588 (see Table 3)

# LRT test comparing the full model on MGGP01 data versus the model with no Tassay effect included
des44 <- DESeqDataSetFromMatrix(countData=M44[,1:20], colData=coldata44[1:20,], design= ~ Regime + Tassay + Regime:Tassay)
des44 <- estimateSizeFactors(des44)
des44 <- estimateDispersions(des44)
des44 <- nbinomLRT(des44,reduced= ~ Regime )
res44<-results(des44)
T_effect<-Cont_all_01[which(as.numeric(res44[,6])<0.05),1]
length(T_effect)
# [1] 3601

des44 <- DESeqDataSetFromMatrix(countData=M44[,1:20], colData=coldata44[1:20,], design= ~ Regime + Tassay + Regime:Tassay)
des44 <- estimateSizeFactors(des44)
des44 <- estimateDispersions(des44)
des44 <- nbinomWaldTest(des44)

# Main effect T23 - T17 averaged over the five selection regimes
res44<-results(des44,contrast=c(0,0,0,0,0,1,0.2,0.2,0.2,0.2))
Cont_all_01[,12]=res44[,6]
Cont_all_01_raw[,12]=res44[,5]
Cont_all_01_l2fc[,12]=res44[,2]
colnames(Cont_all_01)[12]<-"T23_vs_T17"

# Delta T17-T23 in the Stable_17°C lineage vs Delta T17-T23 in the Stable_23°C lineage (interaction testing)
res44<-results(des44,contrast=c(0,0,0,0,0,0,0,-1,1,0))
Cont_all_01[,13]=res44[,6]
Cont_all_01_raw[,13]=res44[,5]
Cont_all_01_l2fc[,13]=res44[,2]

# Delta T17-T23 in the Fluct lineage vs Delta T17-T23 in the FluctBis lineage (interaction testing)
res44<-results(des44,contrast=c(0,0,0,0,0,0,1,0,0,-1))
Cont_all_01[,14]=res44[,6]
Cont_all_01_raw[,14]=res44[,5]
Cont_all_01_l2fc[,14]=res44[,2]

# Delta T17-T23 in the Fluct lineage vs Delta T17-T23 in the Stable_17°C lineage (interaction testing)
res44<-results(des44,contrast=c(0,0,0,0,0,0,-1,1,0,0))
Cont_all_01[,15]=res44[,6]
Cont_all_01_raw[,15]=res44[,5]
Cont_all_01_l2fc[,15]=res44[,2]

# Delta T17-T23 in the FluctBis lineage vs Delta T17-T23 in the Stable_17°C lineage (interaction testing)
res44<-results(des44,contrast=c(0,0,0,0,0,0,0,1,0,-1))
Cont_all_01[,16]=res44[,6]
Cont_all_01_raw[,16]=res44[,5]
Cont_all_01_l2fc[,16]=res44[,2]

# Delta T17-T23 in the Fluct lineage vs Delta T17-T23 in the Stable_23°C lineage (interaction testing)
res44<-results(des44,contrast=c(0,0,0,0,0,0,-1,0,1,0))
Cont_all_01[,17]=res44[,6]
Cont_all_01_raw[,17]=res44[,5]
Cont_all_01_l2fc[,17]=res44[,2]

# Delta T17-T23 in the FluctBis lineage vs Delta T17-T23 in the Stable_23°C lineage (interaction testing)
res44<-results(des44,contrast=c(0,0,0,0,0,0,0,0,1,-1))
Cont_all_01[,18]=res44[,6]
Cont_all_01_raw[,18]=res44[,5]
Cont_all_01_l2fc[,18]=res44[,2]

# Delta T17-T23 in the ancestor vs Delta T17-T23 in the the Stable_17°C lineage (interaction testing)
res44<-results(des44,contrast=c(0,0,0,0,0,0,0,-1,0,0))
Cont_all_01[,19]=res44[,6]
Cont_all_01_raw[,19]=res44[,5]
Cont_all_01_l2fc[,19]=res44[,2]

# Delta T17-T23 in the ancestor vs Delta T17-T23 in the the Stable_23°C lineage (interaction testing)
res44<-results(des44,contrast=c(0,0,0,0,0,0,0,0,-1,0))
Cont_all_01[,20]=res44[,6]
Cont_all_01_raw[,20]=res44[,5]
Cont_all_01_l2fc[,20]=res44[,2]

# Delta T17-T23 in the ancestor vs Delta T17-T23 in the the Fluct lineage (interaction testing)
res44<-results(des44,contrast=c(0,0,0,0,0,0,-1,0,0,0))
Cont_all_01[,21]=res44[,6]
Cont_all_01_raw[,21]=res44[,5]
Cont_all_01_l2fc[,21]=res44[,2]

# Delta T17-T23 in the ancestor vs Delta T17-T23 in the the FluctBis lineage (interaction testing)
res44<-results(des44,contrast=c(0,0,0,0,0,0,0,0,0,-1))
Cont_all_01[,22]=res44[,6]
Cont_all_01_raw[,22]=res44[,5]
Cont_all_01_l2fc[,22]=res44[,2]

# Effet of the Temp. of assay for the Fluct lineage
res44<-results(des44,contrast=c(0,0,0,0,0,1,1,0,0,0))
Cont_all_01[,23]=res44[,6]
Cont_all_01_raw[,23]=res44[,5]
Cont_all_01_l2fc[,23]=res44[,2]

# Effet of the Temp. of assay for the FluctBis lineage
res44<-results(des44,contrast=c(0,0,0,0,0,1,0,0,0,1))
Cont_all_01[,24]=res44[,6]
Cont_all_01_raw[,24]=res44[,5]
Cont_all_01_l2fc[,24]=res44[,2]

# Effet of the Temp. of assay for the Stable_17°C lineage
res44<-results(des44,contrast=c(0,0,0,0,0,1,0,1,0,0))
Cont_all_01[,25]=res44[,6]
Cont_all_01_raw[,25]=res44[,5]
Cont_all_01_l2fc[,25]=res44[,2]

# Effet of the Temp. of assay for the Stable_23°C lineage
res44<-results(des44,contrast=c(0,0,0,0,0,1,0,0,1,0))
Cont_all_01[,26]=res44[,6]
Cont_all_01_raw[,26]=res44[,5]
Cont_all_01_l2fc[,26]=res44[,2]

# Effet of the Temp. of assay for the ancestor
res44<-results(des44,contrast=c(0,0,0,0,0,1,0,0,0,0))
Cont_all_01[,27]=res44[,6]
Cont_all_01_raw[,27]=res44[,5]
Cont_all_01_l2fc[,27]=res44[,2]

des44 <- DESeqDataSetFromMatrix(countData=M44[,1:20], colData=coldata44[1:20,], design= ~ Regime + Tassay + Regime:Tassay)
des44 <- estimateSizeFactors(des44)
des44 <- estimateDispersions(des44)
des44 <- nbinomLRT(des44,reduced= ~ Regime + Tassay)
res44<-results(des44)
RbyT_effect<-Cont_all_01[which(as.numeric(res44[,6])<0.05),1]
length(RbyT_effect)
# [1] 414

A_vs_F_MGGP01<-Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),][which(as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),10])<0.05 & as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),11])<0.05 & abs(as.numeric(Cont_all_01_l2fc[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),10]))>1 & abs(as.numeric(Cont_all_01_l2fc[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),11]))>1),1]
length(A_vs_F_MGGP01)
# [1] Set of 568 genes consistently DE (adjusted p<0.05 & FC >2) in both MGGP01 fluctuating lineages in comparison to the ancestor.

# Preparing data to draw the plot of the Figure 4 (MGGP01)
mt=matrix(ncol=2,nrow=length(A_vs_F_MGGP01)*4)
mt[,2]=c(-c(as.numeric(Cont_all_01_l2fc[which(is.element(Cont_all_01[,1],A_vs_F_MGGP01)==TRUE),8]),as.numeric(Cont_all_01_l2fc[which(is.element(Cont_all_01[,1],A_vs_F_MGGP01)==TRUE),9]),as.numeric(Cont_all_01_l2fc[which(is.element(Cont_all_01[,1],A_vs_F_MGGP01)==TRUE),10]),as.numeric(Cont_all_01_l2fc[which(is.element(Cont_all_01[,1],A_vs_F_MGGP01)==TRUE),11])))
mt[1:length(A_vs_F_MGGP01),1]="St17"
mt[569:1136,1]="St23"
mt[1137:1704,1]="Fluct"
mt[1705:2272,1]="FluctBis"
mt<-as.data.frame(mt)
mt[,2]<-as.numeric(as.character(mt[,2]))
class(mt[,1])
# [1] "factor"
class(mt[,2])
# [1] "numeric"
table(mt[,1])

library(ggplot2)
library("ggbeeswarm")
par(xpd=TRUE, mar=c(6,6,4,4))
ggplot(mt, aes(factor(V1,levels=c("St17","St23","Fluct","FluctBis")), V2)) + geom_beeswarm(size=0.5,cex=0.5,col=c(rep("#8BEAF0",568),rep("#B3DE3A",568),rep("darkviolet",568),rep("brown1",568))) + scale_fill_manual(values=c("darkviolet", "brown1", "#8BEAF0", "#B3DE3A")) + theme(legend.position="none", axis.title.x=element_blank(), axis.text=element_text(size=7,face="bold")) + labs(y="Log2 Fold-Change in comparison to MGGP01 \n \n (for DEG under fluctuations)") +scale_x_discrete(labels = c("Stable 17°C","Stable 23°C","Fluct R_1","Fluct R_2")) + ylim(-10, 10)

# More precisely: 

length(which(as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],A_vs_F_MGGP01)==TRUE),8])<0.05))
# [1] 121 out of the 568 genes consistently DE under fluct (and with FC>2) are also DE at St17 ==> Represents 21 %, see Results in the article
# length(which(is.element(A_vs_F_MGGP01,DE_S17vAnc_01)==TRUE))
# [1] 109 out of the 568 genes are DE at St17 (here a FC>2 is imposed for St17 vs Anc)

length(which(as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],A_vs_F_MGGP01)==TRUE),9])<0.05))
# [1] 232 out of the 568 genes consistently DE under fluct (and with FC>2) are also DE at St23 ==> Represents 41 %, see Results in the article
# length(which(is.element(A_vs_F_MGGP01,DE_S23vAnc_01)==TRUE))
# [1] 197 out of the 568 genes are DE at St23 (here a FC>2 is imposed for St23 vs Anc)

#################
# 2.b Repeatability of evolution under fluctuations: F_1 vs F_2 (Figure 3)
#################

A_vs_F_MGGP01_noFCthreshold<-R_effect_01[which(as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),10])<0.05)[which(is.element(which(as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),10])<0.05),which(as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),11])<0.05))==TRUE)]]
# ==> 1517 genes consistently DE (whatever the fold-change) in both MGGP01 fluct. lineages versus the ancestor

#tiff(file="C:/Users/ajallet/Desktop/ex_14.tiff", height=7000, width=8000, res=1000)
par(mar=c(2,4,2,0)+0.1)
plot(-c(as.numeric(Cont_all_01_l2fc[,10])),-c(as.numeric(Cont_all_01_l2fc[,11])),pch=176,col="grey50",ylab="",xlab="",las=1,cex.axis=1.30,xlim=c(-8,8),ylim=c(-8,8))
# All data (10950 genes) in grey
points(-c(as.numeric(Cont_all_01_l2fc[which(is.element(Cont_all_01[,1],A_vs_F_MGGP01_noFCthreshold)==TRUE),10])),-c(as.numeric(Cont_all_01_l2fc[which(is.element(Cont_all_01[,1],A_vs_F_MGGP01_noFCthreshold)==TRUE),11])),col="red",pch=176)
# The 1517 ones in red
abline(v=1,lty="dashed")
abline(v=c(-1),lty="dashed")
abline(h=1,lty="dashed")
abline(h=c(-1),lty="dashed")
abline(a=0,b=1)
# dev.off()

#################
# 2.c GO terms enrichment analyses
#################

Enrichment(A_vs_F_MGGP01)
go.wall[which(as.numeric(go.wall[,2])<0.05),1:2]
# category over_represented_pvalue
# <0 lignes> 
# For the set of 568 genes consistently DE in both MGGP01 indep. replicates of evolution under fluct.(adjusted p<0.05 & FC >2) ==> No Enrichment

# Genes greyed in Table S2: those specifically DE under fluctuations (i.e. only DE in fluctuating lineages BUT not DE in any of the two stable lineages) + DE in fluct. with FC>2
Specif_F_MGGP01_SDE_above2Fold<-Cont_all_01[which(is.element(Cont_all_01[,1],A_vs_F_MGGP01)==TRUE),1][which(as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],A_vs_F_MGGP01)==TRUE),8])>0.05 & as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],A_vs_F_MGGP01)==TRUE),9])>0.05)]
length(Specif_F_MGGP01_SDE_above2Fold)
#[1] 291 such genes

Enrichment(Specif_F_MGGP01_SDE_above2Fold)
go.wall[which(as.numeric(go.wall[,2])<0.05),1:2]
# category                over_represented_pvalue
# <0 lignes> 
# No enrichment for these 291 genes

# Genes specifically DE under fluctuations (i.e. only DE in fluctuating lineages BUT not DE in any of the two stable lineages) but with FC<2 in at least one fluct lineages
Specif_F_MGGP01_SDE_below2Fold<-Cont_all_01[which(as.numeric(Cont_all_01[,8])>0.05 & as.numeric(Cont_all_01[,9])>0.05),1][which(is.element(Cont_all_01[which(as.numeric(Cont_all_01[,8])>0.05 & as.numeric(Cont_all_01[,9])>0.05),1],A_vs_F_MGGP01_noFCthreshold)==TRUE)][-c(which(is.element(Cont_all_01[which(as.numeric(Cont_all_01[,8])>0.05 & as.numeric(Cont_all_01[,9])>0.05),1][which(is.element(Cont_all_01[which(as.numeric(Cont_all_01[,8])>0.05 & as.numeric(Cont_all_01[,9])>0.05),1],A_vs_F_MGGP01_noFCthreshold)==TRUE)],A_vs_F_MGGP01)==TRUE))]
length(Specif_F_MGGP01_SDE_below2Fold)
#[1] 161 such genes

# Enrichment test on the 291 + 161 (these two sets do not overlap)
Enrichment(c(Specif_F_MGGP01_SDE_above2Fold,Specif_F_MGGP01_SDE_below2Fold))
go.wall[which(as.numeric(go.wall[,2])<0.05),1:2]
# category over_represented_pvalue
# 4 Cellcyclecontrol,celldivision,chromosomepartitioning            0.0159912855
# 9                                         Cytoskeleton            0.0196469753

#################
# 2.d Coexpression analysis (WGCNA clustering)
#################

library(WGCNA)

datExpr0<-t(M44)
datExpr<-datExpr0

# We indicate the Regime of selection and the Temperaturure of assay that correspond to each of the 20 samples
traitData<-matrix(nrow=20,ncol=2)
traitData[1:4,1]=as.character(as.factor("Anc"))
traitData[5:8,1]=as.character(as.factor("St23"))
traitData[9:12,1]=as.character(as.factor("St17"))
traitData[13:16,1]=as.character(as.factor("Fluct"))
traitData[17:20,1]=as.character(as.factor("FluctBis"))
traitData[,2]=rep(c(rep(c(as.character(as.factor("T17"))),2),rep(c(as.character(as.factor("T23"))),2)))
datTraits<-traitData

# Creation of a vector containing the genes for which one of the 21 constrasts (wether main or interaction effects) - and the corresponding LRT test for the tested covariate -
# are significant
Complete<-rownames(des44)[which(is.element(Cont_all_01[,1],unique(c(Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),1][which(as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),2])<0.05 | as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),3])<0.05 | as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),4])<0.05 | as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),5])<0.05 | as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),6])<0.05 | as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),7])<0.05 | as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),8])<0.05 | as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),9])<0.05 | as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),10])<0.05 | as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),11])<0.05)],Cont_all_01[which(is.element(Cont_all_01[,1],T_effect)==TRUE),1][which(as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],T_effect)==TRUE),12])<0.05)],Cont_all_01[which(is.element(Cont_all_01[,1],RbyT_effect)==TRUE),1][which(as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],RbyT_effect)==TRUE),13])<0.05 | as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],RbyT_effect)==TRUE),14])<0.05 | as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],RbyT_effect)==TRUE),15])<0.05 | as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],RbyT_effect)==TRUE),16])<0.05 | as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],RbyT_effect)==TRUE),17])<0.05 | as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],RbyT_effect)==TRUE),18])<0.05 | as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],RbyT_effect)==TRUE),19])<0.05 | as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],RbyT_effect)==TRUE),20])<0.05 | as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],RbyT_effect)==TRUE),21])<0.05 | as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],RbyT_effect)==TRUE),22])<0.05)])))==TRUE)]
length(Complete)
# [1] 5030

datExpr<-datExpr[,which(is.element(colnames(datExpr),Complete)==TRUE)]
dim(datExpr)
# [1]   20 5030

# Normalization step
des44 <- DESeqDataSetFromMatrix(countData=t(datExpr[,]), colData=coldata44[1:20,], design= ~ Regime + Tassay + Regime:Tassay)
des44 <- estimateSizeFactors(des44)
datExpr<-t(counts(des44, normalized=TRUE))

# Log 2 transformation of the data (before that, zeroes has to be removed)				
v=vector()
for (i in 1:20) {v=c(v,as.numeric(which(datExpr[i,]==0)))}
v<-unique(v)
datExpr<-datExpr[,-c(v)]
datExpr<-log2(datExpr)
dim(datExpr)
#   20 4948
# ==> Co-expression clustering will be performed by WGCNA on this dataset (i.e. on log2 values of the normalized counts of the 4947 genes contained in datExpr)

# gsg = goodSamplesGenes(datExpr, verbose = 3);
# gsg$allOK
# [1] TRUE

library(preprocessCore)
library(impute)
library(WGCNA)

# Choose a set of soft-thresholding powers
powers = c(1:20)
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
    main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
	 
-sign(sft$fitIndices[,3])*sft$fitIndices[,2]
 # [1] -0.72142119 -0.38772502 -0.02924733  0.09244198  0.32663473  0.50290155
 # [7]  0.63231796  0.71053986  0.76338481  0.78432446  0.81120724  0.82497943
# [13]  0.84173174  0.86063160  0.85979928  0.86200706  0.86285153  0.86681932
# [19]  0.87520770  0.88592578

# For a soft thresholding power of 12, R²=0.83 for scale-free topology test
softPower = 12;
adjacency = adjacency(datExpr, power = softPower,type="signed");

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
rm(TOM)

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");

# We set the minimum module size at 40:
minModuleSize = 40

# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
deepSplit = 2, pamRespectsDendro = FALSE,
minClusterSize = minModuleSize);

length(dynamicMods)
#[1] 4948

# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)

# Check here: if TRUE, it works
# dynamicColors==MEList$validColors
# unique(dynamicColors)==unique(MEList$validColors)

MEs = MEList$eigengenes
# View(MEs)
# One column = expression profile across the different experimental combinations of {Regime x Temp. of asaay} treatements dispalyed for the summary (i.e. the eignegene) of the genes 
# contained in the module of this column

# Clustering of module eigengenes:
MEDiss = 1-cor(MEs);

# Module merging based on a threshold of similmarity (here 75%)
MEDissThres = 0.25
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

# ==> After merging, 16 MGGP01 modules grouping genes according to their coexpression profiles across the experimental {Regime x Temp. of asaay} conditions
# With WGCNA, modules are automatically named according to their colours

names(mergedMEs)
# For each of the 16 eigengene, the name of the one of the module "colourX" is "MEcolourX"

for (i in 1:16) {print(paste(unique(mergedColors)[i],length(which(mergedColors==unique(mergedColors)[i]))),sep=" ")}
# [1] "darkolivegreen 622"
# [1] "darkorange 151"
# [1] "black 556"
# [1] "blue 760"
# [1] "red 356"
# [1] "lightyellow 125"
# [1] "orange 154"
# [1] "skyblue 456"
# [1] "darkgreen 458"
# [1] "purple 177"
# [1] "white 77"
# [1] "magenta 532"
# [1] "yellowgreen 41"
# [1] "tan 166"
# [1] "darkred 227"
# [1] "darkmagenta 90"
# The sum is equal to 4948 ==> Every genes fall into a module and genes are uniquely assigned to exactly one coexpression module

length(which(is.element(A_vs_F_MGGP01_noFCthreshold,c(rownames(t(datExpr))[which(mergedColors[]==unique(mergedColors)[5])],rownames(t(datExpr))[which(mergedColors[]==unique(mergedColors)[8])],rownames(t(datExpr))[which(mergedColors[]==unique(mergedColors)[2])],rownames(t(datExpr))[which(mergedColors[]==unique(mergedColors)[9])],rownames(t(datExpr))[which(mergedColors[]==unique(mergedColors)[10])]))==TRUE))
# ==> The five MGGP01 "red", "skyblue", "darkorange", "darkgreen" and "purple" modules contain ~63% of the genes that evolved under Fluctuations

# Within each module, genes are sorted by their correlation to the eigengene: the more a gene looks like the average (i.e. to the eigengene) the more it is a hub within the module.
# Hubs are genes that are connected to many other genes of their module and that probably drive the "coexpression behaviour" across conditions of all genes of the module
# //
# When ranked according to their correlation to the eigengene (= a measure of their centrality), are genes involved annotated as "Transcription factors"
#  or as involved in "Signal transduction" more central ? 

################
# 2.e Analysis of expression plasticity in evolved lineages in regard to the ancestor
################

Plastic_MGGP01_Anc<-rownames(M44)[which(is.element(rownames(M44),Cont_all_01[which(as.numeric(Cont_all_01[,27])<0.05 & abs(as.numeric(Cont_all_01_l2fc[,27]))>1),1])==TRUE)]
# [1] 501 genes displaying plasticity in gene expression (at a two-fold threshold) in the MGGP01 ancestor
Robust_MGGP01_Anc<-rownames(M44)[which(is.element(rownames(M44),Plastic_MGGP01_Anc)==FALSE)]

Plastic_MGGP01_S17<-rownames(M44)[which(is.element(rownames(M44),Cont_all_01[which(as.numeric(Cont_all_01[,25])<0.05 & abs(as.numeric(Cont_all_01_l2fc[,25]))>1),1])==TRUE)]
# [1] 2141 genes displaying plasticity in gene expression (at a two-fold threshold) in MGGP01 Stable_17°C lineage
Robust_MGGP01_S17<-rownames(M44)[which(is.element(rownames(M44),Plastic_MGGP01_S17)==FALSE)]

Plastic_MGGP01_S23<-rownames(M44)[which(is.element(rownames(M44),Cont_all_01[which(as.numeric(Cont_all_01[,26])<0.05 & abs(as.numeric(Cont_all_01_l2fc[,26]))>1),1])==TRUE)]
# [1] 609 genes displaying plasticity in gene expression (at a two-fold threshold) in MGGP01 Stable_23°C lineage
Robust_MGGP01_S23<-rownames(M44)[which(is.element(rownames(M44),Plastic_MGGP01_S23)==FALSE)]

Plastic_MGGP01_F1<-rownames(M44)[which(is.element(rownames(M44),Cont_all_01[which(as.numeric(Cont_all_01[,23])<0.05 & abs(as.numeric(Cont_all_01_l2fc[,23]))>1),1])==TRUE)]
# [1] 179 genes displaying plasticity in gene expression (at a two-fold threshold) in MGGP01 Fluct_1 lineage
Robust_MGGP01_F1<-rownames(M44)[which(is.element(rownames(M44),Plastic_MGGP01_F1)==FALSE)]

Plastic_MGGP01_F2<-rownames(M44)[which(is.element(rownames(M44),Cont_all_01[which(as.numeric(Cont_all_01[,24])<0.05 & abs(as.numeric(Cont_all_01_l2fc[,24]))>1),1])==TRUE)]
# [1] 303 genes displaying plasticity in gene expression (at a two-fold threshold) in MGGP01 Fluct_2 lineage
Robust_MGGP01_F2<-rownames(M44)[which(is.element(rownames(M44),Plastic_MGGP01_F2)==FALSE)]

# Genes that lost plasticty in expression during evolution (see Figure 6 Panel B and Table S9)
length(which(is.element(Plastic_MGGP01_Anc,Robust_MGGP01_S17)==TRUE))
# 307 such genes after evolution under the Stable_17°C regime
length(which(is.element(Plastic_MGGP01_Anc,Robust_MGGP01_S23)==TRUE))
# 341 such genes after evolution under the Stable_23°C regime
length(which(is.element(Plastic_MGGP01_Anc,Robust_MGGP01_F1)==TRUE))
# 444 such genes after evolution under the Fluctuating regime (rep. F1)
length(which(is.element(Plastic_MGGP01_Anc,Robust_MGGP01_F2)==TRUE))
# 424 such genes after evolution under the Fluctuating regime (rep. F2)

# Genes which expression became plastic during evolution (see Figure 6 Panel B and Table S9)
length(which(is.element(Robust_MGGP01_Anc,Plastic_MGGP01_S17)==TRUE))
# 1947 such genes after evolution under the Stable_17°C regime
length(which(is.element(Robust_MGGP01_Anc,Plastic_MGGP01_S23)==TRUE))
# 449 such genes after evolution under the Stable_23°C regime
length(which(is.element(Robust_MGGP01_Anc,Plastic_MGGP01_F1)==TRUE)=
# 122 such genes after evolution under the Fluctuating regime (rep. F1)
length(which(is.element(Robust_MGGP01_Anc,Plastic_MGGP01_F2)==TRUE))
# 226 such genes after evolution under the Fluctuating regime (rep. F2)

# Genes which expression remained plastic (i.e. those displaying plasticity in the ancestor and still in the considered evolved lineage)(see Figure 6 Panel B and Table S9)
length(which(is.element(Plastic_MGGP01_Anc,Plastic_MGGP01_S17)==TRUE))
# 194 such genes after evolution under the Stable_17°C regime
length(which(is.element(Plastic_MGGP01_Anc,Plastic_MGGP01_S23)==TRUE))
# 160 such genes after evolution under the Stable_23°C regime
length(which(is.element(Plastic_MGGP01_Anc,Plastic_MGGP01_F1)==TRUE))
# 57 such genes after evolution under the Fluctuating regime (rep. F1)
length(which(is.element(Plastic_MGGP01_Anc,Plastic_MGGP01_F2)==TRUE))
# 77 such genes after evolution under the Fluctuating regime (rep. F2)

###############
# 2.f Density of DE genes - normalized by gene density - along the genome for each of the four lineages of MGGP01 background (Figure 5)
###############

# Creation of the vector "density" that contains the density of annotated genes per 100 kb along the genome, with concatenated positions (from chromosome_1 to chromosome_21)
# Sliding window approach

name_bis<-read.table("C://Users/Arthur/Desktop/name_bis.txt",header=TRUE)
data0<-name_bis[,c(5,3,2)]
data<-data0[,c(3:2)]

X=vector()
for (i in 1:dim(data)[1]){X[i]=strsplit(as.character((data[,1])),split="chr_")[[i]][2]}
data<-cbind(as.numeric(X),data[,2])
data<-cbind(data,rep(NA,dim(data)[1]))

Y<-data

# Creation of data1, containing concatenated positions
for (i in 1:21){
    if (i==1){
        tmp=Y[Y[,1]==i,]
        tmp=tmp[order(tmp[,2]),]
        data1=tmp
        pos1=tmp[,2]
    }
    
    if (i>1){
        tmp=Y[Y[,1]==i,]
        tmp=tmp[order(tmp[,2]),]
        data1<-rbind(data1, tmp)
        pos=tmp[,2]
        max.value<-max(pos1)
        pos=pos+max.value
        pos1<-c(pos1, pos)
        data1[,3]<-pos1
    }        
}

dataA=vector()
for (i in 1:21)
{dataA[i]=tail(data1[which(data1[,1]==i),2])[6]}
# dataA

for (i in 1:21)

{chr<-unique(data0[,3])[i]
list.bin=seq(1000, tail(data1[which(data1[,1]==i),2])[6], by=1000)

tab=cbind(list.bin, rep(NA, length(list.bin)))
colnames(tab)=c("list.bin","count")

# will fill up col 2 with the number of events  
# use data1: concat positions of the events

k=1
v5=as.numeric(data1[which(data1[,1]==i),2])
for (i in list.bin){
    if(i==list.bin[1]){
        bin=v5[v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    if(i!=list.bin[1]){
        bin=v5[v5>lim & v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    
    k=k+1   
}

print(sum(tab[,2]))

require(zoo)
vec <- zoo(tab[,2])
sum(vec)
# [1] 1983

length(vec)
# [1] 6005

assign(paste("vec",chr,sep="_"), vec)

slided.vec<- rollapply(vec, width = 379, by = 5, FUN = mean)
assign(paste("slided.vec",chr,sep="_"), slided.vec)}

# Below: for the mini-chromosomes specifically (from chromosome_14 to chromosome_21)
# Creation of only one contiguous vector for all these eight chromosomes (necessary due to their small sizes and the chosen large window)

rm(X)
rm(data)
rm(data0)
rm(data1)

data0<-name_bis[which(name_bis[,2]=="chr_14" | name_bis[,2]=="chr_15" | name_bis[,2]=="chr_16" | name_bis[,2]=="chr_17" | name_bis[,2]=="chr_18" | name_bis[,2]=="chr_19" | name_bis[,2]=="chr_20" | name_bis[,2]=="chr_21"),c(5,3,2)]
data<-data0[,c(3:2)]

X=vector()
for (i in 1:dim(data)[1]){X[i]=strsplit(as.character((data[,1])),split="chr_")[[i]][2]}
data<-cbind(as.numeric(X),data[,2])
data<-cbind(data,rep(NA,dim(data)[1]))

Y<-data

# Creation of data1, containing concatenated positions
for (i in 14:21){
    if (i==14){
        tmp=Y[Y[,1]==i,]
        tmp=tmp[order(tmp[,2]),]
        data1=tmp
        pos1=tmp[,2]
    }
    
    if (i>14){
        tmp=Y[Y[,1]==i,]
        tmp=tmp[order(tmp[,2]),]
        data1<-rbind(data1, tmp)
        pos=tmp[,2]
        max.value<-max(pos1)
        pos=pos+max.value
        pos1<-c(pos1, pos)
        data1[,3]<-pos1
    }        
}

NN=tail(data1[,3])[6]

list.bin=seq(1000, tail(data1[,3])[6], by=1000)
tab=cbind(list.bin, rep(NA, length(list.bin)))
colnames(tab)=c("list.bin","count")

# will fill up col 2 with the number of events  
# use data1: concat positions of the events

k=1
v5=as.numeric(data1[,3])
for (i in list.bin){
    if(i==list.bin[1]){
        bin=v5[v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    if(i!=list.bin[1]){
        bin=v5[v5>lim & v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    
    k=k+1   
}

vec <- zoo(tab[,2])

slided.vec<- rollapply(vec, width = 379, by = 5, FUN = mean)

density=c(as.numeric(slided.vec_chr_1),as.numeric(slided.vec_chr_2),as.numeric(slided.vec_chr_3),as.numeric(slided.vec_chr_4),as.numeric(slided.vec_chr_5),as.numeric(slided.vec_chr_6),as.numeric(slided.vec_chr_7),as.numeric(slided.vec_chr_8),as.numeric(slided.vec_chr_9),as.numeric(slided.vec_chr_10),as.numeric(slided.vec_chr_11),as.numeric(slided.vec_chr_12),as.numeric(slided.vec_chr_13),as.numeric(slided.vec))
# This vector contains the density (per 100 kb) of annotated genes along the genome, with concatenated positions (from chromosome_1 to chromosome_21)

# Density of DE genes along the genome for the lineage MGGP01 F1 (normalized by the density of annotated genes)

NNN<-data1
end_density<-NN
rm(NN)
NN<-NNN

name_bis<-read.table("C://Users/Arthur/Desktop/name_bis.txt",header=TRUE)
name_bis<-name_bis[which(is.element(name_bis[,5],Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),1][which(as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),10])<0.05 & abs(as.numeric(Cont_all_01_l2fc[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),10]))>1)])==TRUE),]
# Les 666 DE chez MGGP01 F_1
dim(name_bis)
# [1] 666   7

new_NN<-NN[which(is.element(NN[,2],name_bis[which(is.element(name_bis[,3],NN[,2])==TRUE),3])==TRUE),]

data<-name_bis[,c(2,3)]

X=vector()
for (i in 1:dim(data)[1]){X[i]=strsplit(as.character((data[,1])),split="chr_")[[i]][2]}
data<-cbind(as.numeric(X),data[,2])
data<-cbind(data,rep(NA,dim(data)[1]))

Y<-data

for (i in 1:13){
    if (i==1){
        tmp=Y[Y[,1]==i,]
        tmp=tmp[order(tmp[,2]),]
        data1=tmp
        pos1=tmp[,2]
    }
    
    if (i>1){
        tmp=Y[Y[,1]==i,]
        tmp=tmp[order(tmp[,2]),]
        data1<-rbind(data1, tmp)
        pos=tmp[,2]
        max.value<-max(pos1)
        pos=pos+max.value
        pos1<-c(pos1, pos)
        data1[,3]<-pos1
    }        
}
# Contient les positions concaténées tout bien comme y faut pour les 584 (666-82) gènes des CC

for (i in 1:13)
    
{chr<-paste("chr",unique(data1[,1])[i],sep="_")
list.bin=seq(1000, dataA[i], by=1000)

tab=cbind(list.bin, rep(NA, length(list.bin)))
colnames(tab)=c("list.bin","count")

# will fill up col 2 with the number of events  
# use data1: concat positions of the events

k=1
v5=as.numeric(data1[which(data1[,1]==i),2])
for (i in list.bin){
    if(i==list.bin[1]){
        bin=v5[v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    if(i!=list.bin[1]){
        bin=v5[v5>lim & v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    
    k=k+1   
}

print(sum(tab[,2]))

require(zoo)
vec <- zoo(tab[,2])
sum(vec)
# [1] 1983

length(vec)
# [1] 6005

assign(paste("vec",chr,sep="_"), vec)

slided.vec<- rollapply(vec, width = 379, by = 5, FUN = mean)
assign(paste("slided.vec",chr,sep="_"), slided.vec)}

list.bin=seq(1000, end_density, by=1000)
tab=cbind(list.bin, rep(NA, length(list.bin)))
colnames(tab)=c("list.bin","count")

k=1
v5=as.numeric(new_NN[,3])
for (i in list.bin){
    if(i==list.bin[1]){
        bin=v5[v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    if(i!=list.bin[1]){
        bin=v5[v5>lim & v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    
    k=k+1   
}

print(sum(tab[,2]))
# 650

require(zoo)
vec <- zoo(tab[,2])
sum(vec)
# 650

length(vec)
# [1] 4403

slided.vec<- rollapply(vec, width = 379, by = 5, FUN = mean)

dede_F1<-c(as.numeric(slided.vec_chr_1),as.numeric(slided.vec_chr_2),as.numeric(slided.vec_chr_3),as.numeric(slided.vec_chr_4),as.numeric(slided.vec_chr_5),as.numeric(slided.vec_chr_6),as.numeric(slided.vec_chr_7),as.numeric(slided.vec_chr_8),as.numeric(slided.vec_chr_9),as.numeric(slided.vec_chr_10),as.numeric(slided.vec_chr_11),as.numeric(slided.vec_chr_12),as.numeric(slided.vec_chr_13),c(slided.vec))

density<-read.table("C:/Users/Arthur/Desktop/Density.txt")
plot(dede_F1/as.numeric(density[,1]),type="l",col="darkviolet",ylab="",xlab="",xaxt="n", lwd=2)

# Density of DE genes along the genome for the lineage MGGP01 F2 (normalized by the density of annotated genes)

name_bis<-read.table("C://Users/Arthur/Desktop/name_bis.txt",header=TRUE)
name_bis<-name_bis[which(is.element(name_bis[,5],Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),1][which(as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),11])<0.05 & abs(as.numeric(Cont_all_01_l2fc[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),11]))>1)])==TRUE),]
# Les 809 DE chez MGGP01 F_2
dim(name_bis)
# [1] 809   7

new_NN<-NN[which(is.element(NN[,2],name_bis[which(is.element(name_bis[,3],NN[,2])==TRUE),3])==TRUE),]

data<-name_bis[,c(2,3)]

X=vector()
for (i in 1:dim(data)[1]){X[i]=strsplit(as.character((data[,1])),split="chr_")[[i]][2]}
data<-cbind(as.numeric(X),data[,2])
data<-cbind(data,rep(NA,dim(data)[1]))

Y<-data

for (i in 1:13){
    if (i==1){
        tmp=Y[Y[,1]==i,]
        tmp=tmp[order(tmp[,2]),]
        data1=tmp
        pos1=tmp[,2]
    }
    
    if (i>1){
        tmp=Y[Y[,1]==i,]
        tmp=tmp[order(tmp[,2]),]
        data1<-rbind(data1, tmp)
        pos=tmp[,2]
        max.value<-max(pos1)
        pos=pos+max.value
        pos1<-c(pos1, pos)
        data1[,3]<-pos1
    }        
}

for (i in 1:13)
    
{chr<-paste("chr",unique(data1[,1])[i],sep="_")
list.bin=seq(1000, dataA[i], by=1000)

tab=cbind(list.bin, rep(NA, length(list.bin)))
colnames(tab)=c("list.bin","count")

# will fill up col 2 with the number of events  
# use data1: concat positions of the events

k=1
v5=as.numeric(data1[which(data1[,1]==i),2]) 
for (i in list.bin){
    if(i==list.bin[1]){
        bin=v5[v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    if(i!=list.bin[1]){
        bin=v5[v5>lim & v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    
    k=k+1   
}

print(sum(tab[,2]))

require(zoo)
vec <- zoo(tab[,2])
sum(vec)
# [1] 1983

length(vec)
# [1] 6005

assign(paste("vec",chr,sep="_"), vec)

slided.vec<- rollapply(vec, width = 379, by = 5, FUN = mean)
assign(paste("slided.vec",chr,sep="_"), slided.vec)}

list.bin=seq(1000, end_density, by=1000)
tab=cbind(list.bin, rep(NA, length(list.bin)))
colnames(tab)=c("list.bin","count")

k=1
v5=as.numeric(new_NN[,3])
for (i in list.bin){
    if(i==list.bin[1]){
        bin=v5[v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    if(i!=list.bin[1]){
        bin=v5[v5>lim & v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    
    k=k+1   
}

print(sum(tab[,2]))
# 650

require(zoo)
vec <- zoo(tab[,2])
sum(vec)
# 650

length(vec)
# [1] 4403

slided.vec<- rollapply(vec, width = 379, by = 5, FUN = mean)

dede_F2<-c(as.numeric(slided.vec_chr_1),as.numeric(slided.vec_chr_2),as.numeric(slided.vec_chr_3),as.numeric(slided.vec_chr_4),as.numeric(slided.vec_chr_5),as.numeric(slided.vec_chr_6),as.numeric(slided.vec_chr_7),as.numeric(slided.vec_chr_8),as.numeric(slided.vec_chr_9),as.numeric(slided.vec_chr_10),as.numeric(slided.vec_chr_11),as.numeric(slided.vec_chr_12),as.numeric(slided.vec_chr_13),c(slided.vec))

points(dede_F2/as.numeric(density[,1]),type="l",col="brown1",ylab="",xlab="",xaxt="n",lwd=2)

# Density of DE genes along the genome for the lineage MGGP01 Stable_17°C (normalized by the density of annotated genes)

name_bis<-read.table("C://Users/Arthur/Desktop/name_bis.txt",header=TRUE)
name_bis<-name_bis[which(is.element(name_bis[,5],Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),1][which(as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),8])<0.05 & abs(as.numeric(Cont_all_01_l2fc[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),8]))>1)])==TRUE),]
# Les 244 DE chez MGGP01 Stable_17°C
dim(name_bis)
# [1] 244   7

new_NN<-NN[which(is.element(NN[,2],name_bis[which(is.element(name_bis[,3],NN[,2])==TRUE),3])==TRUE),]

data<-name_bis[,c(2,3)]

X=vector()
for (i in 1:dim(data)[1]){X[i]=strsplit(as.character((data[,1])),split="chr_")[[i]][2]}
data<-cbind(as.numeric(X),data[,2])
data<-cbind(data,rep(NA,dim(data)[1]))

Y<-data

for (i in 1:13){
    if (i==1){
        tmp=Y[Y[,1]==i,]
        tmp=tmp[order(tmp[,2]),]
        data1=tmp
        pos1=tmp[,2]
    }
    
    if (i>1){
        tmp=Y[Y[,1]==i,]
        tmp=tmp[order(tmp[,2]),]
        data1<-rbind(data1, tmp)
        pos=tmp[,2]
        max.value<-max(pos1)
        pos=pos+max.value
        pos1<-c(pos1, pos)
        data1[,3]<-pos1
    }        
}

for (i in 1:13)
    
{chr<-paste("chr",unique(data1[,1])[i],sep="_")
list.bin=seq(1000, dataA[i], by=1000)

tab=cbind(list.bin, rep(NA, length(list.bin)))
colnames(tab)=c("list.bin","count")

# will fill up col 2 with the number of events  
# use data1: concat positions of the events

k=1
v5=as.numeric(data1[which(data1[,1]==i),2])
for (i in list.bin){
    if(i==list.bin[1]){
        bin=v5[v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    if(i!=list.bin[1]){
        bin=v5[v5>lim & v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    
    k=k+1   
}

print(sum(tab[,2]))

require(zoo)
vec <- zoo(tab[,2])
sum(vec)
# [1] 1983

length(vec)
# [1] 6005

assign(paste("vec",chr,sep="_"), vec)

slided.vec<- rollapply(vec, width = 379, by = 5, FUN = mean)
assign(paste("slided.vec",chr,sep="_"), slided.vec)}

list.bin=seq(1000, end_density, by=1000)
tab=cbind(list.bin, rep(NA, length(list.bin)))
colnames(tab)=c("list.bin","count")

k=1
v5=as.numeric(new_NN[,3])
for (i in list.bin){
    if(i==list.bin[1]){
        bin=v5[v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    if(i!=list.bin[1]){
        bin=v5[v5>lim & v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    
    k=k+1   
}

print(sum(tab[,2]))
# 650

require(zoo)
vec <- zoo(tab[,2])
sum(vec)
# 650

length(vec)
# [1] 4403

slided.vec<- rollapply(vec, width = 379, by = 5, FUN = mean)

dede_S17<-c(as.numeric(slided.vec_chr_1),as.numeric(slided.vec_chr_2),as.numeric(slided.vec_chr_3),as.numeric(slided.vec_chr_4),as.numeric(slided.vec_chr_5),as.numeric(slided.vec_chr_6),as.numeric(slided.vec_chr_7),as.numeric(slided.vec_chr_8),as.numeric(slided.vec_chr_9),as.numeric(slided.vec_chr_10),as.numeric(slided.vec_chr_11),as.numeric(slided.vec_chr_12),as.numeric(slided.vec_chr_13),c(slided.vec))

points(dede_S17/as.numeric(density[,1]),type="l",col="#8BEAF0",ylab="",xlab="",xaxt="n",lwd=2)

# Density of DE genes along the genome for the lineage MGGP01 Stable_23°C (normalized by the density of annotated genes)

name_bis<-read.table("C://Users/Arthur/Desktop/name_bis.txt",header=TRUE)
name_bis<-name_bis[which(is.element(name_bis[,5],Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),1][which(as.numeric(Cont_all_01[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),9])<0.05 & abs(as.numeric(Cont_all_01_l2fc[which(is.element(Cont_all_01[,1],R_effect_01)==TRUE),9]))>1)])==TRUE),]
# Les 436 DE chez MGGP01 Stable_23°C
dim(name_bis)
# [1] 436   7

new_NN<-NN[which(is.element(NN[,2],name_bis[which(is.element(name_bis[,3],NN[,2])==TRUE),3])==TRUE),]

data<-name_bis[,c(2,3)]

X=vector()
for (i in 1:dim(data)[1]){X[i]=strsplit(as.character((data[,1])),split="chr_")[[i]][2]}
data<-cbind(as.numeric(X),data[,2])
data<-cbind(data,rep(NA,dim(data)[1]))

Y<-data

for (i in 1:13){
    if (i==1){
        tmp=Y[Y[,1]==i,]
        tmp=tmp[order(tmp[,2]),]
        data1=tmp
        pos1=tmp[,2]
    }
    
    if (i>1){
        tmp=Y[Y[,1]==i,]
        tmp=tmp[order(tmp[,2]),]
        data1<-rbind(data1, tmp)
        pos=tmp[,2]
        max.value<-max(pos1)
        pos=pos+max.value
        pos1<-c(pos1, pos)
        data1[,3]<-pos1
    }        
}

for (i in 1:13)
    
{chr<-paste("chr",unique(data1[,1])[i],sep="_")
list.bin=seq(1000, dataA[i], by=1000)

tab=cbind(list.bin, rep(NA, length(list.bin)))
colnames(tab)=c("list.bin","count")

# will fill up col 2 with the number of events  
# use data1: concat positions of the events

k=1
v5=as.numeric(data1[which(data1[,1]==i),2])
for (i in list.bin){
    if(i==list.bin[1]){
        bin=v5[v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    if(i!=list.bin[1]){
        bin=v5[v5>lim & v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    
    k=k+1   
}

print(sum(tab[,2]))

require(zoo)
vec <- zoo(tab[,2])
sum(vec)
# [1] 1983

length(vec)
# [1] 6005

assign(paste("vec",chr,sep="_"), vec)

slided.vec<- rollapply(vec, width = 379, by = 5, FUN = mean)
assign(paste("slided.vec",chr,sep="_"), slided.vec)}

list.bin=seq(1000, end_density, by=1000)
tab=cbind(list.bin, rep(NA, length(list.bin)))
colnames(tab)=c("list.bin","count")

k=1
v5=as.numeric(new_NN[,3])
for (i in list.bin){
    if(i==list.bin[1]){
        bin=v5[v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    if(i!=list.bin[1]){
        bin=v5[v5>lim & v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    
    k=k+1   
}

print(sum(tab[,2]))
# 650

require(zoo)
vec <- zoo(tab[,2])
sum(vec)
# 650

length(vec)
# [1] 4403

slided.vec<- rollapply(vec, width = 379, by = 5, FUN = mean)

dede_S23<-c(as.numeric(slided.vec_chr_1),as.numeric(slided.vec_chr_2),as.numeric(slided.vec_chr_3),as.numeric(slided.vec_chr_4),as.numeric(slided.vec_chr_5),as.numeric(slided.vec_chr_6),as.numeric(slided.vec_chr_7),as.numeric(slided.vec_chr_8),as.numeric(slided.vec_chr_9),as.numeric(slided.vec_chr_10),as.numeric(slided.vec_chr_11),as.numeric(slided.vec_chr_12),as.numeric(slided.vec_chr_13),c(slided.vec))

points(dede_S23/as.numeric(density[,1]),type="l",col="#B3DE3A",ylab="",xlab="",xaxt="n",lwd=2)

# Vertical lines plotted to visualize chromosomes along the concatenated positions
abline(v=length(c(as.numeric(slided.vec_chr_1))),lty=2)
abline(v=length(c(as.numeric(slided.vec_chr_1)))+length(c(as.numeric(slided.vec_chr_2))),lty=2)
abline(v=length(c(as.numeric(slided.vec_chr_1)))+length(c(as.numeric(slided.vec_chr_2)))+length(c(as.numeric(slided.vec_chr_3))),lty=2)
abline(v=length(c(as.numeric(slided.vec_chr_1)))+length(c(as.numeric(slided.vec_chr_2)))+length(c(as.numeric(slided.vec_chr_3)))+length(c(as.numeric(slided.vec_chr_4))),lty=2)
abline(v=length(c(as.numeric(slided.vec_chr_1)))+length(c(as.numeric(slided.vec_chr_2)))+length(c(as.numeric(slided.vec_chr_3)))+length(c(as.numeric(slided.vec_chr_4)))+length(c(as.numeric(slided.vec_chr_5))),lty=2)
abline(v=length(c(as.numeric(slided.vec_chr_1)))+length(c(as.numeric(slided.vec_chr_2)))+length(c(as.numeric(slided.vec_chr_3)))+length(c(as.numeric(slided.vec_chr_4)))+length(c(as.numeric(slided.vec_chr_5)))+length(c(as.numeric(slided.vec_chr_6))),lty=2)
abline(v=length(c(as.numeric(slided.vec_chr_1)))+length(c(as.numeric(slided.vec_chr_2)))+length(c(as.numeric(slided.vec_chr_3)))+length(c(as.numeric(slided.vec_chr_4)))+length(c(as.numeric(slided.vec_chr_5)))+length(c(as.numeric(slided.vec_chr_6)))+length(c(as.numeric(slided.vec_chr_7))),lty=2)
abline(v=length(c(as.numeric(slided.vec_chr_1)))+length(c(as.numeric(slided.vec_chr_2)))+length(c(as.numeric(slided.vec_chr_3)))+length(c(as.numeric(slided.vec_chr_4)))+length(c(as.numeric(slided.vec_chr_5)))+length(c(as.numeric(slided.vec_chr_6)))+length(c(as.numeric(slided.vec_chr_7)))+length(c(as.numeric(slided.vec_chr_8))),lty=2)
abline(v=length(c(as.numeric(slided.vec_chr_1)))+length(c(as.numeric(slided.vec_chr_2)))+length(c(as.numeric(slided.vec_chr_3)))+length(c(as.numeric(slided.vec_chr_4)))+length(c(as.numeric(slided.vec_chr_5)))+length(c(as.numeric(slided.vec_chr_6)))+length(c(as.numeric(slided.vec_chr_7)))+length(c(as.numeric(slided.vec_chr_8)))+length(c(as.numeric(slided.vec_chr_9))),lty=2)
abline(v=length(c(as.numeric(slided.vec_chr_1)))+length(c(as.numeric(slided.vec_chr_2)))+length(c(as.numeric(slided.vec_chr_3)))+length(c(as.numeric(slided.vec_chr_4)))+length(c(as.numeric(slided.vec_chr_5)))+length(c(as.numeric(slided.vec_chr_6)))+length(c(as.numeric(slided.vec_chr_7)))+length(c(as.numeric(slided.vec_chr_8)))+length(c(as.numeric(slided.vec_chr_9)))+length(c(as.numeric(slided.vec_chr_10))),lty=2)
abline(v=length(c(as.numeric(slided.vec_chr_1)))+length(c(as.numeric(slided.vec_chr_2)))+length(c(as.numeric(slided.vec_chr_3)))+length(c(as.numeric(slided.vec_chr_4)))+length(c(as.numeric(slided.vec_chr_5)))+length(c(as.numeric(slided.vec_chr_6)))+length(c(as.numeric(slided.vec_chr_7)))+length(c(as.numeric(slided.vec_chr_8)))+length(c(as.numeric(slided.vec_chr_9)))+length(c(as.numeric(slided.vec_chr_10)))+length(c(as.numeric(slided.vec_chr_11))),lty=2)
abline(v=length(c(as.numeric(slided.vec_chr_1)))+length(c(as.numeric(slided.vec_chr_2)))+length(c(as.numeric(slided.vec_chr_3)))+length(c(as.numeric(slided.vec_chr_4)))+length(c(as.numeric(slided.vec_chr_5)))+length(c(as.numeric(slided.vec_chr_6)))+length(c(as.numeric(slided.vec_chr_7)))+length(c(as.numeric(slided.vec_chr_8)))+length(c(as.numeric(slided.vec_chr_9)))+length(c(as.numeric(slided.vec_chr_10)))+length(c(as.numeric(slided.vec_chr_11)))+length(c(as.numeric(slided.vec_chr_12))),lty=2)
abline(v=length(c(as.numeric(slided.vec_chr_1)))+length(c(as.numeric(slided.vec_chr_2)))+length(c(as.numeric(slided.vec_chr_3)))+length(c(as.numeric(slided.vec_chr_4)))+length(c(as.numeric(slided.vec_chr_5)))+length(c(as.numeric(slided.vec_chr_6)))+length(c(as.numeric(slided.vec_chr_7)))+length(c(as.numeric(slided.vec_chr_8)))+length(c(as.numeric(slided.vec_chr_9)))+length(c(as.numeric(slided.vec_chr_10)))+length(c(as.numeric(slided.vec_chr_11)))+length(c(as.numeric(slided.vec_chr_12)))+length(c(as.numeric(slided.vec_chr_13))),lty=2)
abline(v=5955+733/4403*805,lty=2)
abline(v=5955+1318/4403*805,lty=2)
abline(v=5955+1922/4403*805,lty=2)
abline(v=5955+2474/4403*805,lty=2)
abline(v=5955+3007/4403*805,lty=2)
abline(v=5955+3568/4403*805,lty=2)
abline(v=5955+3985/4403*805,lty=2)
abline(v=c(5955+805),lty=2)

text(c(1126, 1821, 2446, 2943, 3427, 3883, 4329, 4739, 5088, 5345, 5592, 5800, 5955,6760),y=rep(0.20,14),labels=c("6005","9856","13356","16217","19015","21672","24279","26707","28826","30487","32099","33514","34663","39064"),srt=90)

#######################
#######################
# 3. DESeq2 Model 3, data MGGP44
#######################
#######################

#################
# 3.a Model fitting, creation of output data tables containing adj. p-values and fold-changes, and extraction of relevant gene sets 
#################

# library(Hmisc)
source("C:/Users/Arthur/Desktop/Enrichment.r") # Home-made function for GO enrichment analysis (uses the "goseq" R package)

rm(Cont_all)
rm(Cont_all_l2fc)
rm(Cont_all_raw)

library(DESeq2)
Mat<-read.table("C:/Users/Arthur/Desktop/TabCount_11_juin_vrai.txt")
M=matrix(nrow=10972,ncol=40)
rownames(M)<-Mat[,1]
colnames(M)<-colnames(Mat[,2:41])
for (j in 1:40) {for (i in 1:10972){M[i,j]<-as.numeric(Mat[i,j+1])}}
dim(M)
r=unique(c(which(M[,1]>400000),which(M[,2]>400000),which(M[,3]>400000),which(M[,4]>400000),which(M[,5]>400000),which(M[,6]>400000),which(M[,7]>400000),which(M[,8]>400000),which(M[,9]>400000),which(M[,10]>400000),which(M[,11]>400000),which(M[,12]>400000),which(M[,13]>400000),which(M[,14]>400000),which(M[,15]>400000),which(M[,16]>400000),which(M[,17]>400000),which(M[,18]>400000),which(M[,19]>400000),which(M[,20]>400000),which(M[,21]>400000),which(M[,22]>400000),which(M[,23]>400000),which(M[,24]>400000),which(M[,25]>400000),which(M[,26]>400000),which(M[,27]>400000),which(M[,28]>400000),which(M[,29]>400000),which(M[,30]>400000), which(M[,31]>400000),which(M[,32]>400000),which(M[,33]>400000),which(M[,34]>400000),which(M[,35]>400000),which(M[,36]>400000),which(M[,37]>400000),which(M[,38]>400000),which(M[,39]>400000),which(M[,40]>400000)))
Juin<-read.table("C:/Users/Arthur/Desktop/data_11_juin.txt")
Outliers_chr=vector()
for (i in 1:38){Outliers_chr[i]=as.character(Juin[which(Juin[,10]==Mat[r,1][i])[1],1])}
Outliers_start=vector()
for (i in 1:38){Outliers_start[i]=as.character(Juin[which(Juin[,10]==Mat[r,1][i])[1],4])}
Outliers_end=vector()
for (i in 1:38){Outliers_end[i]=as.character(Juin[which(Juin[,10]==Mat[r,1][i])[1],5])}
Outliers_name=vector()
for (i in 1:38){Outliers_name[i]=as.character(Juin[which(Juin[,10]==Mat[r,1][i])[1],10])}
Outliers=cbind(Outliers_chr,Outliers_start,Outliers_end,Outliers_name)
M<-M[-c(which(is.element(Mat[,1],Outliers[(which(Outliers[,1]=="chr_7")),4][grep("67",Outliers[(which(Outliers[,1]=="chr_7")),4])])==TRUE)),]
dim(M)
M44<-M[,21:40] # Because the last 20 columns correspond to MGGP44 data (cols. 21:24 = Anc; cols. 25:28: St23; cols. 29:32= St17; cols. 33:40: both lineages fluct.)
coldata<-read.table("C:/Users/Arthur/Desktop/coldata.txt")
levels(coldata$Regime)<-c("Anc","Fluct","St17","St23","FluctBis")
coldata$Regime[c(17:20,37:40)]="FluctBis" # columns 33 to 36 = first lineage evolved under fluct.; cols 37:40: second lineage evolved under fluct.
coldata44<-coldata[21:40,2:3] # Because the last 20 columns correspond to MGGP44 data
levels(coldata44$Regime)
coldata44$Regime
dim(coldata44)
dim(M44)
# [1] 10950 genes x 20 MGGP44 samples (16 evolved samples + 4 samples corresponding to the data of the ancestor)
# rownames(M44): contains the 10950 gene names (annotation from Grandaubert et al., 2011)

coldata44
              # Regime Tassay
# A44171         Anc    T17
# A44172         Anc    T17
# A44234         Anc    T23
# A44235         Anc    T23
# A44323171     St23    T17
# A44323172     St23    T17
# A44323234     St23    T23
# A44323235     St23    T23
# A44117171     St17    T17
# A44117172     St17    T17
# A44117234     St17    T23
# A44117235     St17    T23
# A441F171     Fluct    T17
# A441F172     Fluct    T17
# A441F234     Fluct    T23
# A441F235     Fluct    T23
# A443F171  FluctBis    T17
# A443F172  FluctBis    T17
# A443F234  FluctBis    T23
# A443F235  FluctBis    T23

# LRT test comparing the full model on MGGP44 data versus the model with no Regime effect included
des44 <- DESeqDataSetFromMatrix(countData=M44[,1:20], colData=coldata44[1:20,], design= ~ Regime + Tassay + Regime:Tassay)
des44 <- estimateSizeFactors(des44)
des44 <- estimateDispersions(des44)
des44 <- nbinomLRT(des44,reduced= ~Tassay )
res44<-results(des44)
Cont_all_44=matrix(nrow=10950,ncol=11)
Cont_all_44_raw=matrix(nrow=10950,ncol=11)
Cont_all_44_l2fc=matrix(nrow=10950,ncol=11)
Cont_all_44[,1]<-rownames(M44)
R_effect_44<-Cont_all_44[which(as.numeric(res44[,6])<0.05),1]
length(R_effect_44)
# [1] 3659

Cont_all_44=matrix(nrow=10950,ncol=11) # Will contain FDR-corrected p-values
Cont_all_44_raw=matrix(nrow=10950,ncol=11) # Will contain raw p-values
Cont_all_44_l2fc=matrix(nrow=10950,ncol=11) # Will contain Log2 Fold-Changes

# After gene names, the first ten columns of Cont_all_44 evaluates the following contrasts between lineages
colnames(Cont_all_44)<-c("Gene_id","St23_v_St17","Fluct_v_FluctBis","Fluct_v_St17","FluctBis_v_St17","Fluct_v_St23","FluctBis_v_St23","Anc_v_St17","Anc_v_St23","Anc_v_Fluct","Anc_v_FluctBis")

Cont_all_44[,1]<-rownames(M44)
Cont_all_44_raw[,1]<-rownames(M44)
Cont_all_44_l2fc[,1]<-rownames(M44)

for (i in 1:36){
    Cont_all_44=cbind(Cont_all_44,rep(NA,10950))
    Cont_all_44_raw=cbind(Cont_all_44_raw,rep(NA,10950))
    Cont_all_44_l2fc=cbind(Cont_all_44_l2fc,rep(NA,10950))}

# The 10 main contrasts ([,2] to [,11]) comparing lineages (based on four samples for each, 2 obtained at Tassay 17°C & 2 at Tassay 23°C) with Wald tests
	
des44 <- DESeqDataSetFromMatrix(countData=M44[,1:20], colData=coldata44[1:20,], design= ~ Regime + Tassay + Regime:Tassay)
des44 <- estimateSizeFactors(des44)
des44 <- estimateDispersions(des44)
des44 <- nbinomWaldTest(des44)

resultsNames(des44)
 # [1] "Intercept"                "Regime_Fluct_vs_Anc"     
 # [3] "Regime_St17_vs_Anc"       "Regime_St23_vs_Anc"      
 # [5] "Regime_FluctBis_vs_Anc"   "Tassay_T23_vs_T17"       
 # [7] "RegimeFluct.TassayT23"    "RegimeSt17.TassayT23"    
 # [9] "RegimeSt23.TassayT23"     "RegimeFluctBis.TassayT23"
 # ==> There will be ten values to provide in combination to calculate the considered contrast (of the Wald test).
 # R has taken "Anc" as the reference of the Regime effect, "T17" as the reference of the Tassay effect

# Main effect St23 - St17
res44<-results(des44,contrast=c(0,0,-1,1,0,0,0,-0.5,0.5,0))
Cont_all_44[,2]=res44[,6]
Cont_all_44_raw[,2]=res44[,5]
Cont_all_44_l2fc[,2]=res44[,2]

# Main effect Fluct - FluctBis
res44<-results(des44,contrast=c(0,1,0,0,-1,0,0.5,0,0,-0.5))
Cont_all_44[,3]=res44[,6]
Cont_all_44_raw[,3]=res44[,5]
Cont_all_44_l2fc[,3]=res44[,2]

# Main effect Fluct - St17
res44<-results(des44,contrast=c(0,1,-1,0,0,0,0.5,-0.5,0,0))
Cont_all_44[,4]=res44[,6]
Cont_all_44_raw[,4]=res44[,5]
Cont_all_44_l2fc[,4]=res44[,2]

# Main effect FluctBis - St17
res44<-results(des44,contrast=c(0,0,-1,0,1,0,0,-0.5,0,0.5))
Cont_all_44[,5]=res44[,6]
Cont_all_44_raw[,5]=res44[,5]
Cont_all_44_l2fc[,5]=res44[,2]

# Main effect Fluct - St23
res44<-results(des44,contrast=c(0,1,0,-1,0,0,0.5,0,-0.5,0))
Cont_all_44[,6]=res44[,6]
Cont_all_44_raw[,6]=res44[,5]
Cont_all_44_l2fc[,6]=res44[,2]

# Main effect FluctBis - St23
res44<-results(des44,contrast=c(0,0,0,-1,1,0,0,0,-0.5,0.5))
Cont_all_44[,7]=res44[,6]
Cont_all_44_raw[,7]=res44[,5]
Cont_all_44_l2fc[,7]=res44[,2]

# Main effect Anc - St17
res44<-results(des44,contrast=c(0,0,-1,0,0,0,0,-0.5,0,0))
Cont_all_44[,8]=res44[,6]
Cont_all_44_raw[,8]=res44[,5]
Cont_all_44_l2fc[,8]=res44[,2]
length(which(as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),8])<0.05 & abs(as.numeric(Cont_all_44_l2fc[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),8]))>1))
# [1] 588 genes significantly DE (and with a FC >2) in Stable 17 lineage vs the MGGP44 ancestor 

DE_S17vAnc_44<-Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),1][which(as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),8])<0.05 & abs(as.numeric(Cont_all_44_l2fc[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),8]))>1)]
median(2^(abs(as.numeric(Cont_all_44_l2fc[which(is.element(Cont_all_44[,1],DE_S17vAnc_44)==TRUE),8]))))
# [1] Median FC for these 588 genes: 2.859428 (see Table 3)

# Main effect Anc - St23
res44<-results(des44,contrast=c(0,0,0,-1,0,0,0,0,-0.5,0))
Cont_all_44[,9]=res44[,6]
Cont_all_44_raw[,9]=res44[,5]
Cont_all_44_l2fc[,9]=res44[,2]
length(which(as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),9])<0.05 & abs(as.numeric(Cont_all_44_l2fc[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),9]))>1))
# [1] 532 genes significantly DE (and with a FC >2) in Stable 23 lineage vs the MGGP44 ancestor 

DE_S23vAnc_44<-Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),1][which(as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),9])<0.05 & abs(as.numeric(Cont_all_44_l2fc[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),9]))>1)]
median(2^(abs(as.numeric(Cont_all_44_l2fc[which(is.element(Cont_all_44[,1],DE_S23vAnc_44)==TRUE),9]))))
# [1] Median FC for these 532 genes: 2.841237 (see Table 3)

# Main effect Anc - Fluct
res44<-results(des44,contrast=c(0,-1,0,0,0,0,-0.5,0,0,0))
Cont_all_44[,10]=res44[,6]
Cont_all_44_raw[,10]=res44[,5]
Cont_all_44_l2fc[,10]=res44[,2]
length(which(as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),10])<0.05 & abs(as.numeric(Cont_all_44_l2fc[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),10]))>1))
# [1] 412 genes significantly DE (and with a FC >2) in Fluct_Rep1 lineage vs the MGGP44 ancestor 

DE_F1vAnc_44<-Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),1][which(as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),10])<0.05 & abs(as.numeric(Cont_all_44_l2fc[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),10]))>1)]
median(2^(abs(as.numeric(Cont_all_44_l2fc[which(is.element(Cont_all_44[,1],DE_F1vAnc_44)==TRUE),10]))))
# [1] Median FC for these 412 genes:  2.780496 (see Table 3)

# Main effect Anc - FluctBis
res44<-results(des44,contrast=c(0,0,0,0,-1,0,0,0,0,-0.5))
Cont_all_44[,11]=res44[,6]
Cont_all_44_raw[,11]=res44[,5]
Cont_all_44_l2fc[,11]=res44[,2]
length(which(as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),11])<0.05 & abs(as.numeric(Cont_all_44_l2fc[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),11]))>1))
# [1] 592 genes significantly DE (and with a FC >2) in Fluct_Rep2 lineage vs the MGGP44 ancestor 

DE_F2vAnc_44<-Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),1][which(as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),11])<0.05 & abs(as.numeric(Cont_all_44_l2fc[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),11]))>1)]
median(2^(abs(as.numeric(Cont_all_44_l2fc[which(is.element(Cont_all_44[,1],DE_F2vAnc_44)==TRUE),11]))))
# [1] Median FC for these 592 genes:  2.760954 (see Table 3)

# LRT test comparing the full model on MGGP01 data versus the model with no Tassay effect included
des44 <- DESeqDataSetFromMatrix(countData=M44[,1:20], colData=coldata44[1:20,], design= ~ Regime + Tassay + Regime:Tassay)
des44 <- estimateSizeFactors(des44)
des44 <- estimateDispersions(des44)
des44 <- nbinomLRT(des44,reduced= ~ Regime )
res44<-results(des44)
T_effect<-Cont_all_44[which(as.numeric(res44[,6])<0.05),1]
length(T_effect)
# [1] 4333

des44 <- DESeqDataSetFromMatrix(countData=M44[,1:20], colData=coldata44[1:20,], design= ~ Regime + Tassay + Regime:Tassay)
des44 <- estimateSizeFactors(des44)
des44 <- estimateDispersions(des44)
des44 <- nbinomWaldTest(des44)

# Main effect T23 - T17 averaged over the five selection regimes
res44<-results(des44,contrast=c(0,0,0,0,0,1,0.2,0.2,0.2,0.2))
Cont_all_44[,12]=res44[,6]
Cont_all_44_raw[,12]=res44[,5]
Cont_all_44_l2fc[,12]=res44[,2]
colnames(Cont_all_44)[12]<-"T23_vs_T17"

# Delta T17-T23 in the Stable_17°C lineage vs Delta T17-T23 in the Stable_23°C lineage (interaction testing)
res44<-results(des44,contrast=c(0,0,0,0,0,0,0,-1,1,0))
Cont_all_44[,13]=res44[,6]
Cont_all_44_raw[,13]=res44[,5]
Cont_all_44_l2fc[,13]=res44[,2]

# Delta T17-T23 in the Fluct lineage vs Delta T17-T23 in the FluctBis lineage (interaction testing)
res44<-results(des44,contrast=c(0,0,0,0,0,0,1,0,0,-1))
Cont_all_44[,14]=res44[,6]
Cont_all_44_raw[,14]=res44[,5]
Cont_all_44_l2fc[,14]=res44[,2]

# Delta T17-T23 in the Fluct lineage vs Delta T17-T23 in the Stable_17°C lineage (interaction testing)
res44<-results(des44,contrast=c(0,0,0,0,0,0,-1,1,0,0))
Cont_all_44[,15]=res44[,6]
Cont_all_44_raw[,15]=res44[,5]
Cont_all_44_l2fc[,15]=res44[,2]

# Delta T17-T23 in the FluctBis lineage vs Delta T17-T23 in the Stable_17°C lineage (interaction testing)
res44<-results(des44,contrast=c(0,0,0,0,0,0,0,1,0,-1))
Cont_all_44[,16]=res44[,6]
Cont_all_44_raw[,16]=res44[,5]
Cont_all_44_l2fc[,16]=res44[,2]

# Delta T17-T23 in the Fluct lineage vs Delta T17-T23 in the Stable_23°C lineage (interaction testing)
res44<-results(des44,contrast=c(0,0,0,0,0,0,-1,0,1,0))
Cont_all_44[,17]=res44[,6]
Cont_all_44_raw[,17]=res44[,5]
Cont_all_44_l2fc[,17]=res44[,2]

# Delta T17-T23 in the FluctBis lineage vs Delta T17-T23 in the Stable_23°C lineage (interaction testing)
res44<-results(des44,contrast=c(0,0,0,0,0,0,0,0,1,-1))
Cont_all_44[,18]=res44[,6]
Cont_all_44_raw[,18]=res44[,5]
Cont_all_44_l2fc[,18]=res44[,2]

# Delta T17-T23 in the ancestor vs Delta T17-T23 in the the Stable_17°C lineage (interaction testing)
res44<-results(des44,contrast=c(0,0,0,0,0,0,0,-1,0,0))
Cont_all_44[,19]=res44[,6]
Cont_all_44_raw[,19]=res44[,5]
Cont_all_44_l2fc[,19]=res44[,2]

# Delta T17-T23 in the ancestor vs Delta T17-T23 in the the Stable_23°C lineage (interaction testing)
res44<-results(des44,contrast=c(0,0,0,0,0,0,0,0,-1,0))
Cont_all_44[,20]=res44[,6]
Cont_all_44_raw[,20]=res44[,5]
Cont_all_44_l2fc[,20]=res44[,2]

# Delta T17-T23 in the ancestor vs Delta T17-T23 in the the Fluct lineage (interaction testing)
res44<-results(des44,contrast=c(0,0,0,0,0,0,-1,0,0,0))
Cont_all_44[,21]=res44[,6]
Cont_all_44_raw[,21]=res44[,5]
Cont_all_44_l2fc[,21]=res44[,2]

# Delta T17-T23 in the ancestor vs Delta T17-T23 in the the FluctBis lineage (interaction testing)
res44<-results(des44,contrast=c(0,0,0,0,0,0,0,0,0,-1))
Cont_all_44[,22]=res44[,6]
Cont_all_44_raw[,22]=res44[,5]
Cont_all_44_l2fc[,22]=res44[,2]

# Effet of the Temp. of assay for the Fluct lineage
res44<-results(des44,contrast=c(0,0,0,0,0,1,1,0,0,0))
Cont_all_44[,23]=res44[,6]
Cont_all_44_raw[,23]=res44[,5]
Cont_all_44_l2fc[,23]=res44[,2]

# Effet of the Temp. of assay for the FluctBis lineage
res44<-results(des44,contrast=c(0,0,0,0,0,1,0,0,0,1))
Cont_all_44[,24]=res44[,6]
Cont_all_44_raw[,24]=res44[,5]
Cont_all_44_l2fc[,24]=res44[,2]

# Effet of the Temp. of assay for the Stable_17°C lineage
res44<-results(des44,contrast=c(0,0,0,0,0,1,0,1,0,0))
Cont_all_44[,25]=res44[,6]
Cont_all_44_raw[,25]=res44[,5]
Cont_all_44_l2fc[,25]=res44[,2]

# Effet of the Temp. of assay for the Stable_23°C lineage
res44<-results(des44,contrast=c(0,0,0,0,0,1,0,0,1,0))
Cont_all_44[,26]=res44[,6]
Cont_all_44_raw[,26]=res44[,5]
Cont_all_44_l2fc[,26]=res44[,2]

# Effet of the Temp. of assay for the ancestor
res44<-results(des44,contrast=c(0,0,0,0,0,1,0,0,0,0))
Cont_all_44[,27]=res44[,6]
Cont_all_44_raw[,27]=res44[,5]
Cont_all_44_l2fc[,27]=res44[,2]

des44 <- DESeqDataSetFromMatrix(countData=M44[,1:20], colData=coldata44[1:20,], design= ~ Regime + Tassay + Regime:Tassay)
des44 <- estimateSizeFactors(des44)
des44 <- estimateDispersions(des44)
des44 <- nbinomLRT(des44,reduced= ~ Regime + Tassay)
res44<-results(des44)
RbyT_effect<-Cont_all_44[which(as.numeric(res44[,6])<0.05),1]
length(RbyT_effect)
# [1] 851

A_vs_F_MGGP44<-Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),][which(as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),10])<0.05 & as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),11])<0.05 & abs(as.numeric(Cont_all_44_l2fc[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),10]))>1 & abs(as.numeric(Cont_all_44_l2fc[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),11]))>1),1]
length(A_vs_F_MGGP44)
# [1] Set of 366 genes consistently DE (adjusted p<0.05 & FC >2) in both MGGP44 fluctuating lineages in comparison to the ancestor.

# Preparing data to draw the plot of the Figure 4 (MGGP44)
mt=matrix(ncol=2,nrow=length(A_vs_F_MGGP44)*4)
mt[,2]=c(-c(as.numeric(Cont_all_44_l2fc[which(is.element(Cont_all_44[,1],A_vs_F_MGGP44)==TRUE),8]),as.numeric(Cont_all_44_l2fc[which(is.element(Cont_all_44[,1],A_vs_F_MGGP44)==TRUE),9]),as.numeric(Cont_all_44_l2fc[which(is.element(Cont_all_44[,1],A_vs_F_MGGP44)==TRUE),10]),as.numeric(Cont_all_44_l2fc[which(is.element(Cont_all_44[,1],A_vs_F_MGGP44)==TRUE),11])))
mt[1:length(A_vs_F_MGGP44),1]="St17"
mt[367:732,1]="St23"
mt[733:1098,1]="Fluct"
mt[1099:1464,1]="FluctBis"
mt<-as.data.frame(mt)
mt[,2]<-as.numeric(as.character(mt[,2]))
class(mt[,1])
# [1] "factor"
class(mt[,2])
# [1] "numeric"
table(mt[,1])

library(ggplot2)
library("ggbeeswarm")
par(xpd=TRUE, mar=c(6,6,4,4))
ggplot(mt, aes(factor(V1,levels=c("St17","St23","Fluct","FluctBis")), V2)) + geom_beeswarm(size=0.5,cex=0.5,col=c(rep("#8BEAF0",366),rep("#B3DE3A",366),rep("darkviolet",366),rep("brown1",366))) + scale_fill_manual(values=c("darkviolet", "brown1", "#8BEAF0", "#B3DE3A")) + theme(legend.position="none", axis.title.x=element_blank(), axis.text=element_text(size=7,face="bold")) + labs(y="Log2 Fold-Change in comparison to MGGP44 \n \n (for DEG under fluctuations)") +scale_x_discrete(labels = c("Stable 17°C","Stable 23°C","Fluct R_1","Fluct R_2"))+ ylim(-10, 10)

# More precisely: 

length(which(as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],A_vs_F_MGGP44)==TRUE),8])<0.05))
# [1] 328 out of the 366 genes consistently DE under fluct (and with FC>2) are also DE at St17 ==> Represents 90 %, see Results in the article
# length(which(is.element(A_vs_F_MGGP44,DE_S17vAnc_44)==TRUE))
# [1] 291 out of the 366 genes are DE at St17 (here a FC>2 is imposed for St17 vs Anc)

length(which(as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],A_vs_F_MGGP44)==TRUE),9])<0.05))
# [1] 222 out of the 366 genes consistently DE under fluct (and with FC>2) are also DE at St23 ==> Represents 61 %, see Results in the Article
# length(which(is.element(A_vs_F_MGGP44,DE_S23vAnc_44)==TRUE))
# [1] 183 out of the 366 genes are DE at St23 (here a FC>2 is imposed for St23 vs Anc)

#################
# 3.b Repeatability of evolution under fluctuations: F_1 vs F_2 (Figure 3)
#################

A_vs_F_MGGP44_noFCthreshold<-R_effect_44[which(as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),10])<0.05)[which(is.element(which(as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),10])<0.05),which(as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),11])<0.05))==TRUE)]]
# ==> 1592 genes consistently DE (whatever the fold-change) in both MGGP44 fluct. lineages versus the ancestor

#tiff(file="C:/Users/ajallet/Desktop/ex_14.tiff", height=7000, width=8000, res=1000)
par(mar=c(2,4,2,0)+0.1)
plot(-c(as.numeric(Cont_all_44_l2fc[,10])),-c(as.numeric(Cont_all_44_l2fc[,11])),pch=176,col="grey50",ylab="",xlab="",las=1,cex.axis=1.30,xlim=c(-8,8),ylim=c(-8,8))
# All data (10950 genes) in grey
points(-c(as.numeric(Cont_all_44_l2fc[which(is.element(Cont_all_44[,1],A_vs_F_MGGP44_noFCthreshold)==TRUE),10])),-c(as.numeric(Cont_all_44_l2fc[which(is.element(Cont_all_44[,1],A_vs_F_MGGP44_noFCthreshold)==TRUE),11])),col="red",pch=176)
# The 1592 ones in red
abline(v=1,lty="dashed")
abline(v=c(-1),lty="dashed")
abline(h=1,lty="dashed")
abline(h=c(-1),lty="dashed")
abline(a=0,b=1)
# dev.off()

#################
# 3.c GO terms enrichment analyses
#################

Enrichment(A_vs_F_MGGP44)
go.wall[which(as.numeric(go.wall[,2])<0.05),1:2]
# category over_represented_pvalue
# 10 Defensemechanisms            1.328182e-02
# For the set of 366 genes consistently DE in both MGGP44 indep. replicates of evolution under fluct.(adjusted p<0.05 & FC >2) ==> Enrichment in "Defensemechanisms" (see Table S2)

# Genes greyed in Table S2: those specifically DE under fluctuations (i.e. only DE in fluctuating lineages BUT not DE in any of the two stable lineages) + DE in fluct. with FC>2
Specif_F_MGGP44_SDE_above2Fold<-Cont_all_44[which(is.element(Cont_all_44[,1],A_vs_F_MGGP44)==TRUE),1][which(as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],A_vs_F_MGGP44)==TRUE),8])>0.05 & as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],A_vs_F_MGGP44)==TRUE),9])>0.05)]
length(Specif_F_MGGP44_SDE_above2Fold)
#[1] 20 such genes

Enrichment(Specif_F_MGGP44_SDE_above2Fold)
go.wall[which(as.numeric(go.wall[,2])<0.05),1:2]
# category                over_represented_pvalue
# <0 lignes> 
# No enrichment for these 20 genes

# Genes specifically DE under fluctuations (i.e. only DE in fluctuating lineages BUT not DE in any of the two stable lineages) but with FC<2 in at least one fluct lineages
Specif_F_MGGP44_SDE_below2Fold<-Cont_all_44[which(as.numeric(Cont_all_44[,8])>0.05 & as.numeric(Cont_all_44[,9])>0.05),1][which(is.element(Cont_all_44[which(as.numeric(Cont_all_44[,8])>0.05 & as.numeric(Cont_all_44[,9])>0.05),1],A_vs_F_MGGP44_noFCthreshold)==TRUE)][-c(which(is.element(Cont_all_44[which(as.numeric(Cont_all_44[,8])>0.05 & as.numeric(Cont_all_44[,9])>0.05),1][which(is.element(Cont_all_44[which(as.numeric(Cont_all_44[,8])>0.05 & as.numeric(Cont_all_44[,9])>0.05),1],A_vs_F_MGGP44_noFCthreshold)==TRUE)],A_vs_F_MGGP44)==TRUE))]
length(Specif_F_MGGP44_SDE_below2Fold)
#[1] 10 such genes

# Enrichment test on the 20 + 10 (these two sets do not overlap)
Enrichment(c(Specif_F_MGGP44_SDE_above2Fold,Specif_F_MGGP44_SDE_below2Fold))
go.wall[which(as.numeric(go.wall[,2])<0.05),1:2]
# category over_represented_pvalue
# No enrichment for these 30 genes





#################
# 3.d Coexpression analysis (WGCNA clustering)
#################

library(WGCNA)

datExpr0<-t(M44)
datExpr<-datExpr0

# We indicate the Regime of selection and the Temperaturure of assay that correspond to each of the 20 samples
traitData<-matrix(nrow=20,ncol=2)
traitData[1:4,1]=as.character(as.factor("Anc"))
traitData[5:8,1]=as.character(as.factor("St23"))
traitData[9:12,1]=as.character(as.factor("St17"))
traitData[13:16,1]=as.character(as.factor("Fluct"))
traitData[17:20,1]=as.character(as.factor("FluctBis"))
traitData[,2]=rep(c(rep(c(as.character(as.factor("T17"))),2),rep(c(as.character(as.factor("T23"))),2)))
datTraits<-traitData

# Creation of a vector containing the genes for which one of the 21 constrasts (wether main or interaction effects) - and the corresponding LRT test for the tested covariate -
# are significant
Complete<-rownames(des44)[which(is.element(Cont_all_44[,1],unique(c(Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),1][which(as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),2])<0.05 | as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),3])<0.05 | as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),4])<0.05 | as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),5])<0.05 | as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),6])<0.05 | as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),7])<0.05 | as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),8])<0.05 | as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),9])<0.05 | as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),10])<0.05 | as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),11])<0.05)],Cont_all_44[which(is.element(Cont_all_44[,1],T_effect)==TRUE),1][which(as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],T_effect)==TRUE),12])<0.05)],Cont_all_44[which(is.element(Cont_all_44[,1],RbyT_effect)==TRUE),1][which(as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],RbyT_effect)==TRUE),13])<0.05 | as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],RbyT_effect)==TRUE),14])<0.05 | as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],RbyT_effect)==TRUE),15])<0.05 | as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],RbyT_effect)==TRUE),16])<0.05 | as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],RbyT_effect)==TRUE),17])<0.05 | as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],RbyT_effect)==TRUE),18])<0.05 | as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],RbyT_effect)==TRUE),19])<0.05 | as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],RbyT_effect)==TRUE),20])<0.05 | as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],RbyT_effect)==TRUE),21])<0.05 | as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],RbyT_effect)==TRUE),22])<0.05)])))==TRUE)]
length(Complete)
# [1] 5401

datExpr<-datExpr[,which(is.element(colnames(datExpr),Complete)==TRUE)]
dim(datExpr)
# [1]   20 5401

# Normalization step
des44 <- DESeqDataSetFromMatrix(countData=t(datExpr[,]), colData=coldata44[1:20,], design= ~ Regime + Tassay + Regime:Tassay)
des44 <- estimateSizeFactors(des44)
datExpr<-t(counts(des44, normalized=TRUE))

# Log 2 transformation of the data (before that, zeroes has to be removed)				
v=vector()
for (i in 1:20) {v=c(v,as.numeric(which(datExpr[i,]==0)))}
v<-unique(v)
datExpr<-datExpr[,-c(v)]
datExpr<-log2(datExpr)
dim(datExpr)
#   20 5305
# ==> Co-expression clustering will be performed by WGCNA on this dataset (i.e. on log2 values of the normalized counts of the 4947 genes contained in datExpr)

# gsg = goodSamplesGenes(datExpr, verbose = 3);
# gsg$allOK
# [1] TRUE

library(preprocessCore)
library(impute)
library(WGCNA)

# Choose a set of soft-thresholding powers
powers = c(1:20)
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
    main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
	 
-sign(sft$fitIndices[,3])*sft$fitIndices[,2]
# [1] -0.590615938 -0.217541707  0.006177922  0.236966207  0.468329864  0.600386513  0.683215098
# [8]  0.732977032  0.757480212  0.786062774  0.804502051  0.803400420  0.810694500  0.805143682
# [15]  0.817505574  0.829891148  0.829319506  0.840294354  0.827566407  0.836354923

# For a soft thresholding power of 13, R²=0.81 for scale-free topology test
softPower = 13;
adjacency = adjacency(datExpr, power = softPower,type="signed");

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
rm(TOM)

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");

# We set the minimum module size at 40:
minModuleSize = 40

# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
deepSplit = 2, pamRespectsDendro = FALSE,
minClusterSize = minModuleSize);

length(dynamicMods)
#[1] 5305

# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)

# Check here: if TRUE, it works
# dynamicColors==MEList$validColors
# unique(dynamicColors)==unique(MEList$validColors)

MEs = MEList$eigengenes
# View(MEs)
# One column = expression profile across the different experimental combinations of {Regime x Temp. of asaay} treatements dispalyed for the summary (i.e. the eignegene) of the genes 
# contained in the module of this column

# Clustering of module eigengenes:
MEDiss = 1-cor(MEs);

# Module merging based on a threshold of similmarity (here 75%)
MEDissThres = 0.25
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

# ==> After merging, 15 MGGP44 modules grouping genes according to their coexpression profiles across the experimental {Regime x Temp. of asaay} conditions
# With WGCNA, modules are automatically named according to their colours

names(mergedMEs)
# For each of the 15 eigengene, the name of the one of the module "colourX" is "MEcolourX"

for (i in 1:15) {print(paste(unique(mergedColors)[i],length(which(mergedColors==unique(mergedColors)[i]))),sep=" ")}
# [1] "grey60 1509"
# [1] "lightcyan 1455"
# [1] "skyblue 93"
# [1] "darkmagenta 181"
# [1] "salmon 313"
# [1] "cyan 224"
# [1] "darkolivegreen 185"
# [1] "tan 160"
# [1] "violet 76"
# [1] "black 409"
# [1] "sienna3 66"
# [1] "green 304"
# [1] "orange 101"
# [1] "darkorange 144"
# [1] "steelblue 85"
# The sum is equal to 5305 ==> Every genes fall into a module and genes are uniquely assigned to exactly one coexpression module

length(which(is.element(A_vs_F_MGGP01_noFCthreshold,c(rownames(t(datExpr))[which(mergedColors[]==unique(mergedColors)[1])],rownames(t(datExpr))[which(mergedColors[]==unique(mergedColors)[2])],rownames(t(datExpr))[which(mergedColors[]==unique(mergedColors)[6])]))==TRUE))
# ==> The three MGGP44 "grey60", "lightcyan", and "cyan" modules contain ~90% of the genes that evolved under Fluctuations

# Within each module, genes are sorted by their correlation to the eigengene: the more a gene looks like the average (i.e. to the eigengene) the more it is a hub within the module.
# Hubs are genes that are connected to many other genes of their module and that probably drive the "coexpression behaviour" across conditions of all genes of the module
# //
# When ranked according to their correlation to the eigengene (= a measure of their centrality), are genes involved annotated as "Transcription factors"
#  or as involved in "Signal transduction" more central ? 

################
# 3.e Analysis of expression plasticity in evolved lineages in regard to the ancestor
################

Plastic_MGGP44_Anc<-rownames(M44)[which(is.element(rownames(M44),Cont_all_44[which(as.numeric(Cont_all_44[,27])<0.05 & abs(as.numeric(Cont_all_44_l2fc[,27]))>1),1])==TRUE)]
# [1] 699 genes displaying plasticity in gene expression (at a two-fold threshold) in the MGGP44 ancestor
Robust_MGGP44_Anc<-rownames(M44)[which(is.element(rownames(M44),Plastic_MGGP44_Anc)==FALSE)]

Plastic_MGGP44_S17<-rownames(M44)[which(is.element(rownames(M44),Cont_all_44[which(as.numeric(Cont_all_44[,25])<0.05 & abs(as.numeric(Cont_all_44_l2fc[,25]))>1),1])==TRUE)]
# [1] 892 genes displaying plasticity in gene expression (at a two-fold threshold) in MGGP44 Stable_17°C lineage
Robust_MGGP44_S17<-rownames(M44)[which(is.element(rownames(M44),Plastic_MGGP44_S17)==FALSE)]

Plastic_MGGP44_S23<-rownames(M44)[which(is.element(rownames(M44),Cont_all_44[which(as.numeric(Cont_all_44[,26])<0.05 & abs(as.numeric(Cont_all_44_l2fc[,26]))>1),1])==TRUE)]
# [1] 900 genes displaying plasticity in gene expression (at a two-fold threshold) in MGGP44 Stable_23°C lineage
Robust_MGGP44_S23<-rownames(M44)[which(is.element(rownames(M44),Plastic_MGGP44_S23)==FALSE)]

Plastic_MGGP44_F1<-rownames(M44)[which(is.element(rownames(M44),Cont_all_44[which(as.numeric(Cont_all_44[,23])<0.05 & abs(as.numeric(Cont_all_44_l2fc[,23]))>1),1])==TRUE)]
# [1] 574 genes displaying plasticity in gene expression (at a two-fold threshold) in MGGP44 Fluct_1 lineage
Robust_MGGP44_F1<-rownames(M44)[which(is.element(rownames(M44),Plastic_MGGP44_F1)==FALSE)]

Plastic_MGGP44_F2<-rownames(M44)[which(is.element(rownames(M44),Cont_all_44[which(as.numeric(Cont_all_44[,24])<0.05 & abs(as.numeric(Cont_all_44_l2fc[,24]))>1),1])==TRUE)]
# [1] 549 genes displaying plasticity in gene expression (at a two-fold threshold) in MGGP44 Fluct_2 lineage
Robust_MGGP44_F2<-rownames(M44)[which(is.element(rownames(M44),Plastic_MGGP44_F2)==FALSE)]

# Genes that lost plasticty in expression during evolution (see Figure 6 - Panel C and Table S9)
length(which(is.element(Plastic_MGGP44_Anc,Robust_MGGP44_S17)==TRUE))
# 553 such genes after evolution under the Stable_17°C regime
length(which(is.element(Plastic_MGGP44_Anc,Robust_MGGP44_S23)==TRUE))
# 504 such genes after evolution under the Stable_23°C regime
length(which(is.element(Plastic_MGGP44_Anc,Robust_MGGP44_F1)==TRUE))
# 563 such genes after evolution under the Fluctuating regime (rep. F1)
length(which(is.element(Plastic_MGGP44_Anc,Robust_MGGP44_F2)==TRUE))
# 592 such genes after evolution under the Fluctuating regime (rep. F2) /!\

# Genes which expression became plastic during evolution (see Figure 6 Panel C and Table S9)
length(which(is.element(Robust_MGGP44_Anc,Plastic_MGGP44_S17)==TRUE))
# 746 such genes after evolution under the Stable_17°C regime
length(which(is.element(Robust_MGGP44_Anc,Plastic_MGGP44_S23)==TRUE))
# 705 such genes after evolution under the Stable_23°C regime
length(which(is.element(Robust_MGGP44_Anc,Plastic_MGGP44_F1)==TRUE))
# 438 such genes after evolution under the Fluctuating regime (rep. F1)
length(which(is.element(Robust_MGGP44_Anc,Plastic_MGGP44_F2)==TRUE))
# 442 such genes after evolution under the Fluctuating regime (rep. F2)

# Genes which expression remained plastic (i.e. those displaying plasticity in the ancestor and still in the considered evolved lineage)(see Figure 6 Panel C and Table S9)
length(which(is.element(Plastic_MGGP44_Anc,Plastic_MGGP44_S17)==TRUE))
# 146 such genes after evolution under the Stable_17°C regime
length(which(is.element(Plastic_MGGP44_Anc,Plastic_MGGP44_S23)==TRUE))
# 195 such genes after evolution under the Stable_23°C regime
length(which(is.element(Plastic_MGGP44_Anc,Plastic_MGGP44_F1)==TRUE))
# 136 such genes after evolution under the Fluctuating regime (rep. F1)
length(which(is.element(Plastic_MGGP44_Anc,Plastic_MGGP44_F2)==TRUE))
# 107 such genes after evolution under the Fluctuating regime (rep. F2)

###############
# 3.f Density of DE genes - normalized by gene density - along the genome for each of the four lineages of MGGP44 background (Figure 5)
###############

# Creation of the vector "density" that contains the density of annotated genes per 100 kb along the genome, with concatenated positions (from chromosome_1 to chromosome_21)
# Sliding window approach

name_bis<-read.table("C://Users/Arthur/Desktop/name_bis.txt",header=TRUE)
data0<-name_bis[,c(5,3,2)]
data<-data0[,c(3:2)]

X=vector()
for (i in 1:dim(data)[1]){X[i]=strsplit(as.character((data[,1])),split="chr_")[[i]][2]}
data<-cbind(as.numeric(X),data[,2])
data<-cbind(data,rep(NA,dim(data)[1]))

Y<-data

# Creation of data1, containing concatenated positions
for (i in 1:21){
    if (i==1){
        tmp=Y[Y[,1]==i,]
        tmp=tmp[order(tmp[,2]),]
        data1=tmp
        pos1=tmp[,2]
    }
    
    if (i>1){
        tmp=Y[Y[,1]==i,]
        tmp=tmp[order(tmp[,2]),]
        data1<-rbind(data1, tmp)
        pos=tmp[,2]
        max.value<-max(pos1)
        pos=pos+max.value
        pos1<-c(pos1, pos)
        data1[,3]<-pos1
    }        
}

dataA=vector()
for (i in 1:21)
{dataA[i]=tail(data1[which(data1[,1]==i),2])[6]}
# dataA

for (i in 1:21)

{chr<-unique(data0[,3])[i]
list.bin=seq(1000, tail(data1[which(data1[,1]==i),2])[6], by=1000)

tab=cbind(list.bin, rep(NA, length(list.bin)))
colnames(tab)=c("list.bin","count")

# will fill up col 2 with the number of events  
# use data1: concat positions of the events

k=1
v5=as.numeric(data1[which(data1[,1]==i),2])
for (i in list.bin){
    if(i==list.bin[1]){
        bin=v5[v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    if(i!=list.bin[1]){
        bin=v5[v5>lim & v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    
    k=k+1   
}

print(sum(tab[,2]))

require(zoo)
vec <- zoo(tab[,2])
sum(vec)
# [1] 1983

length(vec)
# [1] 6005

assign(paste("vec",chr,sep="_"), vec)

slided.vec<- rollapply(vec, width = 379, by = 5, FUN = mean)
assign(paste("slided.vec",chr,sep="_"), slided.vec)}

# Below: for the mini-chromosomes specifically (from chromosome_14 to chromosome_21)
# Creation of only one contiguous vector for all these eight chromosomes (necessary due to their small sizes and the chosen large window)

rm(X)
rm(data)
rm(data0)
rm(data1)

data0<-name_bis[which(name_bis[,2]=="chr_14" | name_bis[,2]=="chr_15" | name_bis[,2]=="chr_16" | name_bis[,2]=="chr_17" | name_bis[,2]=="chr_18" | name_bis[,2]=="chr_19" | name_bis[,2]=="chr_20" | name_bis[,2]=="chr_21"),c(5,3,2)]
data<-data0[,c(3:2)]

X=vector()
for (i in 1:dim(data)[1]){X[i]=strsplit(as.character((data[,1])),split="chr_")[[i]][2]}
data<-cbind(as.numeric(X),data[,2])
data<-cbind(data,rep(NA,dim(data)[1]))

Y<-data

# Creation of data1, containing concatenated positions
for (i in 14:21){
    if (i==14){
        tmp=Y[Y[,1]==i,]
        tmp=tmp[order(tmp[,2]),]
        data1=tmp
        pos1=tmp[,2]
    }
    
    if (i>14){
        tmp=Y[Y[,1]==i,]
        tmp=tmp[order(tmp[,2]),]
        data1<-rbind(data1, tmp)
        pos=tmp[,2]
        max.value<-max(pos1)
        pos=pos+max.value
        pos1<-c(pos1, pos)
        data1[,3]<-pos1
    }        
}

NN=tail(data1[,3])[6]

list.bin=seq(1000, tail(data1[,3])[6], by=1000)
tab=cbind(list.bin, rep(NA, length(list.bin)))
colnames(tab)=c("list.bin","count")

# will fill up col 2 with the number of events  
# use data1: concat positions of the events

k=1
v5=as.numeric(data1[,3])
for (i in list.bin){
    if(i==list.bin[1]){
        bin=v5[v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    if(i!=list.bin[1]){
        bin=v5[v5>lim & v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    
    k=k+1   
}

vec <- zoo(tab[,2])

slided.vec<- rollapply(vec, width = 379, by = 5, FUN = mean)

density=c(as.numeric(slided.vec_chr_1),as.numeric(slided.vec_chr_2),as.numeric(slided.vec_chr_3),as.numeric(slided.vec_chr_4),as.numeric(slided.vec_chr_5),as.numeric(slided.vec_chr_6),as.numeric(slided.vec_chr_7),as.numeric(slided.vec_chr_8),as.numeric(slided.vec_chr_9),as.numeric(slided.vec_chr_10),as.numeric(slided.vec_chr_11),as.numeric(slided.vec_chr_12),as.numeric(slided.vec_chr_13),as.numeric(slided.vec))
# This vector contains the density (per 100 kb) of annotated genes along the genome, with concatenated positions (from chromosome_1 to chromosome_21)

# Density of DE genes along the genome for the lineage MGGP44 F1 (normalized by the density of annotated genes)

NNN<-data1
end_density<-NN
rm(NN)
NN<-NNN

name_bis<-read.table("C://Users/Arthur/Desktop/name_bis.txt",header=TRUE)
name_bis<-name_bis[which(is.element(name_bis[,5],Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),1][which(as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),10])<0.05 & abs(as.numeric(Cont_all_44_l2fc[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),10]))>1)])==TRUE),]
# Les 412 DE chez MGGP44 F_1
dim(name_bis)
# [1] 412   7

new_NN<-NN[which(is.element(NN[,2],name_bis[which(is.element(name_bis[,3],NN[,2])==TRUE),3])==TRUE),]

data<-name_bis[,c(2,3)]

X=vector()
for (i in 1:dim(data)[1]){X[i]=strsplit(as.character((data[,1])),split="chr_")[[i]][2]}
data<-cbind(as.numeric(X),data[,2])
data<-cbind(data,rep(NA,dim(data)[1]))

Y<-data

for (i in 1:13){
    if (i==1){
        tmp=Y[Y[,1]==i,]
        tmp=tmp[order(tmp[,2]),]
        data1=tmp
        pos1=tmp[,2]
    }
    
    if (i>1){
        tmp=Y[Y[,1]==i,]
        tmp=tmp[order(tmp[,2]),]
        data1<-rbind(data1, tmp)
        pos=tmp[,2]
        max.value<-max(pos1)
        pos=pos+max.value
        pos1<-c(pos1, pos)
        data1[,3]<-pos1
    }        
}
# Contient les positions concaténées tout bien comme y faut pour les 584 (666-82) gènes des CC

for (i in 1:13)
    
{chr<-paste("chr",unique(data1[,1])[i],sep="_")
list.bin=seq(1000, dataA[i], by=1000)

tab=cbind(list.bin, rep(NA, length(list.bin)))
colnames(tab)=c("list.bin","count")

# will fill up col 2 with the number of events  
# use data1: concat positions of the events

k=1
v5=as.numeric(data1[which(data1[,1]==i),2])
for (i in list.bin){
    if(i==list.bin[1]){
        bin=v5[v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    if(i!=list.bin[1]){
        bin=v5[v5>lim & v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    
    k=k+1   
}

print(sum(tab[,2]))

require(zoo)
vec <- zoo(tab[,2])
sum(vec)
# [1] 1983

length(vec)
# [1] 6005

assign(paste("vec",chr,sep="_"), vec)

slided.vec<- rollapply(vec, width = 379, by = 5, FUN = mean)
assign(paste("slided.vec",chr,sep="_"), slided.vec)}

list.bin=seq(1000, end_density, by=1000)
tab=cbind(list.bin, rep(NA, length(list.bin)))
colnames(tab)=c("list.bin","count")

k=1
v5=as.numeric(new_NN[,3])
for (i in list.bin){
    if(i==list.bin[1]){
        bin=v5[v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    if(i!=list.bin[1]){
        bin=v5[v5>lim & v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    
    k=k+1   
}

print(sum(tab[,2]))
# 650

require(zoo)
vec <- zoo(tab[,2])
sum(vec)
# 650

length(vec)
# [1] 4403

slided.vec<- rollapply(vec, width = 379, by = 5, FUN = mean)

dede_F1<-c(as.numeric(slided.vec_chr_1),as.numeric(slided.vec_chr_2),as.numeric(slided.vec_chr_3),as.numeric(slided.vec_chr_4),as.numeric(slided.vec_chr_5),as.numeric(slided.vec_chr_6),as.numeric(slided.vec_chr_7),as.numeric(slided.vec_chr_8),as.numeric(slided.vec_chr_9),as.numeric(slided.vec_chr_10),as.numeric(slided.vec_chr_11),as.numeric(slided.vec_chr_12),as.numeric(slided.vec_chr_13),c(slided.vec))

density<-read.table("C:/Users/Arthur/Desktop/Density.txt")
plot(dede_F1/as.numeric(density[,1]),type="l",col="darkviolet",ylab="",xlab="",xaxt="n", lwd=2)

# Density of DE genes along the genome for the lineage MGGP44 F2 (normalized by the density of annotated genes)

name_bis<-read.table("C://Users/Arthur/Desktop/name_bis.txt",header=TRUE)
name_bis<-name_bis[which(is.element(name_bis[,5],Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),1][which(as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),11])<0.05 & abs(as.numeric(Cont_all_44_l2fc[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),11]))>1)])==TRUE),]
# Les 592 DE chez MGGP44 F_2
dim(name_bis)
# [1] 592   7

new_NN<-NN[which(is.element(NN[,2],name_bis[which(is.element(name_bis[,3],NN[,2])==TRUE),3])==TRUE),]

data<-name_bis[,c(2,3)]

X=vector()
for (i in 1:dim(data)[1]){X[i]=strsplit(as.character((data[,1])),split="chr_")[[i]][2]}
data<-cbind(as.numeric(X),data[,2])
data<-cbind(data,rep(NA,dim(data)[1]))

Y<-data

for (i in 1:13){
    if (i==1){
        tmp=Y[Y[,1]==i,]
        tmp=tmp[order(tmp[,2]),]
        data1=tmp
        pos1=tmp[,2]
    }
    
    if (i>1){
        tmp=Y[Y[,1]==i,]
        tmp=tmp[order(tmp[,2]),]
        data1<-rbind(data1, tmp)
        pos=tmp[,2]
        max.value<-max(pos1)
        pos=pos+max.value
        pos1<-c(pos1, pos)
        data1[,3]<-pos1
    }        
}

for (i in 1:13)
    
{chr<-paste("chr",unique(data1[,1])[i],sep="_")
list.bin=seq(1000, dataA[i], by=1000)

tab=cbind(list.bin, rep(NA, length(list.bin)))
colnames(tab)=c("list.bin","count")

# will fill up col 2 with the number of events  
# use data1: concat positions of the events

k=1
v5=as.numeric(data1[which(data1[,1]==i),2]) 
for (i in list.bin){
    if(i==list.bin[1]){
        bin=v5[v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    if(i!=list.bin[1]){
        bin=v5[v5>lim & v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    
    k=k+1   
}

print(sum(tab[,2]))

require(zoo)
vec <- zoo(tab[,2])
sum(vec)
# [1] 1983

length(vec)
# [1] 6005

assign(paste("vec",chr,sep="_"), vec)

slided.vec<- rollapply(vec, width = 379, by = 5, FUN = mean)
assign(paste("slided.vec",chr,sep="_"), slided.vec)}

list.bin=seq(1000, end_density, by=1000)
tab=cbind(list.bin, rep(NA, length(list.bin)))
colnames(tab)=c("list.bin","count")

k=1
v5=as.numeric(new_NN[,3])
for (i in list.bin){
    if(i==list.bin[1]){
        bin=v5[v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    if(i!=list.bin[1]){
        bin=v5[v5>lim & v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    
    k=k+1   
}

print(sum(tab[,2]))
# 650

require(zoo)
vec <- zoo(tab[,2])
sum(vec)
# 650

length(vec)
# [1] 4403

slided.vec<- rollapply(vec, width = 379, by = 5, FUN = mean)

dede_F2<-c(as.numeric(slided.vec_chr_1),as.numeric(slided.vec_chr_2),as.numeric(slided.vec_chr_3),as.numeric(slided.vec_chr_4),as.numeric(slided.vec_chr_5),as.numeric(slided.vec_chr_6),as.numeric(slided.vec_chr_7),as.numeric(slided.vec_chr_8),as.numeric(slided.vec_chr_9),as.numeric(slided.vec_chr_10),as.numeric(slided.vec_chr_11),as.numeric(slided.vec_chr_12),as.numeric(slided.vec_chr_13),c(slided.vec))

points(dede_F2/as.numeric(density[,1]),type="l",col="brown1",ylab="",xlab="",xaxt="n",lwd=2)

# Density of DE genes along the genome for the lineage MGGP44 Stable_17°C (normalized by the density of annotated genes)

name_bis<-read.table("C://Users/Arthur/Desktop/name_bis.txt",header=TRUE)
name_bis<-name_bis[which(is.element(name_bis[,5],Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),1][which(as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),8])<0.05 & abs(as.numeric(Cont_all_44_l2fc[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),8]))>1)])==TRUE),]
# Les 588 DE chez MGGP44 Stable_17°C
dim(name_bis)
# [1] 588   7

new_NN<-NN[which(is.element(NN[,2],name_bis[which(is.element(name_bis[,3],NN[,2])==TRUE),3])==TRUE),]

data<-name_bis[,c(2,3)]

X=vector()
for (i in 1:dim(data)[1]){X[i]=strsplit(as.character((data[,1])),split="chr_")[[i]][2]}
data<-cbind(as.numeric(X),data[,2])
data<-cbind(data,rep(NA,dim(data)[1]))

Y<-data

for (i in 1:13){
    if (i==1){
        tmp=Y[Y[,1]==i,]
        tmp=tmp[order(tmp[,2]),]
        data1=tmp
        pos1=tmp[,2]
    }
    
    if (i>1){
        tmp=Y[Y[,1]==i,]
        tmp=tmp[order(tmp[,2]),]
        data1<-rbind(data1, tmp)
        pos=tmp[,2]
        max.value<-max(pos1)
        pos=pos+max.value
        pos1<-c(pos1, pos)
        data1[,3]<-pos1
    }        
}

for (i in 1:13)
    
{chr<-paste("chr",unique(data1[,1])[i],sep="_")
list.bin=seq(1000, dataA[i], by=1000)

tab=cbind(list.bin, rep(NA, length(list.bin)))
colnames(tab)=c("list.bin","count")

# will fill up col 2 with the number of events  
# use data1: concat positions of the events

k=1
v5=as.numeric(data1[which(data1[,1]==i),2])
for (i in list.bin){
    if(i==list.bin[1]){
        bin=v5[v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    if(i!=list.bin[1]){
        bin=v5[v5>lim & v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    
    k=k+1   
}

print(sum(tab[,2]))

require(zoo)
vec <- zoo(tab[,2])
sum(vec)
# [1] 1983

length(vec)
# [1] 6005

assign(paste("vec",chr,sep="_"), vec)

slided.vec<- rollapply(vec, width = 379, by = 5, FUN = mean)
assign(paste("slided.vec",chr,sep="_"), slided.vec)}

list.bin=seq(1000, end_density, by=1000)
tab=cbind(list.bin, rep(NA, length(list.bin)))
colnames(tab)=c("list.bin","count")

k=1
v5=as.numeric(new_NN[,3])
for (i in list.bin){
    if(i==list.bin[1]){
        bin=v5[v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    if(i!=list.bin[1]){
        bin=v5[v5>lim & v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    
    k=k+1   
}

print(sum(tab[,2]))
# 650

require(zoo)
vec <- zoo(tab[,2])
sum(vec)
# 650

length(vec)
# [1] 4403

slided.vec<- rollapply(vec, width = 379, by = 5, FUN = mean)

dede_S17<-c(as.numeric(slided.vec_chr_1),as.numeric(slided.vec_chr_2),as.numeric(slided.vec_chr_3),as.numeric(slided.vec_chr_4),as.numeric(slided.vec_chr_5),as.numeric(slided.vec_chr_6),as.numeric(slided.vec_chr_7),as.numeric(slided.vec_chr_8),as.numeric(slided.vec_chr_9),as.numeric(slided.vec_chr_10),as.numeric(slided.vec_chr_11),as.numeric(slided.vec_chr_12),as.numeric(slided.vec_chr_13),c(slided.vec))

points(dede_S17/as.numeric(density[,1]),type="l",col="#8BEAF0",ylab="",xlab="",xaxt="n",lwd=2)

# Density of DE genes along the genome for the lineage MGGP44 Stable_23°C (normalized by the density of annotated genes)

name_bis<-read.table("C://Users/Arthur/Desktop/name_bis.txt",header=TRUE)
name_bis<-name_bis[which(is.element(name_bis[,5],Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),1][which(as.numeric(Cont_all_44[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),9])<0.05 & abs(as.numeric(Cont_all_44_l2fc[which(is.element(Cont_all_44[,1],R_effect_44)==TRUE),9]))>1)])==TRUE),]
# Les 532 DE chez MGGP44 Stable_23°C
dim(name_bis)
# [1] 532   7

new_NN<-NN[which(is.element(NN[,2],name_bis[which(is.element(name_bis[,3],NN[,2])==TRUE),3])==TRUE),]

data<-name_bis[,c(2,3)]

X=vector()
for (i in 1:dim(data)[1]){X[i]=strsplit(as.character((data[,1])),split="chr_")[[i]][2]}
data<-cbind(as.numeric(X),data[,2])
data<-cbind(data,rep(NA,dim(data)[1]))

Y<-data

for (i in 1:13){
    if (i==1){
        tmp=Y[Y[,1]==i,]
        tmp=tmp[order(tmp[,2]),]
        data1=tmp
        pos1=tmp[,2]
    }
    
    if (i>1){
        tmp=Y[Y[,1]==i,]
        tmp=tmp[order(tmp[,2]),]
        data1<-rbind(data1, tmp)
        pos=tmp[,2]
        max.value<-max(pos1)
        pos=pos+max.value
        pos1<-c(pos1, pos)
        data1[,3]<-pos1
    }        
}

for (i in 1:13)
    
{chr<-paste("chr",unique(data1[,1])[i],sep="_")
list.bin=seq(1000, dataA[i], by=1000)

tab=cbind(list.bin, rep(NA, length(list.bin)))
colnames(tab)=c("list.bin","count")

# will fill up col 2 with the number of events  
# use data1: concat positions of the events

k=1
v5=as.numeric(data1[which(data1[,1]==i),2])
for (i in list.bin){
    if(i==list.bin[1]){
        bin=v5[v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    if(i!=list.bin[1]){
        bin=v5[v5>lim & v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    
    k=k+1   
}

print(sum(tab[,2]))

require(zoo)
vec <- zoo(tab[,2])
sum(vec)
# [1] 1983

length(vec)
# [1] 6005

assign(paste("vec",chr,sep="_"), vec)

slided.vec<- rollapply(vec, width = 379, by = 5, FUN = mean)
assign(paste("slided.vec",chr,sep="_"), slided.vec)}

list.bin=seq(1000, end_density, by=1000)
tab=cbind(list.bin, rep(NA, length(list.bin)))
colnames(tab)=c("list.bin","count")

k=1
v5=as.numeric(new_NN[,3])
for (i in list.bin){
    if(i==list.bin[1]){
        bin=v5[v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    if(i!=list.bin[1]){
        bin=v5[v5>lim & v5<i]
        if ((length(bin)>0)==FALSE){count=0}
        if (length(bin)>0) {count=length(bin)}
        tab[k,2]<-count
        lim=i
    }
    
    k=k+1   
}

print(sum(tab[,2]))
# 650

require(zoo)
vec <- zoo(tab[,2])
sum(vec)
# 650

length(vec)
# [1] 4403

slided.vec<- rollapply(vec, width = 379, by = 5, FUN = mean)

dede_S23<-c(as.numeric(slided.vec_chr_1),as.numeric(slided.vec_chr_2),as.numeric(slided.vec_chr_3),as.numeric(slided.vec_chr_4),as.numeric(slided.vec_chr_5),as.numeric(slided.vec_chr_6),as.numeric(slided.vec_chr_7),as.numeric(slided.vec_chr_8),as.numeric(slided.vec_chr_9),as.numeric(slided.vec_chr_10),as.numeric(slided.vec_chr_11),as.numeric(slided.vec_chr_12),as.numeric(slided.vec_chr_13),c(slided.vec))

points(dede_S23/as.numeric(density[,1]),type="l",col="#B3DE3A",ylab="",xlab="",xaxt="n",lwd=2)

# Vertical lines plotted to visualize chromosomes along the concatenated positions
abline(v=length(c(as.numeric(slided.vec_chr_1))),lty=2)
abline(v=length(c(as.numeric(slided.vec_chr_1)))+length(c(as.numeric(slided.vec_chr_2))),lty=2)
abline(v=length(c(as.numeric(slided.vec_chr_1)))+length(c(as.numeric(slided.vec_chr_2)))+length(c(as.numeric(slided.vec_chr_3))),lty=2)
abline(v=length(c(as.numeric(slided.vec_chr_1)))+length(c(as.numeric(slided.vec_chr_2)))+length(c(as.numeric(slided.vec_chr_3)))+length(c(as.numeric(slided.vec_chr_4))),lty=2)
abline(v=length(c(as.numeric(slided.vec_chr_1)))+length(c(as.numeric(slided.vec_chr_2)))+length(c(as.numeric(slided.vec_chr_3)))+length(c(as.numeric(slided.vec_chr_4)))+length(c(as.numeric(slided.vec_chr_5))),lty=2)
abline(v=length(c(as.numeric(slided.vec_chr_1)))+length(c(as.numeric(slided.vec_chr_2)))+length(c(as.numeric(slided.vec_chr_3)))+length(c(as.numeric(slided.vec_chr_4)))+length(c(as.numeric(slided.vec_chr_5)))+length(c(as.numeric(slided.vec_chr_6))),lty=2)
abline(v=length(c(as.numeric(slided.vec_chr_1)))+length(c(as.numeric(slided.vec_chr_2)))+length(c(as.numeric(slided.vec_chr_3)))+length(c(as.numeric(slided.vec_chr_4)))+length(c(as.numeric(slided.vec_chr_5)))+length(c(as.numeric(slided.vec_chr_6)))+length(c(as.numeric(slided.vec_chr_7))),lty=2)
abline(v=length(c(as.numeric(slided.vec_chr_1)))+length(c(as.numeric(slided.vec_chr_2)))+length(c(as.numeric(slided.vec_chr_3)))+length(c(as.numeric(slided.vec_chr_4)))+length(c(as.numeric(slided.vec_chr_5)))+length(c(as.numeric(slided.vec_chr_6)))+length(c(as.numeric(slided.vec_chr_7)))+length(c(as.numeric(slided.vec_chr_8))),lty=2)
abline(v=length(c(as.numeric(slided.vec_chr_1)))+length(c(as.numeric(slided.vec_chr_2)))+length(c(as.numeric(slided.vec_chr_3)))+length(c(as.numeric(slided.vec_chr_4)))+length(c(as.numeric(slided.vec_chr_5)))+length(c(as.numeric(slided.vec_chr_6)))+length(c(as.numeric(slided.vec_chr_7)))+length(c(as.numeric(slided.vec_chr_8)))+length(c(as.numeric(slided.vec_chr_9))),lty=2)
abline(v=length(c(as.numeric(slided.vec_chr_1)))+length(c(as.numeric(slided.vec_chr_2)))+length(c(as.numeric(slided.vec_chr_3)))+length(c(as.numeric(slided.vec_chr_4)))+length(c(as.numeric(slided.vec_chr_5)))+length(c(as.numeric(slided.vec_chr_6)))+length(c(as.numeric(slided.vec_chr_7)))+length(c(as.numeric(slided.vec_chr_8)))+length(c(as.numeric(slided.vec_chr_9)))+length(c(as.numeric(slided.vec_chr_10))),lty=2)
abline(v=length(c(as.numeric(slided.vec_chr_1)))+length(c(as.numeric(slided.vec_chr_2)))+length(c(as.numeric(slided.vec_chr_3)))+length(c(as.numeric(slided.vec_chr_4)))+length(c(as.numeric(slided.vec_chr_5)))+length(c(as.numeric(slided.vec_chr_6)))+length(c(as.numeric(slided.vec_chr_7)))+length(c(as.numeric(slided.vec_chr_8)))+length(c(as.numeric(slided.vec_chr_9)))+length(c(as.numeric(slided.vec_chr_10)))+length(c(as.numeric(slided.vec_chr_11))),lty=2)
abline(v=length(c(as.numeric(slided.vec_chr_1)))+length(c(as.numeric(slided.vec_chr_2)))+length(c(as.numeric(slided.vec_chr_3)))+length(c(as.numeric(slided.vec_chr_4)))+length(c(as.numeric(slided.vec_chr_5)))+length(c(as.numeric(slided.vec_chr_6)))+length(c(as.numeric(slided.vec_chr_7)))+length(c(as.numeric(slided.vec_chr_8)))+length(c(as.numeric(slided.vec_chr_9)))+length(c(as.numeric(slided.vec_chr_10)))+length(c(as.numeric(slided.vec_chr_11)))+length(c(as.numeric(slided.vec_chr_12))),lty=2)
abline(v=length(c(as.numeric(slided.vec_chr_1)))+length(c(as.numeric(slided.vec_chr_2)))+length(c(as.numeric(slided.vec_chr_3)))+length(c(as.numeric(slided.vec_chr_4)))+length(c(as.numeric(slided.vec_chr_5)))+length(c(as.numeric(slided.vec_chr_6)))+length(c(as.numeric(slided.vec_chr_7)))+length(c(as.numeric(slided.vec_chr_8)))+length(c(as.numeric(slided.vec_chr_9)))+length(c(as.numeric(slided.vec_chr_10)))+length(c(as.numeric(slided.vec_chr_11)))+length(c(as.numeric(slided.vec_chr_12)))+length(c(as.numeric(slided.vec_chr_13))),lty=2)
abline(v=5955+733/4403*805,lty=2)
abline(v=5955+1318/4403*805,lty=2)
abline(v=5955+1922/4403*805,lty=2)
abline(v=5955+2474/4403*805,lty=2)
abline(v=5955+3007/4403*805,lty=2)
abline(v=5955+3568/4403*805,lty=2)
abline(v=5955+3985/4403*805,lty=2)
abline(v=c(5955+805),lty=2)

text(c(1126, 1821, 2446, 2943, 3427, 3883, 4329, 4739, 5088, 5345, 5592, 5800, 5955,6760),y=rep(0.20,14),labels=c("6005","9856","13356","16217","19015","21672","24279","26707","28826","30487","32099","33514","34663","39064"),srt=90)
