############################################
##################### R script supplemental File XX Article. Jallet et al.
##################### Goal: evidence variant frequency differences between temperatures for Plastic Transcripts
##################### 
##################### last version January 2020

### Contents

1- Calculate allele frequency and Confidence Intervals within each RNAseq sample with variable sites and multiple alleles
2- Compare Replicates



	# 1- Allele Frequency and Confidence Interval using binomiale

# estimation of read number col 6 7 >> replicate 1
# estimation of read number col 12 13 >> replicate 2

N=c('Sample_01_17', 'Sample_01_23')

N=c('Sample_44_17', 'Sample_44_23')

## call function
freqCI(N)



##  function to calc freq and CI (binomiale)
freqCI=function(N){}
    for (n in N){

name.input=gsub("()", "",paste(n, "_2alleles_AD_cov10_dup.txt",sep='') )
ad=read.table(file=name.input, header=T)

xx<-matrix(NA, nrow=nrow(ad), ncol=6)

colnames(xx)<- c("freq1","low_CI1", "high_CI1","freq2", "low_CI2", "high_CI2")

# calc allele freq and Confidence Intervals
# allele freq are reference allele frequency

    for (i in 1: nrow(ad)){
        
ci1<-binom.test(ad[i,6], (ad[i,6]+ad[i,7]))
f1<-ad[i,6]/(ad[i,6]+ad[i,7])
xx[i,1]<-f1
xx[i,2:3]<-ci1[[4]]

ci2<-binom.test(ad[i,12], (ad[i,12]+ad[i,13]))
f2<-ad[i,12]/(ad[i,12]+ad[i,13])
xx[i,4]<-f2
xx[i,5:6]<-ci2[[4]]

}
                            
   AD=cbind (ad, xx)
   
   name.output=gsub("()", "",paste(n, ".freq.txt",sep='') )
   write.table(AD, file=name.output)
   
 }
 } #end function
 
 
  

 # reopen both files:

list.name=c('Sample_01_17', 'Sample_01_23')
list.name=c('Sample_44_17', 'Sample_44_23')

for (i in list.name){
name.input=gsub("()", "",paste(i, ".freq.txt",sep='') )
 
AD = read.table(name.input, header=TRUE)

chr=unique(AD[,'chrom.x'])
### Need a loop per chromosome
         for (j in chr){

         
            X=AD[AD[,'chrom.x']==j,]
            L=1:nrow(X)
            X=cbind(L, X)
            # remove site positions with high frequency where non should occur (in ancestors)

        ### FILTER 1 : 'het like' sites in ancestors need to be filtered out
        
            high<- X[X[,'freq1']>0.04 & X[,'freq2'] >0.04,]
            check=high$pos #OK
            pos<-high$L
                
                newX = X[-pos,]
                mis1<-X[pos,]
                
                if (j==chr[1]){filtX=newX}
                if (j!=chr[1]){filtX<-rbind(filtX, newX)}
                
                if (j==chr[1]){M1=mis1}
                if (j!=chr[1]){M1<-rbind(M1, mis1)}
        
        ### FILTER 2: when there is an overlap between the two calls:  max(A) < min(B) || min(A) < max(B) ; need to filter out the opposite pattern 
        # L is ready
                
            dis= filtX[filtX[,'high_CI1']< filtX[,'low_CI2'] & filtX[,'high_CI2']< filtX[,'low_CI1'],]


            check=dis$pos #OK
            pos2<-dis$L
                
            
            if(length(pos2)==0){
                                    endX=filtX
                                    mis1=mis1
                                    }
                k=1
                if(length(pos2)>0){
                
                
                                    endX = filtX[-pos2,]
                                    mis2<-filtX[pos2,]
                k=k+1                    
                                    }
                if (k==1){mis=mis1}
                if (k>1){mis=rbind(mis1, mis2)}
            
                    
                if (j==chr[1]){
                                    filteredX=endX
                                    
                                    MM <-mis
                                    
                                }
                if (j!=chr[1]){
                                    filteredX<-rbind(filteredX, endX)
                                    MM <-rbind(MM, mis)
                                }
        
        
  }# end chr
  name.out=gsub("()", "",paste(i, ".filtered_freq.txt",sep='') )
  
  write.table(filteredX, file=name.out)
  
  name.rm_out=gsub("()", "",paste(i, "_remove_freq.txt",sep='') )
  
  write.table(MM, file=name.rm_out)
  
  
  
    }#end genotypes
    



   


## NOW open the data file for the evolved lineages:
# do it separately for MGGP01 and MGGP44: 

rm(i)
names=c("Sample_12F_17", "Sample_12F_23","Sample_13F_17", "Sample_13F_23","Sample_1123_17", "Sample_1123_23", "Sample_1217_17", "Sample_1217_23")
 
names=c("Sample_441F_17", "Sample_441F_23","Sample_443F_17", "Sample_443F_23","Sample_44323_17", "Sample_44323_23", "Sample_44117_17", "Sample_44117_23")

 
 
 ############# identify the positions of non plastic genes, to be removed for this manuscript, apply filter to all evolved lineages 
 
 
 m1=read.table(file="Sample_01_17.MGGP01_remove_freq.txt", header=T)
 m2=read.table(file="Sample_01_23.MGGP01_remove_freq.txt", header=T)

m1=read.table(file="Sample_44_17_remove_freq.txt", header=T)
 m2=read.table(file="Sample_44_23_remove_freq.txt", header=T)


# remove pos of m1 and m2 from all samples
 
    for (i in names){


#open files and remove known "non plastic" variant positions on each chrom from the ancestor RNAseq data/SNP calling

    name.in=gsub("()", "",paste(i,"_2alleles_AD_cov10_dup.txt",sep='') )

M <- read.table(name.in, header=T)


chrom=unique(M[,'chrom.x'])
for (j in chr){

	tmp=M[M[,'chrom.x']==j,]

		mm1=m1[m1[,'chrom.x']==j,]
			pos1=mm1[,'pos']
		mm2=m2[m2[,'chrom.x']==j,]
			pos2=mm2[,'pos']
		tmp1=tmp[tmp[,'pos']!=pos1,]

 			tmp2=tmp1[tmp1[,'pos']!=pos2,] 
 			
 			if (j == chrom[1]){file=tmp2}
 			if (j != chrom[1]){file=rbind(file, tmp2)}
 			

  


}
 name.out=gsub("()", "",paste(i,"_2alleles_AD_cov10_dup_FILTERED.txt",sep='') )

file=cbind(1:nrow(file), file) # add a vector to avoid problem when I reopen the files
 file=na.omit(file)
 
    write.table(file, file=name.out)


}

	# 2- compare the snp freq CI and calling between temperatures for each evolved lineage :


# 2-1- calc CI in filtered files

rm(i)

N=c("Sample_12F_17", "Sample_12F_23","Sample_13F_17", "Sample_13F_23","Sample_1123_17", "Sample_1123_23", "Sample_1217_17", "Sample_1217_23")
N=c("Sample_441F_17", "Sample_441F_23","Sample_443F_17", "Sample_443F_23","Sample_44323_17", "Sample_44323_23", "Sample_44117_17", "Sample_44117_23")



    for (n in N){

name.input=gsub("()", "",paste(n, "_2alleles_AD_cov10_dup_FILTERED.txt",sep='') )
ad=read.table(file=name.input, header=T)

xx<-matrix(NA, nrow=nrow(ad), ncol=6)

colnames(xx)<- c("freq1","low_CI1", "high_CI1","freq2", "low_CI2", "high_CI2")

# calc allele freq and Confidence Intervals
# allele freq are reference allele frequency

    for (i in 1: nrow(ad)){
        
ci1<-binom.test(ad[i,7], (ad[i,7]+ad[i,8]))
f1<-ad[i,7]/(ad[i,7]+ad[i,8])
xx[i,1]<-f1
xx[i,2:3]<-ci1[[4]]

ci2<-binom.test(ad[i,13], (ad[i,13]+ad[i,14]))
f2<-ad[i,13]/(ad[i,13]+ad[i,14])
xx[i,4]<-f2
xx[i,5:6]<-ci2[[4]]

}
                            
   AD=cbind (ad, xx)
   
   name.output=gsub("()", "",paste(n, ".freq.txt",sep='') )
   write.table(AD, file=name.output)
   
 }
 
 
 

# 2-2-compare between temperatures : 

	# call function 'check_calling_temp' using pairs of samples = same genotype at 2 temperatures
	
	
check_calling_temp(c("Sample_12F_17", "Sample_12F_23"))	
check_calling_temp(c("Sample_13F_17", "Sample_13F_23"))	
check_calling_temp(c("Sample_1123_17", "Sample_1123_23"))		
check_calling_temp(c("Sample_1217_17", "Sample_1217_23"))	

check_calling_temp(c("Sample_441F_17", "Sample_441F_23"))	
check_calling_temp(c("Sample_443F_17", "Sample_443F_23"))	
check_calling_temp(c("Sample_44323_17", "Sample_44323_23"))		
check_calling_temp(c("Sample_44117_17", "Sample_44117_23"))	

	
list.name=	c("Sample_44117_17", "Sample_44117_23")
	#### function ####
check_calling_temp<-function(list.name){

for (i in list.name){

cat(i, "\n")
name.input=gsub("()", "",paste(i, ".freq.txt",sep='') )
 
if (i==list.name[1]) {AD17 = read.table(name.input, header=TRUE)}
if (i==list.name[2]) {AD23 = read.table(name.input, header=TRUE)}

		}
		
		
		### keep only the common positions between temperatures
		
		
		
chr=unique(AD17[,'chrom.x'])
### Need a loop per chromosome
         for (j in chr){

         
            X=AD17[AD17[,'chrom.x']==j,]
            L=1:nrow(X)
            X=cbind(L, X)

	   		Y=AD23[AD23[,'chrom.x']==j,]
            L=1:nrow(Y)
            Y=cbind(L, Y)

        ### FILTER 1 : keep common positions only, of course (where the transcripts are expressed in both sides)
        
        XY = merge(X, Y, by.x='pos', by.y='pos')

       
### FILTER 2 : search for absence of overlap between the two temperature :  max(A) < min(B) || min(A) < max(B) 
### and both repeats, 

                
            D= XY[min(XY[,'low_CI1.x'], XY[,'low_CI2.x']) > max(XY[,'high_CI1.y'], XY[,'high_CI2.y']),]
            
            
            Dp= XY[min(XY[,'low_CI1.y'], XY[,'low_CI2.y']) > max(XY[,'high_CI1.x'], XY[,'high_CI2.x']),]
       
                       
                if (j==chr[1]){
                                  XY1=XY  
                                    differ1=rbind(D, Dp)
                                    
                                }
                if (j==chr[1] & nrow(differ1)==0){
                                    cat(i, "\n",j, 'no difference', "\n")
                                    
                                }
                if (j!=chr[1]){XY1<-rbind(XY1, XY)}                                
                if (j!=chr[1] & nrow(differ1)==0){
                                    cat(j, 'no difference', "\n")
                                   
                                }
                                
                if (j!=chr[1] & nrow(differ1)>0){
                                    differ1<-rbind(differ1, differ)
                                   
                                }
        
        	### VERIFY ALLELES ARE THE SAME BETWEEN THE TWO TEMPERATURES
        	
        	tt=cbind(as.character(XY[,'alt.x.x']),as.character(XY[,'alt.y.y']))
        	tt=cbind(data.frame(tt), tt[,1]==tt[,2])
        	allele=tt[tt[,3]=='FALSE',]
        	
        	      if(dim(allele)[1]==0) 	{cat("same alleles, ok", "\n\n")}
        	
        	
        	}# end chr
  
   name.output=gsub("()", "",paste(i, "_check_alleles.txt",sep='') )
   write.table(XY1, file=name.output)
 
  
  
}### end-of-function
	
 
	
	# merge these positions when in plastic genes
























   





















