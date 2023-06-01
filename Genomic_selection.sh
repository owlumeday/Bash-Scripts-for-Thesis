#!/bin/bash


cd  ~/thesis/GWAS_DATA

mkdir Genomic_Selection
cd Genomic_Selection
#______________________excluding Families with less than 10 sibs__________________________#

cat ~/thesis/GWAS_DATA/Mowi_Lice_2017.txt | awk 'NR > 1 {print}' > Mowi_Lice_2017_noheader.txt

module load R/4.2.1-foss-2022a
echo 'Mowi <- read.table("Mowi_Lice_2017_noheader.txt", header = F)
#____________________Data Preparation_______________________________#
Mowi <- Mowi[,1:7]
colnames(Mowi) <- c("Individual","Sire", "Dam" ,"Counter","Lice_Count","Body_Weight", "Tank")

#I only need animals with phenotypes so i exclude those without phenotype
Mowi <- Mowi[complete.cases(Mowi$Lice_Count, Mowi$Body_Weight),] #2935 animals made it through
Mowi$Family <- paste0(Mowi$Sire,"_",Mowi$Dam)
Mowi$Counter_x_Tank <- paste0(Mowi$Counter,"_",Mowi$Tank)
#unique(Mowi$Family)

frq <- data.frame(table(Mowi$Family)) #Determining the amount of sibs per family
frq <- frq[order(frq$Freq),]

remove <- which(frq$Freq < 10) #checking if there are families with less than 10 siblings
frq <- frq[-c(remove),] #removing families with less than 10 siblings

Mowi <- subset(Mowi, Family %in% frq$Var1)
#2878 individuals met this req but i need a number divisible by 5 for my 5 fold
#therefore i expunge 3 individuals at random
 
a <- sample(1:nrow(Mowi),3, replace=F) #I need a number divisible by 5 so i remove 3 random obs
Mowi <- Mowi[-a,] #Now i have 2875 individuals which is divisible into 5 groups
write.table(Mowi[,c("Individual","Individual")], "individual.txt", col.names = F, row.names = F,quote = F)
# I will need this individual file created in order to extract the GEno of this 2875 individuals


#___________________________________Splitting to 5-Fold Validation set_____________________#
Mowi <- Mowi[order(Mowi$Family),] #ordering by family so as to divide into folds
Mowi$Fold <- as.factor(c(1:5)) #making the n-fold based on family
Mowi <- Mowi[order(Mowi$Individual),] #Ordering back by individual so as to match with geno file
Pheno <-Mowi
rownames(Pheno)<-Mowi[,1]
write.table(Pheno, "Pheno.txt", col.names = T, row.names = F,quote = F)' > Individuals_and_Pheno.R

Rscript Individuals_and_Pheno.R

#Extracting the Families with 10sibs or more from the IMPUTED data and making it a genotype format

#______________________________GRM ARRAY DATA_____________________________________
#NOTE: i noticed array snps had missing values.. this will cause a problem in GRM.. Therefore i extract array data from the imputed data
#since they were not imputed for.. apart for the missings which were imputed for

cat ~/thesis/Array_MLA/ArrayNew/ArrayCorrect.bim | cut -f 2 > Array_SNPs.txt 
#extracting the individuals with families more than 10sibs and full genotype for array

../plink --bfile ~/thesis/GWAS_DATA/FinalOut1_imp_ACGT_Filt --make-rel square \
--keep individual.txt --extract Array_SNPs.txt --hwe 1e-25 \
--chr-set 31 --out FinalOut1_Array_ACGT_FiltHet

#______________________________GRM IMPUTED DATA_____________________________________


../plink --bfile ~/thesis/GWAS_DATA/FinalOut1_imp_ACGT_Filt --make-rel square \
--nonfounders --allow-no-sex \
--keep individual.txt --hwe 1e-25 \
--chr-set 31 --out FinalOut1_imp_ACGT_FiltHet

### R

#_________________________________R SCRIPT GENOMIC PREDICTION____________________

echo '
#________________________________PHENOTYPE AND OTHER DATA_____________________

Pheno <- read.table("Pheno.txt", header = T)

#I need to adjust for Fixed effects
model <- lm(log1p(Pheno$Lice_Count) ~ Pheno$Body_Weight + Pheno$Counter_x_Tank)
Pheno$Adjusted_pheno <- model$residuals
#The adjusted phenotype here is the residual of linear model corrected for fixed effects.
#Therefore the adjusted pheno consist of genetic and unexplained variance.
#_____________________________GRM data_______________________________________

##Array GRM
grm.Array <- read.table('FinalOut1_Array_ACGT_FiltHet.rel',stringsAsFactors = FALSE,header=FALSE)
animID <-read.table('FinalOut1_Array_ACGT_FiltHet.rel.id',header=FALSE,stringsAsFactors=FALSE)
colnames(grm.Array)<-animID$V2
rownames(grm.Array)<-animID$V2
grm.Array <- as.matrix(grm.Array)

grm.Imp <- read.table('FinalOut1_imp_ACGT_FiltHet.rel',stringsAsFactors = FALSE,header=FALSE)
animID <-read.table('FinalOut1_imp_ACGT_FiltHet.rel.id',header=FALSE,stringsAsFactors=FALSE)
colnames(grm.Imp)<-animID$V2
rownames(grm.Imp)<-animID$V2
grm.Imp <- as.matrix(grm.Imp)
#write.table(grm.Imp,'grm.Imp.txt',quote=F,col.names=T,row.names=T,sep='\t')

#_________________5-fold Cross Validation_______________________#

n_fold <- 5 #number of folds.. since i am doing a 5-fold cross validation
Family_freq <- matrix(nrow=length(unique(Pheno$Family)), ncol = n_fold)
rownames(Family_freq) <- unique(Pheno$Family)
colnames(Family_freq) <- paste0("Freq_",1:ncol(Family_freq))#This tells the number of individuals in a family belonging to each fold


test_individuals <- matrix(nrow=nrow(Pheno)/n_fold, ncol = n_fold)
colnames(test_individuals) <- paste0("Fold_",1:ncol(test_individuals))

for (i in 1:n_fold){
  #i =1
  fold <- subset(Pheno, Fold == i)
  test_individuals[,i] <- rownames(fold) #stores the rowname of the test_individuals
  tmp_freq <- data.frame(table(fold$Family))
  Family_freq[,i] <- tmp_freq$Freq
}
Family_freq<-as.data.frame(Family_freq)
Family_freq$No_of_Sibs <- rowSums(Family_freq) #Having an idea of the sum of individuals per family

#write.table(Family_freq, "Family_freq.txt", col.names = T, row.names = T,quote = F)

test_individuals <- as.data.frame(test_individuals)

#________________________________BGLR___________________________________________#

#install.packages("BGLR")
#install.packages("BGLR")
library(BGLR)

fold = 5

ETA_Array <- list(FIXED = list(~factor(Counter_x_Tank)+Body_Weight,
                         data=Pheno,model="FIXED"), #fixed effect of counter*tank and covariate BW
                G=list(K=grm.Array, model = "RKHS"))#Marker effect
			
ETA_Imputed <- list(FIXED = list(~factor(Counter_x_Tank)+Body_Weight,
                         data=Pheno,model="FIXED"), #fixed effect of counter*tank and covariate BW
				G=list(K=grm.Imp, model = "RKHS"))#Marker effect
			
uhat_Array <- matrix(ncol = fold, nrow = nrow(test_individuals)) #The estimated breeding values
Array_Accuracy <- matrix(nrow=1,ncol=fold)

uhat_Imputed <- matrix(ncol = fold, nrow = nrow(test_individuals)) #The estimated breeding values
Imputed_Accuracy <- matrix(nrow=1,ncol=fold)

y_adjusted <- matrix(ncol = fold, nrow = nrow(test_individuals)) #The y_adjusted is same for array and imputed

for (i in 1:fold){
  y <- as.matrix(log1p(Pheno$Lice_Count))
  ywithNA <- y
  tst <- as.numeric(test_individuals[,i])
  ywithNA[tst] <- NA
  
 #_______________________FITTING MODEL__________________________#
  GBLUP_Array <- BGLR(y = ywithNA, ETA=ETA_Array, nIter=1500, burnIn=100)
  uhat_Array[,i] <- GBLUP_Array$ETA$G$u[tst]
  y_adjusted[,i] <- Pheno$Adjusted_pheno[tst] #The y_adjusted is same for array and imputed
  Array_Accuracy[1,i] <- cor(GBLUP_Array$ETA$G$u[tst],Pheno$Adjusted_pheno[tst])/sqrt(0.2059) 
  #we check the correlation between the EBV and phenotype and divide by 0.2059 which is the heritability using array data
  
  GBLUP_Imputed <- BGLR(y = ywithNA, ETA=ETA_Imputed, nIter=1500, burnIn=100)
  #_____________________________ACCURACY_______________________#
  uhat_Imputed[,i] <- GBLUP_Imputed$ETA$G$u[tst]
  Imputed_Accuracy[1,i] <- cor(GBLUP_Imputed$ETA$G$u[tst],Pheno$Adjusted_pheno[tst])/sqrt(0.2181)  
  #we check the correlation between the EBV and phenotype and divide by 0.2181 which is the heritability using imputed data
}

#_________________________________ALL FOLD ACCURACY_________________________
#The idea here is to merge together all test individuals and compare their predicted bv with true phenotype

all_uhat_Array <- tidyr::gather(as.data.frame(uhat_Array))  #This is the predicted bv for test individuals across all validations for array and imputed
all_uhat_Imputed <- tidyr::gather(as.data.frame(uhat_Imputed))
all_yadj <- tidyr::gather(as.data.frame(y_adjusted))

All_fold-accuracy_Array <- cor(all_uhat_Array[,2],all_yadj[,2])/sqrt(0.2059) 
All_fold-accuracy_Imputed <- cor(all_uhat_Imputed[,2],all_yadj[,2])/sqrt(0.2181) 

#Standard error of correlation
sqrt((1-cor(all_uhat_Array[,2],all_yadj[,2])^2)/(2875-2))  #Standard error of correlation for array

sqrt((1-cor(all_uhat_Imputed[,2],all_yadj[,2])^2)/(2875-2)) 
##################
#(GBLUP$ETA$MRK$R2)#Explained variance of marker' > BGLR_Script_Array.R
#  #plot(GBLUP$yHat[tst],y[tst])
Rscript BGLR_Script_Array.R

