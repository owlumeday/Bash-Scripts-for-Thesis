#!/bin/bash

echo ""
#-----------------------SCRIPT FOR IMPUTATION ACCURACY--------------------------#

#First its necessary to carry out a minor imputation... 
#To fill up the missing values in the reference data

echo "Since the whole genome 4million SNPs for Salmon is spread out on 29 chromosomes, \
which chromosome will you like to impute for?"

#------------------------INPUT FILES------------------------------------------#

chromosome=$1  #This should be a number between 1 to 29
#-----------------------------------------------------------------------------#

################################################
# Create temporary folder for analysis
FOLDER=Imputation_Accuracy_chr_${chromosome}
mkdir ${FOLDER}
cp FImpute3 plink IndividualSampling.R calc.R updatedFAM_SampleID.fam ArrayDataFinal* ${FOLDER}/.
cd ${FOLDER}
#################################################

#The IndividualSampling.R is a pre-writen R script that splits the individuals into batches of test and validation
#The updatedFAM_SampleID.fam file has the long name and order of the individuals.


echo "-------------------------CONVERTING VCF TO BINARY & RECODING------------------------------"

./plink --vcf ../FiltSNPs_${chromosome}.vcf \
--make-bed --double-id \
--out ref_FiltSNPs_${chromosome} \
--autosome-num ${chromosome} \
--recode A

#This step here is to change the names in FAM file to the long names.

cp updatedFAM_SampleID.fam ref_FiltSNPs_${chromosome}.fam

./plink --bfile ref_FiltSNPs_${chromosome} \
--make-bed --double-id \
--out ref_FiltSNPs_${chromosome} \
--autosome-num ${chromosome} --recode A

rm *.log
rm *.nosex

echo "---------------------------FILE CONVERSION DONE--------------------------------"
echo ""


echo "------------------------FILE MANIPULATION FOR IMPUTATION-------------------------"

awk 'NR>1 {print $2,1}' ref_FiltSNPs_${chromosome}.raw > IDs.txt
awk 'NR>1' ref_FiltSNPs_${chromosome}.raw | cut -d' ' -f7- | awk '{gsub(/NA/,5); print}' | 
awk 'BEGIN {FS=" ";OFS=""} {$1=$1; print}' | 
paste -d' ' IDs.txt - > Fimp.geno
echo 'IID Chip Call.........' > header
cat header Fimp.geno > ref.geno_chr_${chromosome}

cat ref_FiltSNPs_${chromosome}.bim | awk '{print $2,$1,$4,NR}' > tmp
echo 'SNP_ID Chr Pos Chip1' > chipheader
cat chipheader tmp > ref.snpinfo_chr_${chromosome}

rm IDs.txt Fimp.geno chipheader header tmp *.raw

echo ""
echo "------------------------------DONE WITH FILE MANUPULATION------------------------"
echo ""
echo ""


echo "----------------------------PREPARING FIMPUTE FILE FOR MINOR/FIRST IMPUTATION-------------------------------"
echo ""

#The aim of this imputation is to fill up the missing genotypes in the sequence data.

echo "title=*population based imputation*;
genotype_file=*ref.geno_chr_${chromosome}*;
snp_info_file=*ref.snpinfo_chr_${chromosome}*;
output_folder=*Minor_imputation${chromosome}*;
save_genotype;
njob=5;" > minor_imputation.ctr
sed -i 's/*/"/g' minor_imputation.ctr

echo " "
echo "data processing eneded, FImpute will start soon ............"
echo " "

#***** run FImpute *****#
./FImpute3 minor_imputation.ctr

echo " "
echo " "
echo "**********************************************************"
echo "******           Imputation finished             *********"
echo "******                                           *********"
echo "**********************************************************"


#---------------------- Extract the imputed data and make a PLINK file--------------------------#

cp ./Minor_imputation${chromosome}/genotypes_imp.txt .

cat genotypes_imp.txt | awk 'NR>1 {print $3}' | 
awk 'BEGIN {FS="";OFS=" "} {$1=$1; print $0}' | 
awk '{for (i=1;i<=NF;i++) { if($i==0) $i="1 1"; else if($i==1) $i="1 2"; else if($i==2) $i="2 2"; else if($i==5) $i="0 0"}  print}' > geno
cat genotypes_imp.txt | awk 'NR>1 {print $1,$1,0,0,0,-9}' > ids.txt
paste -d' ' ids.txt geno > file.ped
cat ref.snpinfo_chr_${chromosome} | awk 'NR>1 {print $2,$1,0,$3}' > file.map
rm ids.txt geno


#The file ped and map file is the complete reference genotype for all the individuals

./plink --file file --make-bed --autosome-num ${chromosome} --out file 

#This is the complete plink format of reference population(sequence data) - 4million SNPs

rm genotypes_imp.txt 

#---------------------------------SECOND STAGE ---------------------------------#

echo " "

cut -f 1-2 -d " " ref_FiltSNPs*.fam > All_individuals.txt

module load R/4.1.2-foss-2021b

Rscript IndividualSampling.R 

#Rscript that generates the test and training individuals systematically.

rm IndividualSampling.R All_individuals.txt


#creating a file of the SNPS on array for the particular chromosome

./plink --bfile ArrayDataFinal --make-bed --chr-set 29 --chr ${chromosome} --allow-extra-chr --out ArrayDataFinal_chr_${chromosome}

cut -f 2 ArrayDataFinal_chr_${chromosome}.bim > Array_SNPs.txt #Array SNPs (part of 50k on this particular chromosome) to be extracted for test individuals

#---------------------Extracting the Training and Test individual dataset---------------------#

icnt=0 #iteration counter
for i in {1..10}
do
let "icnt=icnt+1"
echo "${icnt} Test_Individuals_${i}.txt"

./plink --bfile file  --remove Test_Individuals_${i}.txt --make-bed --recode A --autosome-num ${chromosome} --out TrainingSet_${i}

echo " "
#we keep 10sets of systematically selected 177 individuals for the training by removing 20 test_individuals

./plink --bfile file \
--keep Test_Individuals_${i}.txt \
--make-bed \
--autosome-num ${chromosome} \
--out Temp_TestSet_${i}
#we keep 10sets of 20 systematically selected individuals for the test

echo " "

./plink --bfile Temp_TestSet_${i} \
--extract Array_SNPs.txt --make-bed \
--recode A --autosome-num ${chromosome} \
--out TestSet_${i}

cat TestSet_${i}.bim | cut -f 1-4 > TestSet_${i}.map #This are the Low Density SNPs

echo""
echo " "
done
rm Temp_TestSet*.bim Test_Individuals* *.log *.nosex 

#------------------------THIRD PART-------------------------#

echo "-------------------------Reference/Training Genotype----------------------------"

icnt=0 #iteration counter
for i in {1..10}
do
let "icnt=icnt+1"
echo "${icnt} TrainingSet_${i}.raw"

awk 'NR>1 {print $2,1}' TrainingSet_${i}.raw > IDs_TrainingSet_${i}.txt
awk 'NR>1' TrainingSet_${i}.raw | cut -d ' ' -f 7- | awk '{gsub(/NA/,5); print}' |
awk 'BEGIN {FS=" ";OFS=""} {$1=$1; print}' |
paste -d' ' IDs_TrainingSet_${i}.txt - > Fimp.geno_TrainingSet_${i}
echo 'IID Chip Call.........' > header

cat header Fimp.geno_TrainingSet_${i} > geno_TrainingSet_${i}
echo ""
echo ""


echo "-----------------START REFERNCE/TRAININGSET SNP_INFORMATION-----------------------"
echo ""
echo "${icnt} TrainingSet_${i}.bim"
echo ""

cat TrainingSet_${i}.bim | awk '{print $2,$1,$4,NR}' > tmp_${i}
echo 'SNP_ID Chr Pos Chip1' > chipheader
cat chipheader tmp_${i} > snpinfo_TrainingSet_${i}


echo ""
echo ""
echo "------------------DONE WITH REFERENCE/TRAININGSET SNP_INFO---------------------------"
done
rm IDs*.txt header Fimp* ref* chipheader tmp*
wc -l snpinfo_TrainingSet_*



###________________________________________________________________________________#

echo "-------------------------VALIDATION/TEST Genotype----------------------------"

icnt=0 #iteration counter
for i in {1..10}
do
let "icnt=icnt+1"
echo "${icnt} TestSet_${i}.raw"

awk 'NR>1 {print $2,2}' TestSet_${i}.raw > IDs_TestSet_${i}.txt
awk 'NR>1' TestSet_${i}.raw | cut -d ' ' -f 7- | awk '{gsub(/NA/,5); print}' |
awk 'BEGIN {FS=" ";OFS=""} {$1=$1; print}' |
paste -d' ' IDs_TestSet_${i}.txt - > Fimp.geno_TestSet_${i}
echo 'IID Chip Call.........' > header
cat header Fimp.geno_TestSet_${i} > geno_TestSet_${i}
echo ""
echo ""


echo "-----------------START VALIDATION/TESTSET SNP_INFORMATION-----------------------"
echo "${icnt} TestSet_${i}.bim"
echo ""

cat TestSet_${i}.bim | awk '{print $2,$1,$4,NR}' > tmp_${i}
echo 'SNP_ID Chr Pos Chip1' > chipheader
cat chipheader tmp_${i} > snpinfo_TestSet_${i}

echo ""
echo ""
echo "------------------DONE WITH REFERENCE/TRAININGSET SNP_INFO---------------------------"
done
rm IDs*.txt header Fimp* chipheader tmp* *.raw
wc -l snpinfo_TestSet_*

#-------------------------------FOURTH PART----------------------------------#

echo "*********** Create a single MAPfile for FIMPUTE **********"

icnt=0
for i in {1..10}
do
let "icnt=icnt+1"
echo "${icnt} snpinfo_TestSet_${i} and snpinfo_TrainingSet_${i}"
echo ""

cp snpinfo_TestSet_${i} tmpval_${i} 
cp snpinfo_TrainingSet_${i} tmpref_${i}
done



#---------------------------Rscript to merge the SNPS-----------------------------#

echo "for (i in 1:10){
refmap <- read.table(paste0('tmpref_',i),header=T)
valmap <- read.table(paste0('tmpval_',i),header=T)
refvalmap <- merge(refmap,valmap,by=1,all=T,sort=F)
refvalmap <- refvalmap[order(refvalmap[,2],refvalmap[,3]),]
refvalmap[,5] <- 0
refvalmap[which(refvalmap[,7]!=0),5] <-seq(1,length(which(refvalmap[,7]!=0)),1)
refvalmap <- refvalmap[,c(1,2,3,4,5)]
colnames(refvalmap) <- c('SNP_ID','Chr','Pos','Chip1','Chip2')
write.table(refvalmap,paste0('snpinfo_',i),quote=F,row.names=F,col.names=T)
}" > mergeSNPs.R

Rscript mergeSNPs.R

rm mergeSNPs.R

echo "#******** ----------Merging the Geno Files Together-----------  **********#"

icnt=0
for i in {1..10}
do
let "icnt=icnt+1"
echo "${icnt} geno_TestSet_${i} and geno_TrainingSet_${i}"
echo ""

awk 'NR>1' geno_TestSet_${i} > tmp_${i}
cat geno_TrainingSet_${i} tmp_${i} > finaloutfile.geno_${i}
cp snpinfo_${i} finaloutfile.snpinfo_${i}
done
rm tmpref* snpinfo* geno_*


#---------------------------------create .ctr file for FIMPUTE-----------------------#

icnt=0
for i in {1..10}
do
echo "title=*population based imputation_${i}*;
genotype_file=*finaloutfile.geno_${i}*;
snp_info_file=*finaloutfile.snpinfo_${i}*;
output_folder=*imputation${i}*;
save_genotype;
njob=5;" > imputation_${i}.ctr
sed -i 's/*/"/g' imputation_${i}.ctr
./FImpute3 imputation_${i}.ctr

rm imputation_${i}.ctr
echo ""
echo "-----------------------DONE WITH IMPUTATION_${i}------------------------"

echo "**********************************************************"
echo "******           Imputation finished             *********"
echo "******                                           *********"
echo "**********************************************************"



#--------------------------ACCURACY OF IMPUTATION --------------------#

#EXTRACTING THE IMPUTED GENO OF CHIP 2 INDIVIDUALS

#******** Extract the imputed data and make a PLINK file ***********#

cp ./imputation${i}/genotypes_imp.txt .

#I need to extract the individuals that were imputed for i.e chip2

cp ./imputation${i}/stat_anim_imp.txt .

cat stat_anim_imp.txt | awk '$2==2 {print $1,$1}' > individuals_imputed_${i}

cat genotypes_imp.txt | awk 'NR>1 {print $3}' | 
awk 'BEGIN {FS="";OFS=" "} {$1=$1; print $0}' | 
awk '{for (i=1;i<=NF;i++) { if($i==0) $i="1 1"; else if($i==1) $i="1 2"; else if($i==2) $i="2 2"; else if($i==5) $i="0 0"}  print}' > geno
cat genotypes_imp.txt | awk 'NR>1 {print $1,$1,0,0,0,-9}' > ids.txt
paste -d' ' ids.txt geno > Imputed_file.ped
cat finaloutfile.snpinfo_${i} | awk 'NR>1 {print $2,$1,0,$3}' > Imputed_file.map

rm ids.txt geno

echo ""
echo ""

./plink --file Imputed_file --make-bed --autosome-num ${chromosome} --keep individuals_imputed_${i} --recode 12 --out Imputed_Geno_${i} 

echo " "
echo " "

./plink --bfile Imputed_Geno_${i} --autosome-num ${chromosome} --recode transpose --out Imputed_Geno_${i} #to get the tped

#This shud give us the plink format of the imputed genotype for the imputed snps of the validation/test animals


#EXTRACTING THE TRUE GENO OF THE CHIP 2 INDIVIDUALS

#NOTE - here we extract the true genotype of the Test/Validation individuals and we exclude the non-imputed SNP positions 

./plink --file file --make-bed --keep individuals_imputed_${i} --recode 12 --autosome-num ${chromosome} --out TrueGeno_${i}

echo " "
echo " "

./plink --bfile TrueGeno_${i} --recode transpose --out TrueGeno_${i} --autosome-num ${chromosome}


rm individuals_imputed_${i} *.log *.nosex 
echo " "
echo ".................................................................................. "
done

#rm finaloutfile* Imputed_file* TrueGeno*.bim TrueGeno*.fam Imputed_Geno*.bim Imputed_Geno*.fam stat_anim_imp.txt minor_imputation.ctr genotypes_imp.txt Array_SNPs.txt


#--------------------------ACCURACIES------------------------

echo 'iter=0
source("calc.R")
for (i in 1:10){
iter=iter + 1
cat("Processing Set_",iter)
Result <- calc.imp.accuracy(paste0("TrueGeno_",i),paste0("Imputed_Geno_",i),paste0("TrueGeno_",i,".map"),paste0("TestSet_",i,".map"),50,format=c("tped","tped"),"simple",c(0,0))
write.table(Result[1], paste0("SNP_specific_",i),row.names = F, col.names = F, quote = F)
write.table(Result[2], paste0("Sample_specific.ALLSNP_",i), row.names = F, col.names = F, quote = F)
write.table(Result[3], paste0("Sample_specific.impSNP_",i), row.names = F, col.names = F, quote = F)
cat("Done with set_",iter)
}' > Accuracy.R

Rscript Accuracy.R

#------------------ANIMAL ACCURACY---------------------#


cat Sample_specific.impSNP_* > AnimalBased_ImputedSNPs_Accuracy_${chromosome}.txt
cat Sample_specific.ALLSNP_* > AnimalBased_AllSNPs_Accuracy_${chromosome}.txt

#rm Sample_specific.impSNP_* Sample_specific.ALLSNP_* TrueGeno* Imputed_Geno* TestSet*



#-------------------SNP ACCURACY--------------------------#

cat SNP_specific_1 | cut -f 1,2,4 -d " " > SNP_specific_new_1
for i in {2..10}
do
cat SNP_specific_${i} | cut -f 4 -d " " > SNP_specific_new_${i}
done

paste -d " " SNP_specific_new_* > SNPs_Based_Accuracy.txt  #This calculates for all snps not just imputed ones
rm SNP_specific* 


echo 'temp <- read.table("SNPs_Based_Accuracy.txt", header = F)
temp$Mean_Accuracy <- rowSums(temp[,-c(1:2)], na.rm = T)/10
write.table(temp[,c(1,2,13)], "SNPs_Based_Accuracy.txt", row.names = F, col.names = F, quote = F)' > SNPs_Accuracy.R 

Rscript SNPs_Accuracy.R

#rm *.R



#seperate array snps before applying threshold.. since they were not imputed for.. tmpval_1 consist of the Array SNPs that are on the sequence..


cat tmpval_1 | cut -f 1 -d ' ' | awk 'NR > 1'  > Exclude_ArrayDataFinal_chr_${chromosome}_SNPs.txt

#we first exclude the array SNPs bcos they were not imputed for and therefore should not be subjected to accuracy threshold

grep -vf Exclude_ArrayDataFinal_chr_${chromosome}_SNPs.txt SNPs_Based_Accuracy.txt | awk '$3 >= 0.6 {print $0}' | cat Exclude_ArrayDataFinal_chr_${chromosome}_SNPs.txt - > SNPs_for_Imputation_Chr${chromosome}.txt

#we then add the array SNPs to the imputed snps that met the specified threshold.. so we do not lose any array SNPs

wc -l SNPs_for_Imputation_Chr*.txt 


echo "#__________________________AVERAGE OF IMPUTATION ACCURACIES________________________________#"

#SUmmary of chromosome accuracy for SNPs... The mean of the chromosome... 

echo "Average SNP_Based_Accuracy chr${chromosome}"

grep -vf Exclude_ArrayDataFinal_chr_${chromosome}_SNPs.txt SNPs_Based_Accuracy.txt | awk '{ total += $3 } END { print total/NR }' > SNP_Based_Imputation_Average_chr${chromosome}

# This calculates the average of imputation accuracy of each chromosome after excluding non imputed array SNPs


#SUmmary of chromosome accuracy for SNPs... The mean of the chromosome... 

echo "Average AnimalBased_Accuracy_chr${chromosome}"

cat AnimalBased_AllSNPs_Accuracy_${chromosome}.txt | awk '{ total += $2 } END { print total/NR }' > AnimalBased_AllSNPs_Imputation_Average_chr${chromosome}

cat AnimalBased_ImputedSNPs_Accuracy_${chromosome}.txt | awk '{ total += $2 } END { print total/NR }' > AnimalBased_ImputedSNPs_Imputation_Average_chr${chromosome}



echo "#__________________________DONE WITH AVERAGE OF IMPUTATION ACCURACIES________________________________#"

echo " "
echo " "


#Some needed Results to be copied to a folder..

wc -l Exclude_ArrayDataFinal_chr_${chromosome}_SNPs.txt > Number_of_ArraySNPs_on_chr${chromosome}

grep -vf Exclude_ArrayDataFinal_chr_${chromosome}_SNPs.txt SNPs_Based_Accuracy.txt | awk '$3 >= 0.6 {print $0}' | wc -l > Number_of_ImputedSNPs_greater_than_0.6_on_chr${chromosome}

grep -vf Exclude_ArrayDataFinal_chr_${chromosome}_SNPs.txt SNPs_Based_Accuracy.txt | awk '$3 >= 0.6 {print $0}' | awk '{ total += $3 } END { print total/NR }' > SNP_Based_Imputation_Average_greater_than_0.6


mkdir ./Accuracies

cp AnimalBased_* SNP_Based_Imputation_Average_* Number_of_ArraySNPs_on_chr${chromosome} Number_of_ImputedSNPs_greater_than_0.6_on_chr${chromosome} ./Accuracies


echo "__________________________________________________________________________________________________________________________________________________________"

echo ' '
echo ' '


echo "#__________________________PROCEEDING TO REFERENCE TO TARGET POPULATION IMPUTATION________________________________#"



################################################
# Create temporary folder for analysis
FOLDER=Imputation_for_Target_Population
mkdir ${FOLDER}
cp FImpute3 plink file.* SNPs_for_Imputation_Chr${chromosome}.txt ArrayDataFinal.* Exclude_ArrayDataFinal_chr_${chromosome}_SNPs.txt ${FOLDER}/.
cd ${FOLDER}
#################################################




echo "#__________________________EXTRACTING SNPS THAT MET THE >= 0.6 THRESHOLD for TARGET POPULATION IMPUTATION________________________________#"

#This extracted SNPs would be used to carry out imputation for the TARGET POPULATION.....

./plink --bfile file \
--make-bed --double-id \
--out Reference_chr_${chromosome} \
--autosome-num ${chromosome} --extract SNPs_for_Imputation_Chr${chromosome}.txt \
--recode A

rm *.log
rm *.nosex


echo "---------------------------FILE CONVERSION DONE--------------------------------"
echo ""


echo "------------------------FILE MANIPULATION FOR IMPUTATION-------------------------"
echo ' '
echo ' '


echo "-------------------------Reference Genotype----------------------------"

awk 'NR>1 {print $2,1}' Reference_chr_${chromosome}.raw > IDs.txt
awk 'NR>1' Reference_chr_${chromosome}.raw | cut -d' ' -f7- | awk '{gsub(/NA/,5); print}' | 
awk 'BEGIN {FS=" ";OFS=""} {$1=$1; print}' | 
paste -d' ' IDs.txt - > Fimp.geno
echo 'IID Chip Call.........' > header
cat header Fimp.geno > ref.geno_chr_${chromosome}

echo "-----------------START reference SNP_INFORMATION-----------------------"

cat Reference_chr_${chromosome}.bim | awk '{print $2,$1,$4,NR}' > tmp
echo 'SNP_ID Chr Pos Chip1' > chipheader
cat chipheader tmp > ref.snpinfo_chr_${chromosome}

rm IDs.txt Fimp.geno chipheader header tmp *.raw

wc -l ref.snpinfo_chr_*
echo "------------------DONE WITH reference SNP_INFO---------------------------"
echo ""
echo "------------------------------DONE WITH FILE MANUPULATION FOR REFERENCE CHROMOSOME ${chromosome}------------------------"
echo ""
echo ""

###############################################################################################
###________________________________________________________________________________#

echo "-------------------------Target/samples Genotype----------------------------"

#splitting the array data into chromosomes

./plink --bfile ArrayDataFinal --make-bed --recode A --chr-set 29 --chr ${chromosome} --allow-extra-chr --extract Exclude_ArrayDataFinal_chr_${chromosome}_SNPs.txt --out Targetpopulation_chr_${chromosome}


echo "Targetpopulation_chr_${chromosome}.raw"

awk 'NR>1 {print $2,2}' Targetpopulation_chr_${chromosome}.raw > IDs_TestSet.txt
awk 'NR>1' Targetpopulation_chr_${chromosome}.raw | cut -d ' ' -f 7- | awk '{gsub(/NA/,5); print}' |
awk 'BEGIN {FS=" ";OFS=""} {$1=$1; print}' |
paste -d' ' IDs_TestSet.txt - > Fimp.geno
echo 'IID Chip Call.........' > header
cat header Fimp.geno > Target.geno_chr_${chromosome}
echo ""
echo ""


echo "-----------------START TARGET SNP_INFORMATION-----------------------"
echo "Targetpopulation_chr_${chromosome}.bim"
echo ""

cat Targetpopulation_chr_${chromosome}.bim | awk '{print $2,$1,$4,NR}' > tmp
echo 'SNP_ID Chr Pos Chip1' > chipheader
cat chipheader tmp > Target.snpinfo_chr_${chromosome}

echo ""
echo ""
echo "------------------DONE WITH Target population SNP_INFO---------------------------"

rm IDs*.txt header Fimp* chipheader tmp*
wc -l Target.snpinfo_* 


#-------------------------------------------------------------------------------------------------------------------###

echo "*********** Create a single MAPfile for FIMPUTE **********"

echo "Reference.snpinfo_${chromosome} and Target.snpinfo_${chromosome}"
echo ""

cp Target.snpinfo_chr_${chromosome} tmpval
cp ref.snpinfo_chr_${chromosome} tmpref



#---------------------------Rscript to merge the SNPS-----------------------------#

filtering the sequence removed some SNPs 

echo "refmap <- read.table('tmpref',header=T)
valmap <- read.table('tmpval',header=T)
refvalmap <- merge(refmap,valmap,by=1,all=T,sort=F)
refvalmap <- refvalmap[order(refvalmap[,2],refvalmap[,3]),]
refvalmap[,5] <- 0
refvalmap[which(refvalmap[,7]!=0),5] <-seq(1,length(which(refvalmap[,7]!=0)),1)
refvalmap <- refvalmap[,c(1,2,3,4,5)]
refvalmap <- na.omit(refvalmap)
colnames(refvalmap) <- c('SNP_ID','Chr','Pos','Chip1','Chip2')
write.table(refvalmap,'snpinfo',quote=F,row.names=F,col.names=T)" > mergeSNPs.R

#NOTE - Including na.omit here has removed SNPs that are on array and not on seq

module load R/4.1.2-foss-2021b
Rscript mergeSNPs.R


rm mergeSNPs.R


echo "#******** ----------Merging the Geno Files Together-----------  **********#"


echo "Target.geno_chr_${chromosome} and ref.geno_chr_${chromosome}"
echo ""

awk 'NR>1' Target.geno_chr_${chromosome} > tmp
cat ref.geno_chr_${chromosome} tmp > finaloutfile.geno_${chromosome}
cp snpinfo finaloutfile.snpinfo_${chromosome}

rm tmp* 



#---------------------------------create .ctr file for FIMPUTE-----------------------#


echo "title=*Reference-Target-Population-Imputation_${chromosome}*;
genotype_file=*finaloutfile.geno_${chromosome}*;
snp_info_file=*finaloutfile.snpinfo_${chromosome}*;
output_folder=*Final_imputation${chromosome}*;
save_genotype;
njob=5;" > imputation.ctr
sed -i 's/*/"/g' imputation.ctr

./FImpute3 imputation.ctr

rm imputation.ctr
echo ""
echo "-----------------------DONE WITH IMPUTATION-----------------------"
echo ""
echo "**********************************************************"
echo "******           Imputation finished             *********"
echo "******                                           *********"
echo "**********************************************************"

