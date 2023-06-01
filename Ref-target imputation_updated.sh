#!/bin/bash

chromosome=$1
rm -r ~/thesis/Accuracy/Imputation_Accuracy_chr_${chromosome}/Imputation_for_Target_Population
cd ~/thesis/Accuracy/Imputation_Accuracy_chr_${chromosome}/

################################################
# Create temporary folder for analysis
FOLDER=Imputation_for_Target_Population
mkdir ${FOLDER}
cp FImpute3 plink SNPs_for_Imputation_Chr${chromosome}.txt Array_SNPs.txt ${FOLDER}/.
cd ${FOLDER}
#################################################

#__________________________________________________________-_-__-__-__-__-__-__-__-__-__-__-__-__-__-_____________________________________________________-
#Reference_allle ACTG -- #The new dataset is all merged.. I splitted to chromosomes
./plink --file ~/thesis/Seq_MLA/SeqNew/SeqCorrect --make-bed \
--chr ${chromosome} --allow-extra-chr \
--chr-set 29 --double-id \
--out tmp_ref_FiltSNPs_${chromosome}

cat tmp_ref_FiltSNPs_${chromosome}.bim | awk '{print $2,$5}' > reference_allele_seq_chr${chromosome}.txt

#Array ACTG
grep -wf Array_SNPs.txt reference_allele_seq_chr${chromosome}.txt > reference_allele_array_chr${chromosome}.txt

#__________________________________________________________-_-__-__-__-__-__-__-__-__-__-__-__-__-__-_____________________________________________________-

echo "-------------------------CONVERTING TO BINARY & RECODING------------------------------"

#To update the name in the fam file and recode
cp ~/thesis/updatedFAM_SampleID.fam tmp_ref_FiltSNPs_${chromosome}.fam

./plink --bfile tmp_ref_FiltSNPs_${chromosome} --aec \
--recode 12 --chr-set 29 --double-id \
--reference-allele reference_allele_seq_chr${chromosome}.txt \
--out ref_FiltSNPs_${chromosome} #The reference-allele flag keeps the ACTG format before converting to 12 format....


./plink --file ref_FiltSNPs_${chromosome} --aec \
--chr-set 29 --double-id \
--make-bed --recode A \
--out ref_FiltSNPs_${chromosome}
#rm *.log *.map *.ped


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


#__________________________________________________________-_-__-__-__-__-__-__-__-__-__-__-__-__-__-_____________________________________________________-

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


echo "------------------------FILE MANIPULATION FOR REF-TARGET IMPUTATION-------------------------"
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

# we are only using 2017 Target Population Individuals... about 3000

cat ~/thesis/Array_MLA/ArrayNew/ArrayCorrect.fam | grep ^2017 | cut -f 1-2 -d ' ' > TargetIndividuals_2017_GWAS.txt


#splitting the array data into chromosomes and 
## change array format into plink 12

./plink --file ~/thesis/Array_MLA/ArrayNew/ArrayCorrect --recode 12 \
--chr-set 29 --nonfounders \
--reference-allele reference_allele_array_chr${chromosome}.txt \
--chr ${chromosome} --allow-extra-chr \
--out Targetpopulation_chr_${chromosome} --keep TargetIndividuals_2017_GWAS.txt

./plink --file Targetpopulation_chr_${chromosome} --make-bed \
--extract Array_SNPs.txt \
--chr-set 29 --recode A --out Targetpopulation_chr_${chromosome}


echo "Targetpopulation_chr_${chromosome}.raw" #This is for the 3185 2017 individuals

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

rm IDs*.txt header Fimp* chipheader tmp
wc -l Target.snpinfo_* 


#-------------------------------------------------------------------------------------------------------------------###

echo "*********** Create a single MAPfile for FIMPUTE **********"

echo "Reference.snpinfo_${chromosome} and Target.snpinfo_${chromosome}"
echo ""

cp Target.snpinfo_chr_${chromosome} tmpval
cp ref.snpinfo_chr_${chromosome} tmpref



#---------------------------Rscript to merge the SNPS-----------------------------#

#filtering the sequence removed some SNPs 

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


echo "title=*Reference-Target Population Imputation_${chromosome}*;
genotype_file=*finaloutfile.geno_${chromosome}*;
snp_info_file=*finaloutfile.snpinfo_${chromosome}*;
output_folder=*Final_imputation${chromosome}*;
save_genotype;
njob=5;" > imputation.ctr
sed -i 's/*/"/g' imputation.ctr

./FImpute3 imputation.ctr

rm imputation.ctr
echo ""
echo ""
echo "**********************************************************"
echo "******           Imputation finished             *********"
echo "******                                           *********"
echo "**********************************************************"

echo "-----------------------DONE WITH IMPUTATION-----------------------"


echo "#------------------------EXTRACT IMPUTED SNPS FOR THE INDIVIDUALS-------------------#"

#EXTRACTING THE IMPUTED GENO OF CHIP 2 INDIVIDUALS

#******** Extract the imputed data and make a PLINK file ***********#

cp ./Final_imputation${chromosome}/genotypes_imp.txt .

#I need to extract the individuals that were imputed for i.e chip2 Target Population #3000+

cp ./Final_imputation${chromosome}/stat_anim_imp.txt .

cat stat_anim_imp.txt | awk '$2==2 {print $1,$1}' > individuals_imputed_${chromosome}

cat genotypes_imp.txt | awk 'NR>1 {print $3}' | 
awk 'BEGIN {FS="";OFS=" "} {$1=$1; print $0}' | 
awk '{for (i=1;i<=NF;i++) { if($i==0) $i="1 1"; else if($i==1) $i="1 2"; else if($i==2) $i="2 2"; else if($i==5) $i="0 0"}  print}' > geno
cat genotypes_imp.txt | awk 'NR>1 {print $1,$1,0,0,0,-9}' > ids.txt
paste -d' ' ids.txt geno > GWAS_chr${chromosome}.ped


cat finaloutfile.snpinfo_${chromosome} | awk 'NR>1 {print $2,$1,0,$3}' > GWAS_chr${chromosome}.map

rm ids.txt geno

echo ""
echo ""

./plink --file GWAS_chr${chromosome} --make-bed --autosome-num ${chromosome} \
--keep individuals_imputed_${chromosome} --out GWAS_Data_chr${chromosome} 

#This creates the file needed for GWAS.
