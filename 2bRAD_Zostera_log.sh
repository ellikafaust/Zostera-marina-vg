# 2bRAD pipline log for Zostera marina
# Based of reference-based walktrhough https://github.com/z0on/2bRAD_denovo/blob/master/2bRAD_README.sh
# Last edited by Ellika Faust Nov 2024


#####################################
###     DOWNLOAD 2bRAD SCRIPTS    ###
#####################################

# downloading and installing all 2bRAD scripts in $HOME/bin (or change to whatever directory you want)
cd
mkdir bin
cd ~/bin
# cloning github repositories
git clone https://github.com/z0on/2bRAD_denovo.git
#git clone https://github.com/z0on/2bRAD_GATK.git
# move scripts to ~/bin from sub-directories
#mv 2bRAD_GATK/* .
mv 2bRAD_denovo/* .
# remove now-empty directory
rm -rf 2bRAD_denovo
rm -rf 2bRAD_GATK

# designating all .pl and .py files (perl and python scripts) as executable
chmod +x *.pl
chmod +x *.py
chmod +x *.R

# adding ~/bin to your $PATH
cd
nano .bashrc
# paste this where appropriate (note: .bashrc configuration might be specific to your cluster, consult your sysadmin if in doubt)
export PATH=$HOME/bin:$PATH

# Ctl-o, Ctl-x  (to save and exit in nano)
# log out and re-login to make sure .bashrc changes took effect

# does it work?
# try running a script from $HOME:
cd
2bRAD_trim_launch.pl
# if you get "command not found" something is wrong

#####################################
### DOWNLOAD READS FROM BASESPACE ###
#####################################

#organize fastqs
mkdir fastqs
cd fastqs

# copy across all data using rsync 

#check that you have all files
ls ./*L001_R1_001.fastq.gz | wc -l


############################
####### INDEX GENOME #######
############################

module load bioinfo-tools
module load bowtie2
module load samtools
module load picard/2.23.4

export REFERENCE_GENOME="/proj/genome/Zostera_marina.mainGenome.fasta"
export REFERENCE_GENOME_DICT="/proj/genome/Zostera_marina.mainGenome.dict"

bowtie2-build $REFERENCE_GENOME $REFERENCE_GENOME
samtools faidx $REFERENCE_GENOME
java -jar $PICARD_ROOT/picard.jar CreateSequenceDictionary R=$REFERENCE_GENOME  O=$REFERENCE_GENOME_DICT


######################################
#### CONCATENTATE LANE DUPLICATES ####
######################################

# This is a special case for VF-3360 data where there was an issue with demultiplexing and we have 6 files per samples rather than 2.
for file in ./*_barcode_Unknown_L00*_R1_001.fastq.gz
do echo `cat ${file} ${file/_barcode_Unknown_/_CGGGCT_} ${file/_barcode_Unknown_/_AGATCT_} > ${file/_barcode_Unknown_/_}`
done

# concatenate the rest
for file in ./*L001_R1_001.fastq.gz
do echo `cat ${file} ${file/_L001_/_L002_} > ${file/_L001_R1_001.fastq.gz/}_R1.fq.gz`
done

# This is a special case for VF-3360 data where there was an issue with demultiplexing and we have 6 files per samples rather than 2.
rm *_barcode_Unknown_*fq.gz
rm *_CGGGCT_*fq.gz
rm *_AGATCT_*fq.gz


#CHECK CONCATENATENATION RESULTS MAKE SENSE
ls *R1.fq.gz | wc -l




######################################
############## RUN FASTQ #############
######################################

#RUN FASTQC
module load bioinfo-tools
module load FastQC
mkdir Fastqc_Results/
> runFQC
echo '#!/bin/bash -l' > runFQC
for file in *.fq
do echo "fastqc -o Fastqc_Results -f fastq $file " >> runFQC
done

sh runFQC

#check for completion
ls Fastqc_Results/*zip | wc -l
ls Fastqc_Results/*html | wc -l



############################
###     PREP READS       ###
############################

# Trimming low-quality bases and remove pcr duplicates

# ============================
#This method is aimed to run in parallel, but still under testing on the cluster

echo '#!/bin/bash -l

#SBATCH --ntasks=10
#SBATCH -J trim
#SBATCH -o trim.%J.out
#SBATCH -e trim.%J.err
#SBATCH --time=01:00:00




for file in ./*R1.fq
do
  srun -n1 --exclusive trim2bRAD_dedup.pl input=$file > ${file}.tr0 &
done
wait' > trim.sh


# do we have expected number of *.tr0 files created?
ls -l *R1.fq.tr0 | wc -l
# NOTE: if you use the srun script the files will e created first empty so need to quality check that too


########## HOR-11 ##############
# NOTE: UJ-3099-HOR-11_S63_R1.fq is a large file and keeps running out of memory and killing the jobs
# so had to split it up
split -l 42200000 UJ-3099-HOR-11_S63_R1.fq HOR-11_

# This makes 8 files of *HOR-11_a*

# and
echo '#!/bin/bash -l

#SBATCH --ntasks=8
#SBATCH -J trimming
#SBATCH -o trimming.%J.out
#SBATCH -e trimmming.%J.err
#SBATCH --time=00:20:00


for file in ./HOR*
do
  srun -n1 --exclusive trim2bRAD_dedup.pl input=$file > ${file}.tr0 &
done
wait' > trimming.sh

rm UJ-3099-HOR-11_S63_R1.fq.tr0 # to remove empty file
####### end of HOR-11#######


# ===============================
> trimse.sh
echo '#!/bin/bash -l

#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 3:00:00
#SBATCH -J trimse
#SBATCH -o trimse.%J.out
#SBATCH -e trimse.%J.err


module load bioinfo-tools
module load cutadapt' > trimse.sh

# for reference-based analysis: trimming poor quality bases off ends:
# removing reads with qualities at ends less than Q15
# you can use -m 25 (minimum read length) if you are aligning to genome otherwise use 36

for file in *.tr0; do
echo "cutadapt -q 15,15 -m 25 -o ${file/.tr0/}.trim $file > ${file}_trimlog.txt" >> trimse.sh;
done

# execute all commands in trimse file (serial or parallel using Launcher, if your system allows)

####### HOR-11 concat ############
cat HOR-11_aa.trim HOR-11_ab.trim HOR-11_ac.trim HOR-11_ad.trim HOR-11_ae.trim HOR-11_af.trim HOR-11_ag.trim HOR-11_ah.trim > UJ-3099-HOR-11_S63_R1.fq.trim
gzip HOR-11_a*
##################################


# do we have expected number of *.trim files created?
ls -l *.trim | wc -l

#==============
# Mapping reads to reference and formatting bam files

# for reference-based:
export REFERENCE_GENOME="/proj/genome/Zostera_marina.mainGenome.fasta"
export REFERENCE_GENOME_DICT="/proj/genome/Zostera_marina.mainGenome.dict"

# mapping with --local option, enables clipping of mismatching ends (guards against deletions near ends of RAD tags)
> maps.sh
echo '#!/bin/bash -l

#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 8:00:00
#SBATCH -J maps
#SBATCH -o maps.%J.out
#SBATCH -e maps.%J.err


module load bioinfo-tools
module load bowtie2' > maps.sh


2bRAD_bowtie2_launch.pl '\.trim$' $REFERENCE_GENOME >> maps.sh

# takes a long time, so split them in separate scripts

head -n 200 maps.sh > mapshead.sh
head -n 12 maps.sh > mapstail.sh
tail -n 187 maps.sh >> mapstail.sh

# execute all commands written to maps

cat mapshead.30843205.err mapstail.30843204.err > maps
ls *.sam | wc -l  # number should match number of trim files
grep 'unpaired' maps | wc -l # just make sure this is also the same as .sam files are written part by part


# check alignment rates
>alignmentRates
for F in `ls *trim`; do
M=`grep -E '^[ATGCN]+$' $F | wc -l | grep -f - maps -A 4 | tail -1 | perl -pe 's/maps\.e\d+-|% overall alignment rate//g'` ;
echo "$F.sam $M">>alignmentRates;
done


# next stage is compressing, sorting and indexing the SAM files, so they become BAM files:
# adding read groups, validating bams
#-------------------------------

# Add read groups to sam files before converting to bam
echo '#!/bin/bash -l '> addrg.sh # 1 min for 10 files
for i in *.trim.bt2.sam; do
	newfile="$(basename $i .fq.trim.bt2.sam)"
	echo "addReadGroup.pl -i $i -o ${newfile}.fq.trim.bt2.rg.sam -r ${newfile} -s ${newfile} -l ${newfile} -u ${i:0:7} -p Illumina -c NGI" >>addrg.sh;
done

sbatch  -p core -n 1 -t 1:00:00 -o addrg.%J.out -e addrg.%J.err  addrg.sh

ls *bt2.sam | wc -l
ls *rg.sam | wc -l


# convert, sort and index from sam to bam file

export REFERENCE_GENOME="/proj/genome/Zostera_marina.mainGenome.fasta"
export REFERENCE_GENOME_DICT="/proj/genome/Zostera_marina.mainGenome.dict".dict

echo '#!/bin/bash -l
module load bioinfo-tools
module load samtools'> s2b.sh 

for file in *rg.sam; do
echo "samtools sort --reference $REFERENCE_GENOME -O bam -o ${file/.sam/}.bam $file && samtools index ${file/.sam/}.bam">>s2b.sh;
done

sbatch  -p core -n 1 -t 1:00:00 -o s2b.%J.out -e s2b.%J.err  s2b.sh


#check bams
ls *rg.bam | wc -l

# validate bams
echo '#!/bin/bash -l
module load bioinfo-tools
module load picard'>validateBams.sh # 1 min 30 sec for 10 samples
>bam_validationsummary.out
>bam_validationsummary.err
for i in *rg.bam; do
	echo "java -jar \$PICARD_ROOT/picard.jar ValidateSamFile -I $i -R $REFERENCE_GENOME -MODE SUMMARY 1>>bam_validationsummary.out 2>>bam_validationsummary.err" >>validateBams.sh;
done

sbatch  -p core -n 1 -t 1:00:00 -o validateBams.%J.out -e validateBams.%J.err validateBams.sh
cat bam_validationsummary.out

# If everything looks good, remove all sam files
# rm *bt2.sam
# rm *rg.sam


# BAM files are the input into various genotype calling / popgen programs, this is the main interim result of the analysis. Archive them.

# Next we check the average depth to see if we proceed with GATK (>10x) or ANGSD (<10x)
# if your coverage is >10x, go to GATK section below

module load bioinfo-tools
module load samtools

> average_depth
for file in *.bam; do
echo "working on $file"
D=`samtools depth $file | awk '{sum+=$3} END { print "Average = ",sum/NR}'`
echo "$file $D" >> average_depth;
done

# listing all bam filenames
ls *bam >bams
module load ANGSD

#----------- assessing base qualities and coverage depth

# angsd settings:
# -minMapQ 20 : only highly unique mappings (prob of erroneous mapping = 1%)
# -baq 1 : realign around indels (not terribly relevant for 2bRAD reads mapped with --local option)
# -maxDepth : highest total depth (sum over all samples) to assess; set to 10x number of samples
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -maxDepth 3570 -minInd 10"

# T O   D O :
TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"

angsd -b bams -r Chr01 -GL 1 $FILTERS $TODO -P 1 -out dd
# summarizing results (using modified script by Matteo Fumagalli)
Rscript ~/bin/plotQC.R prefix=dd
# proportion of sites covered at >5x:
cat quality.txt

# scp dd.pdf to laptop to look at distribution of base quality scores, fraction of sites in each sample passing coverage thresholds, and fraction of sites passing genotyping rates cutoffs. Use these to guide choices of -minQ,  -minIndDepth and -minInd filters in subsequent ANGSD runs




#############################
#        G  A  T  K         #
#############################

#        G  A  T  K
# ("hard-call" genotyping, use only for high-coverage data, >10x after deduplication)

#-------------------------------------------------
# renaming individuals
#-------------------------------------------------
# NOTE: some individual bam files that were mislabeled

# Wrong name	Correct name
# YST-19 	    STE-18-rep
# BJO-4-Rep 	ALA-05-rep
# FUR-20-Rep 	HOR-11-rep
# GRO-8-Rep 	HOG-12-rep
# HOG-12-Rep 	GRO-08-rep
# HOR-11-Rep 	FUR-20-rep
# ALA-5-Rep 	BJO-04-rep
# STE-18-Rep 	YST-19
# -----------------------------------------------


export REFERENCE_GENOME="/proj/genome/Zostera_marina.mainGenome.fasta"
export REFERENCE_GENOME_DICT="/proj/genome/Zostera_marina.mainGenome.dict".dict

ls *.bam > bams

# writing command script with SLURM header (some fields might be different on your cluster, contact your IT people!)
>unigt.sh
echo '#!/bin/bash
#SBATCH -n 20
#SBATCH -o unigt.%J.out
#SBATCH -e unigt.%J.err
#SBATCH -t 6:00:00 # took 5 hours
module load bioinfo-tools
module load GATK/3.8-0
java -Djava.io.tmpdir=$TMPDIR -jar $GATK_HOME/GenomeAnalysisTK.jar -T UnifiedGenotyper \
-R /proj/genome/Zostera_marina.mainGenome.fasta -nct 20 \
--genotype_likelihoods_model SNP \' >unigt.sh
cat bams | perl -pe 's/(\S+\.bam)/-I $1 \\/' >> unigt.sh
echo '-o primary.vcf ' >> unigt.sh


sbatch unigt.sh


# Implementing a minDP3 filter here prior to VQSR
module load bioinfo-tools
module load vcftools
vcftools --vcf primary.vcf --minDP 3 --recode --recode-INFO-all --out primary_DP3
# After filtering, kept 732 out of 732 Individuals
# Outputting VCF file...
# After filtering, kept 25761 out of a possible 25761 Sites

#----------
# Variant quality score recalibration (VQSR)

# making a tab-delimited table of clone (replicate) sample pairs
cat clonepairs.tab
# paste names of bams that are clone pairs, tab delimited, one pair per line; for example
# Ctl-O , enter, Ctl-X

# extracting "true snps" subset (reproducible across replicates)
# parameter hetPairs can vary depending on replication scheme (3 is good when you have triplicates)

# had to run with default hetPairs =1 rather than 2 as we got to few variants and GATK threw an error
# also tried som other adjustments but none did the trick...
replicatesMatch.pl vcf=primary_DP3.recode.vcf replicates=clonepairs.tab altPairs=0 hetPairs=0  > vqsr.vcf
# 25761 total SNPs
# 8951 pass hets and match filters
# 2627 show non-reference alleles
# 8951 have alterantive alleles in at least 0 replicate pair(s)
# 8951 have matching heterozygotes in at least 0 replicate pair(s)
# 2597 polymorphic
# 2627 written

# actually 8951 are written

replicatesMatch.pl vcf=primary_DP3.recode.vcf replicates=clonepairs.tab hetPairs=1  > vqsr1.vcf
# 25761 total SNPs
# 8951 pass hets and match filters
# 2627 show non-reference alleles
# 2627 have alterantive alleles in at least 1 replicate pair(s)
# 2573 have matching heterozygotes in at least 1 replicate pair(s)
# 2573 polymorphic
# 2573 written

replicatesMatch.pl vcf=primary_DP3.recode.vcf replicates=clonepairs.tab hetPairs=2  > vqsr2.vcf
# 25761 total SNPs
# 8951 pass hets and match filters
# 2627 show non-reference alleles
# 2627 have alterantive alleles in at least 1 replicate pair(s)
# 753 have matching heterozygotes in at least 2 replicate pair(s)
# 753 polymorphic
# 753 written


# determining transition-transversion ratio for true snps (will need it for tranche calibration)
module load bioinfo-tools
module load vcftools
vcftools --vcf vqsr.vcf --TsTv-summary
# Ts/Tv ratio: 2.016
vcftools --vcf vqsr1.vcf --TsTv-summary
# Ts/Tv ratio: 2.542
vcftools --vcf vqsr2.vcf --TsTv-summary
#Ts/Tv ratio: 2.576

# creating recalibration models
export REFERENCE_GENOME="/proj/genome/Zostera_marina.mainGenome.fasta"
module load bioinfo-tools
module load GATK/3.8-0

# put your actual number into the next code chunk, --target_titv

java -Djava.io.tmpdir=$HOME -jar $GATK_HOME/GenomeAnalysisTK.jar  -T VariantRecalibrator \
-R $REFERENCE_GENOME -input primary_DP3.recode.vcf -nt 1 \
-resource:repmatch,known=true,training=true,truth=true,prior=15  vqsr1.vcf \
-an QD -an DP -an FS -an MQ -an SOR -mode SNP --maxGaussians 4 \
--target_titv 2.542 -tranche 90.0 -tranche 95.0 -tranche 99.0 -tranche 100 \
-recalFile primary.recal -tranchesFile primary.recalibrate.tranches -rscriptFile primary.recalibrateSNPs.R
# -tranche 70.0 -tranche 75.0 -tranche 80.0 -tranche 85.0

# now copy all recalibrate* files to your laptop, run the R script, examine the resulting plot and tranches.pdf
# Rscript recalibrateSNPs.R

# Tested the following:
# Lowered the Gausian distributions from 6 to 4 becaue of the low number of SNPs
# Lowered prior from 30 to 10
# Testing all -an QD -an DP -an FS -an MQ -an SOR -an ReadPosRankSum -an MQRankSum and checking distributions
# displaying every 5 tranches from 55 - 100
# Testing with 660 SNPs in training (vqsr2.vcf) and 2640 (vqsr.vcf)

# The pos SNPs (green) are those which were found in the training sets passed into the VariantRecalibrator step,
# while the neg SNPs (pruple) are those which were found to be furthest away from the learned Gaussians
# and thus given the lowest probability of being true

# Do the annotation dimensions provide a clear separation between the known SNPs (most of which are true)
# and the novel SNPs (most of which are false)?

# ==================== Visualisation of annotations =============================
# This is a guide for visualising annotation in vcf files
# here I have only done it for 6 annotations FS_SOR_MQRS_RPRS_QD_MQ_DP
# for comparision see https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
module load bcftools
bcftools query primary_DP3.recode.vcf -f '%FS\t%SOR\t%MQRankSum\t%ReadPosRankSum\t%QD\t%MQ\t%DP\n' > primary_FS_SOR_MQRS_RPRS_QD_MQ_DP.txt

bgzip -c vqsr1.vcf > vqsr1.vcf.gz
tabix -p vcf vqsr1.vcf.gz
bcftools query vqsr1.vcf.gz -f '%FS\t%SOR\t%MQRankSum\t%ReadPosRankSum\t%QD\t%MQ\t%DP\n' > vqsr1_FS_SOR_MQRS_RPRS_QD_MQ_DP.txt

# To plot it run this for all files in R
'
library("ggplot2")
library("tidyr")

data <- read.delim("primary_FS_SOR_MQRS_RPRS_QD_MQ_DP.txt", header = F, col.names = c("FS", "SOR", "MQRS", "RPRS", "QD", "MQ", "DP"), na.strings = ".")
data_long <- data %>%                          # Apply pivot_longer function
  pivot_longer(colnames(data)) %>%
  as.data.frame()
ggp3 <- ggplot(data_long, aes(x = value)) +    # Draw histogram & density
  geom_histogram(aes(y = ..density..)) +
  geom_density(col = "#1b98e0", size = 1) +
  facet_wrap(~ name, scales = "free")
ggp3'


# ================================

# applying recalibration:
java -Djava.io.tmpdir=$HOME -jar $GATK_HOME/GenomeAnalysisTK.jar -T ApplyRecalibration \
-R $REFERENCE_GENOME -input primary_DP3.recode.vcf -nt 1 \
--ts_filter_level 99.0 -mode SNP \
-recalFile primary.recal -tranchesFile primary.recalibrate.tranches -o primary.recal.vcf


# Applying filters
module load bioinfo-tools
module load vcftools

vcftools --vcf primary.recal.vcf --remove-filtered-all --min-alleles 2 --max-alleles 2 --recode-INFO-all --recode --out primary.filt
# After filtering, kept 732 out of 732 Individuals
# Outputting VCF file...
# After filtering, kept 11636 out of a possible 25761 Sites

# identifying poorly genotyped individuals
vcftools --vcf primary.filt.recode.vcf --missing-indv --out primary
# look at number of sites genotyped per individual (4th column):
head primary.imiss
# see if some samples are much lower in the number of sites than others
# for example, if you want to remove samples showing less than 40000 sites:

awk '{ sum += $5 } END { if (NR > 0) print sum / NR }' primary.imiss # average 0.124936
cat primary.imiss | awk '$5>0.25' | cut -f 1 > primary.underSequenced
wc -l primary.underSequenced
# 105
grep -F -f primary.underSequenced clonepairs.tab | wc -l
# 19
grep -vF -f primary.underSequenced clonepairs.tab > primary.clonepairs_good.tab
# applying filter and selecting polymorphic biallelic loci genotyped in 90% or more individuals
# (harsh genotyping rate cutoff is strongly recommended for best quality and to avoid RAD loci affected by null
# alleles because of mutations in restriction site)


vcftools --vcf primary.filt.recode.vcf --remove primary.underSequenced --max-missing 0.9 --recode-INFO-all --recode --out primary.filt2
# After filtering, kept 628 out of 732 Individuals
# Outputting VCF file...
# After filtering, kept 9311 out of a possible 11636 Sites


# selecting only polymorphic sites (they all are in denovo pipeline!) and sites with no excess heterozygosity
grep -E "#|0/1|0/0.+1/1|1/1.+0/0" primary.filt2.recode.vcf > primary.polymorphs.vcf
hetfilter.pl vcf=primary.polymorphs.vcf maxhet=0.5 >primary.best.vcf
# 9082 total loci
# 0 dropped because fraction of missing genotypes exceeded 0.5
# 24 dropped because fraction of heterozygotes exceeded 0.5
# 9058 written

# genotypic match between pairs of replicates (the most telling one is the last one, HetsDiscoveryRate - fraction of correctly called heterozygotes; if it is under 90% perhaps use fuzzy genotyping with ANGSD - see above)
repMatchStats.pl vcf=primary.best.vcf replicates=primary.clonepairs_good.tab > primary.hetMatch.txt

#==========================================================
#==========================================================
# SOME CLONEPAIRS WHERE BAD SO HERE I AM REDOING ABOVE VQSR USING ONLY GOOD CLONEPAIRS
#==========================================================
#==========================================================


# many had low HetsDiscoveryRate, so created 3 different filters to test redo replicatesMatch.pl
cat primary.hetMatch.txt | awk -F '\t' '$10>=0.9' | wc -l # 60
cat primary.hetMatch.txt | awk -F '\t' '$10>=0.8' | wc -l # 74
cat primary.hetMatch.txt | awk -F '\t' '$10>=0.7' | wc -l # 78
cat primary.hetMatch.txt | awk -F '\t' '$10>=0.9' | cut -d ':' -f 1  > keep_clonepairs9.tab


grep -F -f keep_clonepairs9.tab clonepairs.tab > clonepairs9.tab

replicatesMatch.pl vcf=primary_DP3.recode.vcf replicates=clonepairs9.tab hetPairs=1  > vqsr9.vcf
# 25761 total SNPs
# 11657 pass hets and match filters
# 2923 show non-reference alleles
# 2923 have alterantive alleles in at least 1 replicate pair(s)
# 2847 have matching heterozygotes in at least 1 replicate pair(s)
# 2847 polymorphic
# 2847 written

module load bcftools
bgzip -c vqsr9.vcf > vqsr9.vcf.gz
tabix -p vcf vqsr9.vcf.gz
bcftools query vqsr9.vcf.gz -f '%FS\t%SOR\t%MQRankSum\t%ReadPosRankSum\t%QD\t%MQ\t%DP\n' > vqsr9_FS_SOR_MQRS_RPRS_QD_MQ_DP.txt

vcftools --vcf vqsr9.vcf --TsTv-summary
# Ts/Tv ratio: 2.362

#creating recalibration models
export REFERENCE_GENOME="/proj/genome/Zostera_marina.mainGenome.fasta"
module load bioinfo-tools
module load GATK/3.8-0

java -Djava.io.tmpdir=$HOME -jar $GATK_HOME/GenomeAnalysisTK.jar  -T VariantRecalibrator \
-R $REFERENCE_GENOME -input primary_DP3.recode.vcf -nt 1 \
-resource:repmatch,known=true,training=true,truth=true,prior=20  vqsr9.vcf \
-an QD -an DP -an FS -an MQ -an SOR -mode SNP --maxGaussians 4 \
--target_titv 2.362 -tranche 90.0 -tranche 95.0 -tranche 99.0 -tranche 100 \
-recalFile second.recal -tranchesFile recalibrate.tranches -rscriptFile recalibrateSNPs.R
# -tranche 70.0 -tranche 75.0 -tranche 80.0 -tranche 85.0
# ================================

# applying recalibration:
java -Djava.io.tmpdir=$HOME -jar $GATK_HOME/GenomeAnalysisTK.jar -T ApplyRecalibration \
-R $REFERENCE_GENOME -input primary_DP3.recode.vcf -nt 1 \
--ts_filter_level 99.0 -mode SNP \
-recalFile second.recal -tranchesFile recalibrate.tranches -o recal.vcf

# Applying filters
module load bioinfo-tools
module load vcftools

# applying filter and selecting polymorphic biallelic loci genotyped in 90% or more individuals
# (harsh genotyping rate cutoff is strongly recommended for best quality and to avoid RAD loci affected by null
# alleles because of mutations in restriction site)

vcftools --vcf recal.vcf --remove-filtered-all --min-alleles 2 --max-alleles 2 --recode-INFO-all --recode --out filt
# After filtering, kept 732 out of 732 Individuals
# Outputting VCF file...
# After filtering, kept 12444 out of a possible 25761 Sites

# identifying poorly genotyped individuals
vcftools --vcf filt.recode.vcf --missing-indv --out out2
# look at number of sites genotyped per individual (4th column):
head out2.imiss
# see if some samples are much lower in the number of sites than others
# for example, if you want to remove samples showing less than 40000 sites:

awk '{ sum += $5 } END { if (NR > 0) print sum / NR }' out2.imiss # average 0.152477
cat out2.imiss | awk '$5>0.40' | cut -f 1 > underSequenced2
wc -l underSequenced2
# 64 - 1
grep -F -f underSequenced2 clonepairs.tab | wc -l
# 15
grep -vF -f underSequenced2 clonepairs.tab > clonepairs_good2.tab

vcftools --vcf filt.recode.vcf --remove underSequenced2 --max-missing 0.9 --recode-INFO-all --recode --out filt2
# After filtering, kept 669 out of 732 Individuals
# Outputting VCF file...
# After filtering, kept 8900 out of a possible 12444 Sites

# selecting only polymorphic sites (they all are in denovo pipeline!) and sites with no excess heterozygosity
grep -E "#|0/1|0/0.+1/1|1/1.+0/0" filt2.recode.vcf >polymorphs.vcf
hetfilter.pl vcf=polymorphs.vcf maxhet=0.5 >best.vcf

# 8834 total loci
# 0 dropped because fraction of missing genotypes exceeded 0.5
# 42 dropped because fraction of heterozygotes exceeded 0.5
# 8792 written

#---------------
# Final touches

#-------------------
# KEEPING replicates
#-----------------------
thinner.pl infile=best.vcf criterion=maxAF >thinMaxaf.vcf
# 8792 total loci
# 5849 loci selected

# genotypic match between pairs of replicates (the most telling one is the last one, HetsDiscoveryRate - fraction of correctly called heterozygotes; if it is under 90% perhaps use fuzzy genotyping with ANGSD - see above)
repMatchStats.pl vcf=thinMaxaf.vcf replicates=clonepairs_good2.tab > hetMatch.txt
awk -F '\t' '$10<0.9' hetMatch.txt | wc -l
# 29 -7 = 22

#-------------------------------------------------
# renaming individuals
#-------------------------------------------------

# Rename
grep -n '#C' thinMaxaf.vcf | cut -f10-800 | sed 's/\t/\n/g' > indv
# first renaming the mislabeled replicate samples
sed -i.backup 's/YST-19/STE-18-rep/g; s/BJO-4-Rep/ALA-05-rep/g; s/FUR-20-Rep/HOR-11-rep/g; s/GRO-8-Rep/HOG-12-rep/g; s/HOG-12-Rep/GRO-08-rep/g; s/HOR-11-Rep/FUR-20-rep/g; s/ALA-5-Rep/BJO-04-rep/g; s/STE-18-Rep/YST-19/g; sx-\([1-9]_\)x-0\1xg; sx-\([1-9]-\)x-0\1xg; s/_S.*$//g; s/a$//g; s/_R1//g; s/KO\-DON/KOD/g; s/repl/rep/g; s/Rep/rep/g' indv
# create file for renaming
paste indv.backup indv > rename_indv
#rename using bcftools

module load bcftools
bcftools reheader -s rename_indv thinMaxaf.vcf > thinMaxaf_renamed.vcf

vcftools --vcf thinMaxaf_renamed.vcf --missing-indv --out minDP

awk '{ sum += $5 } END { if (NR > 0) print sum / NR }' minDP.imiss # average 0.022


#----------
# 25% MISSINGNES
# -----------
cat minDP.imiss | awk '$5>0.25' | cut -f 1  > ind_miss_25
wc -l ind_miss_25
# 2 -1 = 1

vcftools --vcf thinMaxaf_renamed.vcf --remove ind_miss_25 --max-missing 0.9 --recode-INFO-all --recode --out all_miss25
# After filtering, kept 668 out of 669 Individuals
# Outputting VCF file...
# After filtering, kept 5849 out of a possible 5849 Sites


# het match between replicates
# renaming clone pairs files 
sed 's/YST-19/STE-18-rep/g; s/BJO-4-Rep/ALA-05-rep/g; s/FUR-20-Rep/HOR-11-rep/g; s/GRO-8-Rep/HOG-12-rep/g; s/HOG-12-Rep/GRO-08-rep/g; s/HOR-11-Rep/FUR-20-rep/g; s/ALA-5-Rep/BJO-04-rep/g; s/STE-18-Rep/YST-19/g; sx-\([1-9]_\)x-0\1xg; sx-\([1-9]-\)x-0\1xg; s/_S.*\t/\t/g; s/a\t/\t/g; s/_S.*$//g; s/a$//g; s/_R1//g; s/KO\-DON/KOD/g; s/repl/rep/g; s/Rep/rep/g' clonepairs_good2.tab > clonepairs_renamed.tab

grep -F -f indv clonepairs_renamed.tab | grep -vF -f ind_miss_25 > clonepairs_renamed_25.tab
repMatchStats.pl vcf=all_miss25.recode.vcf replicates=clonepairs_renamed_25.tab > hetMatch_all_miss25.txt
awk -F '\t' '$10<0.9' hetMatch.txt | wc -l
# 29 -7 = 21
awk -F '\t' '$10<0.9' hetMatch_all_miss25.txt | wc -l
# 29 -7 = 21


# here only extracting individuals used for this study
grep -E "VH|KOD|STE" indv | sort > ellika_indv_all

# ind list without replicates
grep 'rep' ellika_indv_all > ellika_indv_uniq

vcftools --vcf all_miss25.recode.vcf --keep ellika_indv_all --recode --recode-INFO-all --out zostera_300123_ellika_miss25 # ellikas
# After filtering, kept 404 out of 668 Individuals
# Outputting VCF file...
# After filtering, kept 5849 out of a possible 5849 Sites


#SORT
# sort individuals according to a list
# pop order
echo "KOD
GAS
KAV
JOR
HAV
LIN
FIS
VHA
VFL
SVA
VIK
STE
KAK
WAL
NIN
NOV
SBR
HAL" > pop_order

# reorders the file
:> ellika_indv_all_sorted
while read TESTLINE; do
  grep "${TESTLINE}" ellika_indv_all >> ellika_indv_all_sorted
done < pop_order

# a pop list of individuals pop code
sed 's/\-/\t/g' ellika_indv_all_sorted | cut -f 3 > ellika_pop_all_sorted


bcftools view -S ellika_indv_all_sorted --force-samples zostera_300123_ellika_miss25.recode.vcf -O v -o zostera_300123_ellika_miss25_sorted.vcf

export VCF_FILE="zostera_300123_ellika_miss25_sorted.vcf"
vcftools --vcf $VCF_FILE --freq --out ${VCF_FILE/.vcf/}
vcftools --vcf $VCF_FILE --counts --out ${VCF_FILE/.vcf/}
vcftools --vcf $VCF_FILE --depth --out ${VCF_FILE/.vcf/}
vcftools --vcf $VCF_FILE --site-depth --out ${VCF_FILE/.vcf/}
vcftools --vcf $VCF_FILE --site-mean-depth --out ${VCF_FILE/.vcf/}
vcftools --vcf $VCF_FILE --missing-indv --out ${VCF_FILE/.vcf/}
vcftools --vcf $VCF_FILE --missing-site --out ${VCF_FILE/.vcf/}
vcftools --vcf $VCF_FILE --singletons --out ${VCF_FILE/.vcf/}
vcftools --vcf $VCF_FILE --plink --out ${VCF_FILE/.vcf/}
vcftools --vcf $VCF_FILE --hardy --out ${VCF_FILE/.vcf/}
vcftools --vcf $VCF_FILE --counts --out ${VCF_FILE/.vcf/}
#vcftools --vcf $VCF_FILE --geno-r2 --out ${VCF_FILE/.vcf/}
vcftools --vcf $VCF_FILE --relatedness --out ${VCF_FILE/.vcf/}
vcftools --vcf $VCF_FILE --site-pi --out ${VCF_FILE/.vcf/}
vcftools --vcf $VCF_FILE --het --out ${VCF_FILE/.vcf/}


#------------------------
# FINAL
#-----------------------
# 5849 Sites

# get list of individuals and of positions (removing any scaffold positions)
module load bcftools 
bcftools query -l zostera_300123_ellika_miss25_sorted.vcf > inds
cut -c9- inds > ind_rename

# rename samples names to not includ sequencing run 
bcftools reheader -s ind_rename zostera_300123_ellika_miss25_sorted.vcf  -o zostera_5849.vcf 

# remove monomorphic sites
bcftools view --exclude 'AC==0 || AC==AN' zostera_5849.vcf  -o zostera_4864.vcf 




# ----------------------------
# ----Nucleotide diveristy----
# ----------------------------
 

# ---------- pi.sh ----------
module load vcftools
VCF_FILE="zostera_filt_25_named.vcf"


# for meadows 

for i in KOD \
GAS \
KAV \
JOR \
HAV \
LIN \
FIS \
VHA \
VFL \
SVA \
VIK \
STE \
KAK \
WAL \
NIN \
NOV \
SBR \
HAL; do
    awk -F ',' -v grp="$i" '$3 == grp {print $8}' clusters.csv > keep$i
    vcftools --vcf $VCF_FILE --site-pi --keep keep$i --out $i
done

# for ref and impacted
for i in Reference Impacted ; do
    awk -F ',' -v grp="$i" '$4 == grp {print $8}' clusters.csv > keep$i
    vcftools --vcf $VCF_FILE --site-pi --keep keep$i --out $i
done

# for cluster 
for i in C1 C2 C3 C4 C5 ; do
    awk -F ',' -v grp="$i" '$5 == grp {print $8}' clusters.csv > keep$i
    vcftools --vcf $VCF_FILE --site-pi --keep keep$i --out $i
done

# for north south
for i in North South ; do
    awk -F ',' -v grp="$i" '$6 == grp {print $8}' clusters.csv > keep$i
    vcftools --vcf $VCF_FILE --site-pi --keep keep$i --out $i
done

# for size
for i in Small Large; do
    awk -F ',' -v grp="$i" '$7 == grp {print $8}' clusters.csv > keep$i
    vcftools --vcf $VCF_FILE --site-pi --keep keep$i --out $i
done



# concatenate the files
head -n1 KOD.sites.pi | awk -v grp="GROUP" '{print grp "\t" $0}' > group_pi.txt
for i in KOD \
GAS \
KAV \
JOR \
HAV \
LIN \
FIS \
VHA \
VFL \
SVA \
VIK \
STE \
KAK \
WAL \
NIN \
NOV \
SBR \
HAL \
C1 C2 C3 C4 C5 North South Reference Impacted Small Large; do
	tail -n +2 $i.sites.pi | awk -v grp="$i" '{print grp "\t" $0}' >> group_pi.txt
done


# mean pi 
echo -e 'site\tmean' >sum_pi.txt
for i in KOD \
GAS \
KAV \
JOR \
HAV \
LIN \
FIS \
VHA \
VFL \
SVA \
VIK \
STE \
KAK \
WAL \
NIN \
NOV \
SBR \
HAL \
C1 C2 C3 C4 C5 North South Reference Impacted Small Large; do
	awk -F '\t' -v grp="$i"  'NR > 1 { sum += $3; count++ } END { print grp "\t" sum / count }' $i.sites.pi  >>sum_pi.txt
done







