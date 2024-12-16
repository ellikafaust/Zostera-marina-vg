


# get list of individuals and of positions (removing any scaffold positions)
module load bcftools 
bcftools query -l zostera_filt_25.recode.vcf > inds
cut -c9- inds > ind_rename

bcftools query -f '%CHROM\t%POS\n' zostera_filt_25.recode.vcf | grep "Chr" >pos

# rename samples names to not includ sequencing run 
bcftools reheader -s ind_rename zostera_filt_25.recode.vcf  -o zostera_filt_25_named.vcf 

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
    vcftools --vcf $VCF_FILE --keep keep$i --positions pos --recode --out group$i
    ./currentNe/currentNe -t 8 group${i}.recode.vcf 6
    rm group${i}.recode.vcf
    rm group${i}.log
done

# for ref and impacted
for i in Reference Impacted ; do
    awk -F ',' -v grp="$i" '$4 == grp {print $8}' clusters.csv > keep$i
    vcftools --vcf $VCF_FILE --keep keep$i --positions pos --recode --out group$i
    ./currentNe/currentNe -t 8 group${i}.recode.vcf 6
    rm group${i}.recode.vcf
    rm group${i}.log
done

# for cluster 
for i in C1 C2 C3 C4 C5 ; do
    awk -F ',' -v grp="$i" '$5 == grp {print $8}' clusters.csv > keep$i
    vcftools --vcf $VCF_FILE --keep keep$i --positions pos --recode --out group$i
    ./currentNe/currentNe -t 8 group${i}.recode.vcf 6
    rm group${i}.recode.vcf
    rm group${i}.log
done

# for north south
for i in North South ; do
    awk -F ',' -v grp="$i" '$6 == grp {print $8}' clusters.csv > keep$i
    vcftools --vcf $VCF_FILE --keep keep$i --positions pos --recode --out group$i
    ./currentNe/currentNe -t 8 group${i}.recode.vcf 6
    rm group${i}.recode.vcf
    rm group${i}.log
done

# for size
for i in Small Large; do
    awk -F ',' -v grp="$i" '$7 == grp {print $8}' clusters.csv > keep$i
    vcftools --vcf $VCF_FILE --keep keep$i --positions pos --recode --out group$i
    ./currentNe/currentNe -t 8 group${i}.recode.vcf 6
    rm group${i}.recode.vcf
    rm group${i}.log
done



####################
# currentNe
####################

# Define the header
echo "SNPs miss sibs F He Ho Ne CI_dw CI_up Ne_corrected CI_dw_corrected CI_up_corrected" > ne_sum 

# Extract the desired lines from each file and append them to the file
for file in group*.recode_currentNe_OUTPUT.txt; do
	  echo $file
    SNPs=$(awk 'NR==17' "$file")
    miss=$(awk 'NR==25' "$file")
    sibs=$(awk 'NR==29' "$file")
    F=$(awk 'NR==33' "$file")
    Ho=$(awk 'NR==43' "$file")
    He=$(awk 'NR==45' "$file")
    Ne=$(awk 'NR==50' "$file")
    CI_dw=$(awk 'NR==56' "$file")
    CI_up=$(awk 'NR==58' "$file")
    Ne_corrected=$(awk 'NR==66' "$file")
    CI_dw_corrected=$(awk 'NR==72' "$file")
    CI_up_corrected=$(awk 'NR==74' "$file")
    echo "$SNPs $miss $sibs $F $He $Ho $Ne $CI_dw $CI_up $Ne_corrected $CI_dw_corrected $CI_up_corrected" >> ne_sum
done


####### GONE #####



module load bcftools
module load vcftools
VCF_FILE="zostera_filt_25_named.vcf"


# To get the map file, we need to create a chromosome map in which we assume number for chromosomes.
bcftools view -H $VCF_FILE | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom-map45
# replacing contig names with numbers, as required by GONE:
nl chrom-map45 > cc
cut -f 1 cc > c1
cut -f 2 cc > c2
paste c2 c1 > chrom-map45
head -n 6 chrom-map45 >chrom-map6

# make list of individuals to extract
# use list to create ped file
# extract the ped data from column 7 and forwards
# make GONE appropriate first 6 columns
# paste the two files together and overwrite the ped file
# remove temporary files

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
    vcftools --vcf $VCF_FILE --keep keep$i --positions pos --plink --chrom-map chrom-map6 --out $i
    cut -f7- ${i}.ped > temp${i}
    awk '{ printf "1\tIND%s\t0\t0\t1\t-9\n", $1 }' keep${i} > ind${i}
    paste ind${i} temp${i} > ${i}.ped
    rm temp${i}*
    rm ind${i}
done


#change input parameter so that cm/mb is 2.736054043, 
# this is equal to one recombination per generation for the average genome siez 36.54898567 and 100cm per chromosome



# and finally run GONE

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
 bash script_GONE.sh $i
done

#WARNING: GONE overwrites temporary files and the output file "outfileHWD", and therefore 
# (1) multiple instances of GONE cannot be run in the same directory; 
# (2) running GONE using a loop will only overwrite the output file "outfileHWD" and time and seed file but HW can be found in the output and time is in the slurm.

# average number of SNPs used per chromosome
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
grep "NSNP_calculations=" OUTPUT_$i | cut -d "=" -f 2 
done | awk 'BEGIN{s=0;}{s+=$1;}END{print s/NR;}'  
# 199.448

# Ne estimates
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
    ne_value=$(awk 'NR==3' Output_Ne_$i | cut -f 2)  # Extract the value from the 2nd column of line 3
    echo "$i $ne_value"  # Print the population identifier and its corresponding ne value
done > recent_ne.txt


# last 100 generations Ne
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
awk 'NR>=3 && NR<=502' Output_Ne_$i | awk -v grp="$i" '{print grp "\t" $0}' >> group_ne.txt
done