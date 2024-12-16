# log for calculation of runs of homozygosity



# ----------------------------
# ----------- ROH ------------
# ----------------------------

## Genome Coverage ##

# First we want to figure out which parameters are best suited for our data.
# For this we test parameters on simulated individuals with 100% homzygosity.
# The idea is with 100% homozygot individuals you should have a close to 100% coverage of the genome

module load plink
module load vcftools
module load bcftools




VCF_FILE="zostera_300123_ellika_miss25_sorted.vcf"
#renaming the ind files so I can remove replicates
bcftools query -l zostera_300123_ellika_miss25_sorted.vcf > inds
paste inds <(cut -c 9- inds) > inds_rename

cut -c 9- inds | grep "rep" > reps
cut -c 9- inds | grep -v "rep" > no_reps
grep -vf no_reps reps >> no_reps

# rename samples names to not includ sequencing run and rep
bcftools reheader --sample inds_rename $VCF_FILE -o zostera_300123_ellika_miss25_sorted_named.vcf

bcftools view -S no_reps zostera_300123_ellika_miss25_sorted_named.vcf -o zoster_341_named.vcf

VCF_FILE="zoster_341_named.vcf"

# Testing the individual with the most missingnes first


# make all individuals 100% homozygous 
bcftools view -h $VCF_FILE > ${VCF_FILE}.temp
#bcftools view -H $VCF_FILE | sed s'/\.\/\./0\/0/ ; s/0\/1/0\/0/'g >> ${VCF_FILE}.temp
bcftools view -H $VCF_FILE | sed s'/0\/1/0\/0/'g >> ${VCF_FILE}.temp
mv ${VCF_FILE}.temp ${VCF_FILE/.vcf/_hom.vcf}




VCF_FILE="zoster_341_named_hom.vcf"

# this is a for loop which calculates roh for a range of max missing values and then summaries the individual output in a file called rod_stats.txt

bcftools query -l ${VCF_FILE} > inds

echo "miss ID nROH length nSNP fracSNP averagekb coverage n500 l500 " >roh_stats.txt
sites="5849"
genomesize="260492"

for miss in 1 2 3 4 5
do
    plink --vcf $VCF_FILE --allow-extra-chr --homozyg --homozyg-window-missing $miss --homozyg-snp 10 --homozyg-window-snp 7 --homozyg-kb 300 --homozyg-gap 400 --out roh_miss$miss
    for ID in $(cut -f1 inds)
    do
        echo "working with $ID and missingnes $miss"
        nROH=$(grep $ID roh_miss$miss.hom  |sed -e 's/^[ \t]*//;s/ \+/\t/g'| awk -F"\t" '{ sum+=$9>1} END {OFMT="%.2f";print sum}')
        length=$(grep $ID roh_miss$miss.hom  |sed -e 's/^[ \t]*//;s/ \+/\t/g'| awk -F"\t" '{ sum+=$9} END {OFMT="%.2f";print sum}')
        nSNP=$(grep $ID roh_miss$miss.hom  |sed -e 's/^[ \t]*//;s/ \+/\t/g'| awk -F"\t" '{ sum+=$10} END {OFMT="%.2f";print sum}')
        averagekb=$(bc -l <<< "scale=2;$length/$nROH")
        coverage=$(bc -l <<< "scale=4;$length/$genomesize")
        fracSNP=$(bc -l <<< "scale=4;$nSNP/$sites")
        n500=$(grep $ID roh_miss$miss.hom  |sed -e 's/^[ \t]*//;s/ \+/\t/g'| awk -F"\t" '{ sum+=$9>500} END {print sum}')
        l500=$(grep $ID roh_miss$miss.hom  |sed -e 's/^[ \t]*//;s/ \+/\t/g'| awk -F"\t" '{ if ($9 > 500) sum += $9 } END {OFMT="%.2f";print sum}')
        #n1000=$(grep $ID roh_miss$miss.hom  |sed -e 's/^[ \t]*//;s/ \+/\t/g'| awk -F"\t" '{ sum+=$9>1000} END {print sum}')
        #l1000=$(grep $ID roh_miss$miss.hom  |sed -e 's/^[ \t]*//;s/ \+/\t/g'| awk -F"\t" '{ if ($9 > 1000) sum += $9 } END {print sum}')
        #n5000=$(grep $ID roh_miss$miss.hom  |sed -e 's/^[ \t]*//;s/ \+/\t/g'| awk -F"\t" '{ sum+=$9>5000} END {print sum}')
        #l5000=$(grep $ID roh_miss$miss.hom  |sed -e 's/^[ \t]*//;s/ \+/\t/g'| awk -F"\t" '{ if ($9 > 5000) sum += $9 } END {print sum}')
        echo "$miss $ID $nROH $length $nSNP $fracSNP $averagekb $coverage $n500 $l500" >>roh_stats.txt
    done
done


#---------------------------

# Now the testing different levels of gaps


VCF_FILE="zoster_341_named.vcf"

for gap in 200 300 400 500 1000
do
    plink --vcf $VCF_FILE --allow-extra-chr --homozyg --homozyg-window-missing 2 --homozyg-snp 10 --homozyg-window-snp 7 --homozyg-kb 300 --homozyg-gap $gap --out roh_gap$gap
done 
# Calculate total length and number of roh
for gap in 200 300 400 500 1000
do
    echo "gap$gap"
    nROH=$(tail -n+2 roh_gap$gap.hom  |sed -e 's/^[ \t]*//;s/ \+/\t/g'| awk -F"\t" '{ sum+=$9>1} END {OFMT="%.2f";print sum}')
    length=$(tail -n+2 roh_gap$gap.hom  |sed -e 's/^[ \t]*//;s/ \+/\t/g'| awk -F"\t" '{ sum+=$9} END {OFMT="%.2f";print sum}')
    echo "$nROH $length"
done



# testing het levels
for het in $(seq 10)
do
    plink --vcf $VCF_FILE --allow-extra-chr --homozyg --homozyg-snp 10 --homozyg-het $het --homozyg-window-snp 7 --homozyg-window-missing 2 --homozyg-kb 300 --homozyg-gap 400 --out roh_het$het
done





# Calculate total length and number of roh
for het in $(seq 10 )
do
    echo "het$het"
    nROH=$(tail -n+2 roh_het$het.hom  |sed -e 's/^[ \t]*//;s/ \+/\t/g'| awk -F"\t" '{ sum+=$9>1} END {OFMT="%.2f";print sum}')
    length=$(tail -n+2 roh_het$het.hom  |sed -e 's/^[ \t]*//;s/ \+/\t/g'| awk -F"\t" '{ sum+=$9} END {OFMT="%.2f";print sum}')
    echo "$nROH $length"
done




echo "het ID nROH length nSNP fracSNP averagekb FROH n500 l500 FROH500 n1000 l1000 FROH1000 n5000 l5000 FROH1000" >roh_het_g400_stats.txt
sites="5849"
L_300="240209"
L_500="237092"
L_1000="229626"
L_5000="228409"

for het in $(seq 10 )
do
    for ID in $(cut -f1 inds)
    do
        echo "working with $ID and het $het"
        nROH=$(grep $ID roh_het$het.hom  |sed -e 's/^[ \t]*//;s/ \+/\t/g'| awk -F"\t" '{ sum+=$9>1} END {OFMT="%.2f";print sum}')
        length=$(grep $ID roh_het$het.hom  |sed -e 's/^[ \t]*//;s/ \+/\t/g'| awk -F"\t" '{ sum+=$9} END {OFMT="%.2f";print sum}')
        nSNP=$(grep $ID roh_het$het.hom  |sed -e 's/^[ \t]*//;s/ \+/\t/g'| awk -F"\t" '{ sum+=$10} END {OFMT="%.2f";print sum}')
        averagekb=$(bc -l <<< "scale=2;$length/$nROH")
        FROH=$(bc -l <<< "scale=4;$length/$L_300")
        fracSNP=$(bc -l <<< "scale=4;$nSNP/$sites")
        n500=$(grep $ID roh_het$het.hom  |sed -e 's/^[ \t]*//;s/ \+/\t/g'| awk -F"\t" '{ sum+=$9>500} END {print sum}')
        l500=$(grep $ID roh_het$het.hom  |sed -e 's/^[ \t]*//;s/ \+/\t/g'| awk -F"\t" '{ if ($9 > 500) sum += $9 } END {OFMT="%.2f";print sum}')
        FROH500=$(bc -l <<< "scale=4;$l500/$L_500")
        n1000=$(grep $ID roh_het$het.hom  |sed -e 's/^[ \t]*//;s/ \+/\t/g'| awk -F"\t" '{ sum+=$9>1000} END {print sum}')
        l1000=$(grep $ID roh_het$het.hom  |sed -e 's/^[ \t]*//;s/ \+/\t/g'| awk -F"\t" '{ if ($9 > 1000) sum += $9 } END {print sum}')
        FROH1000=$(bc -l <<< "scale=4;$l1000/$L_1000")
        n5000=$(grep $ID roh_het$het.hom  |sed -e 's/^[ \t]*//;s/ \+/\t/g'| awk -F"\t" '{ sum+=$9>5000} END {print sum}')
        l5000=$(grep $ID roh_het$het.hom  |sed -e 's/^[ \t]*//;s/ \+/\t/g'| awk -F"\t" '{ if ($9 > 5000) sum += $9 } END {print sum}')
        FROH5000=$(bc -l <<< "scale=4;$l1000/$L_1000")
        echo "$het $ID $nROH $length $nSNP $fracSNP $averagekb $FROH $n500 $l500 $FROH500 $n1000 $l1000 $FROH1000 $n5000 $l5000 $FROH1000" >>roh_het_g400_stats.txt
    done
done


# final run
VCF_FILE="zoster_341_named.vcf"

plink --vcf $VCF_FILE --allow-extra-chr --homozyg group --homozyg-snp 10 --homozyg-het 2 --homozyg-window-missing 2 --homozyg-window-snp 7 --homozyg-kb 300 --homozyg-gap 400 --out roh_final

