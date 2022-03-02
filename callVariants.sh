###############################################################################
# Nugen pipeline for calling variants from targeted sequencing assay
# Author: Stefan Gavriliuc May - June 2021 (Poissant lab)
# 
# Calls variants using the suggested Nugen pipeline
# Keep track of number of reads retained after each step
# Convert .vcfs to table using GATK VariantsToTable utility
###############################################################################

conda activate variant_analysis

# Specify reference genome and reference genome build for bowtie2
# If you do not have the EquCab2 reference genome, it can be obtained with the following
# The file is large (7.7 gigabytes compressed) 
#wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Equus_caballus/Ensembl/EquCab2/Equus_caballus_Ensembl_EquCab2.tar.gz
#tar -xzvf Equus_caballus_Ensembl_EquCab2.tar.gz
BOWTIE2_REFERENCE="/home/stefan/Documents/Research/msc/EquCab2/Sequence/Bowtie2Index/genome"
FASTA_REFERENCE="/home/stefan/Documents/Research/msc/EquCab2/Sequence/WholeGenomeFasta/genome.fa"

###############################################################################
# 0. Preliminary analysis: create vcf from plink files, get raw read counts
###############################################################################
plink --file data/SNPchip_genotypes_in_Plink_format/Sable_Nugen_Illumina_Affy_combined \
    --horse --recode vcf # makes plink.vcf
# Change chromosome 32 to X, then sort VCF by chromosome & position
sed 's/^32/X/g' plink.vcf | ./extra_scripts/vcfsort.sh - > ref.vcf
rm plink.* # remove unnecessary files from plink

# Create sample statistics file, print out raw read counts
rm sampleStatistics.csv
for sample in seqs/*fastq.gz
do
    base=`basename -s "_L001_R1_001.fastq.gz" $sample`
    zcat $sample | echo $base $((`wc -l`/4)) >> sampleStatistics.csv
done

###############################################################################
# 1. Use trimgalore with relaxed values 
# (length = 15, stringency = 7)
###############################################################################
for sample in seqs/*R1_001.fastq.gz; do
    echo $sample
    base=`basename -s "_L001_R1_001.fastq.gz" $sample`
    echo $base
    trim_galore --length 15 --stringency 7 --clip_R1 40 \
        $sample -o intermediates
done

# read counts using relaxed parameters
rm tmp.txt
for sample in intermediates/*_trimmed.fq.gz; do
    base=`basename -s "_L001_R1_001_trimmed.fastq.gz" $sample`
    #zcat $sample | echo $base $((`wc -l`/4)) | cut -f2 -d' '
    zcat $sample | echo $base $((`wc -l`/4)) | cut -f2 -d' ' >> tmp.trim.txt
done

# Append to sample statistics file, create a temporary file for renaming
paste -d " " sampleStatistics.csv tmp.trim.txt > tmp.txt
mv tmp.txt sampleStatistics.csv
rm tmp.trim.txt

###############################################################################
# 2. Single end read alignment
# Relaxed trimgalore, local alignment
###############################################################################
for sample in intermediates/*trimmed.fq.gz ; do
    echo $sample
    base=`basename -s "_trimmed.fq.gz" $sample`
    echo $base
    bowtie2 -p 2 --very-sensitive-local -U $sample \
        -x $BOWTIE2_REFERENCE | samtools sort -o $base".bam"
done

mv *.bam intermediates

rm tmp.map.txt
# Number of reads that mapped to reference genome (EquCab2)
for sample in intermediates/*.bam; do
    base=`basename -s ".bam" $sample`
    echo "`samtools view -c -F 260 $sample` " >> tmp.map.txt
done

paste -d " " sampleStatistics.csv tmp.map.txt > tmp.txt
mv tmp.txt sampleStatistics.csv
rm tmp.map.txt

# Add read group information to bam files to preserve sample names
for sample in intermediates/*.bam ; do
    echo $sample
    base=`basename -s ".bam" $sample`
    echo $base

    gatk AddOrReplaceReadGroups \
        I=$sample \
        O=$base"_rg.bam" \
        RGID=1 \
        RGLB=lib1 \
        RGPL=nugen \
        RGPU=unit1 \
        RGSM=$sample

    # Index alignment files
    samtools index $base"_rg.bam"
done

mv *_rg.bam* intermediates

###############################################################################
# 3. call variants using mpileup
###############################################################################
# Compress and index reference variants
bgzip ref.vcf
tabix -p vcf ref.vcf.gz

# Create targets file from reference variants
bcftools query -f '%CHROM\t%POS\t%REF,%ALT\n' ref.vcf.gz | \
    bgzip -c > targets.tsv.gz && tabix -s1 -b2 targets.tsv.gz

gunzip ref.vcf.gz

# Merge bams and run mpileup on all bams at once
# Note: there is no difference between running on merged bams or individually
# Specify targets file, and for output specify depth at each locus (-a DP)
samtools merge merged.bam intermediates/*_rg.bam
bcftools mpileup -d1000000 --skip-indels -Ov merged.bam \
    --fasta-ref $FASTA_REFERENCE -a FORMAT/DP \
    --targets-file targets.tsv.gz -o merged.mpileup 

# Calling variants with -v parameter will remove null alleles
bcftools call merged.mpileup -v -f GQ -t targets.tsv.gz -O v -m > merged.vcf

# Some additional quality filtering:
bcftools filter -S . -e "FMT/GQ<20" merged.vcf -o merged.filt.vcf

###############################################################################
# 3b. merge vcfs, extract genotype and depth information
###############################################################################
bgzip merged.mpileup ; tabix -p vcf merged.mpileup.gz ; gunzip merged.mpileup
bgzip merged.vcf ; tabix -p vcf merged.vcf.gz ; gunzip merged.vcf
bgzip merged.filt.vcf ; tabix -p vcf merged.filt.vcf.gz ; gunzip merged.filt.vcf

# Merge all samples, clean up sample names, then sort
#bcftools merge outputs/*.vcf.gz -Ov -o nugen.vcf
#sed -i 's;intermediates/;;g' nugen.vcf
#sed -i 's/_L001_R1_001.bam//g' nugen.vcf

#./vcfsort.sh nugen.vcf > nugen.sort.vcf
./extra_scripts/vcfsort.sh merged.mpileup > merged.sort.mpileup
./extra_scripts/vcfsort.sh merged.vcf > merged.sort.vcf
./extra_scripts/vcfsort.sh merged.filt.vcf > merged.filt.sort.vcf
rm merged.mpileup ; rm merged.vcf ; rm merged.filt.vcf

# Just some text formatting
sed -i 's;intermediates/;;g' merged.sort.mpileup
sed -i 's/_L001_R1_001.bam//g' merged.sort.mpileup
sed -i 's;intermediates/;;g' merged.sort.vcf
sed -i 's/_L001_R1_001.bam//g' merged.sort.vcf
sed -i 's;intermediates/;;g' merged.filt.sort.vcf
sed -i 's/_L001_R1_001.bam//g' merged.filt.sort.vcf

# Get output statistics: depth table for each snp/sample (for updated names + original names)
gatk VariantsToTable -V merged.sort.mpileup -F CHROM -F POS \
    -F REF -F ALT -GF DP -O merged.depths.mpileup.table

gatk VariantsToTable -V merged.sort.vcf -F CHROM -F POS \
    -F REF -F ALT -GF DP -O merged.depths.table

# For the final table, get the genotypes as well 
gatk VariantsToTable -V merged.filt.sort.vcf -F CHROM -F POS \
    -F REF -F ALT -GF DP -O merged.filt.depths.table
gatk VariantsToTable -V merged.filt.sort.vcf -F CHROM -F POS \
    -F REF -F ALT -GF GT -O merged.filt.genotypes.table
    
# do the same for reference vcf genotypes
gatk VariantsToTable -V ref.vcf -F CHROM -F POS -F REF -F ALT -GF GT \
    -O ref.genotypes.table

gatk VariantsToTable -V ref.vcf -F CHROM -F POS -F REF -F ALT -GF DP \
    -O ref.depths.table


###############################################################################
# 4. Run Rscript to clean output files 
# & then remove unnecessary files (2nd is optional)
###############################################################################
cd extra_scripts ; Rscript tableCleaningAndFormatting.R ; cd ..

mv *table data
mv sampleStatistics.csv data

rm -rf intermediates
rm targets*
rm ref.vcf*
rm merged.*
    
conda deactivate
