###############################################################################
# Nugen pipeline for calling variants from targeted sequencing assay
# Abbreviated version: raw .fastq to .vcf
###############################################################################

conda activate variant_analysis

# Requires:
#   - paired end .fastqreads ending with "R1_001.fastq.gz" for fwd reads
#   in a directory called "seqs"
#   - a .vcf file specifying target sites ; generated from plink .map and .ped using
#   "plink --file <file_name> --horse --recode vcf"
#   - EquCab2.0 reference genome with Bowtie2 formatted files + whole genome fasta file

# Specify reference genome and reference genome build for bowtie2
BOWTIE2_REFERENCE="/home/stefan/Documents/Research/msc/EquCab2/Sequence/Bowtie2Index/genome"
FASTA_REFERENCE="/home/stefan/Documents/Research/msc/EquCab2/Sequence/WholeGenomeFasta/genome.fa"

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

# Merge bams and run mpileup on all bams at once for convenience
# Specify targets file, and for output specify depth at each locus (-a DP)
samtools merge merged.bam intermediates/*_rg.bam
bcftools mpileup -d1000000 --skip-indels -Ov merged.bam \
    --fasta-ref $FASTA_REFERENCE -a FORMAT/DP \
    --targets-file targets.tsv.gz -o merged.mpileup 

# Calling variants with -v parameter will remove null alleles
bcftools call merged.mpileup -v -f GQ -t targets.tsv.gz -O v -m > merged.vcf

# Some additional quality filtering:
bcftools filter -S . -e "FMT/GQ<20" merged.vcf -o merged.filt.vcf
    
conda deactivate
