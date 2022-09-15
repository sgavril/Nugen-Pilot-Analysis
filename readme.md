Allegro Targeted Genotyping variant calling
================

# Overview

This workflow calls variants using sequence data obtained from the
Allegro (formerly Nugen) Targeted Sequencing (ATG) assay. This repo
contains the code used to call variants in our publication 
available at https://link.springer.com/article/10.1007/s12686-022-01259-2.
The repo is kept public here for posterity. 

### Set up

First, download miniconda3 (documentation available at
<https://docs.conda.io/en/latest/miniconda.html>). Then create and
activate the environment using the commands below.

``` zsh
conda env create --file environment.yml
conda activate variant_analysis
```

Next we will need to download the EquCab2 reference genome. I am sure an
online version could be accessed, but I found it handy to have it
available locally. The file is large (about 7.7g compressed), so beware.
The file can be downloaded using:

``` zsh
wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Equus_caballus/Ensembl/EquCab2/Equus_caballus_Ensembl_EquCab2.tar.gz
tar -xzvf Equus_caballus_Ensembl_EquCab2.tar.gz
```

And lastly, before getting started with the workflow, I will create a
VCF file of the reference genotype files (which are currently in PLINK
format), along with changing chromosome 32 to X and sorting the VCF file
using a script I found on github.

``` zsh
plink --file data/SNPchip_genotypes_in_Plink_format/Sable_Nugen_Illumina_Affy_combined \
    --horse --recode vcf # makes plink.vcf
# Change chromosome 32 to X, then sort VCF by chromosome & position
sed 's/^32/X/g' plink.vcf | ./extra_scripts/vcfsort.sh - > ref.vcf
rm plink.* # remove unnecessary files from plink
```

    ## PLINK v1.90b6.18 64-bit (16 Jun 2020)          www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to plink.log.
    ## Options in effect:
    ##   --file data/SNPchip_genotypes_in_Plink_format/Sable_Nugen_Illumina_Affy_combined
    ##   --horse
    ##   --recode vcf
    ## 
    ## 15427 MB RAM detected; reserving 7713 MB for main workspace.
    ## Scanning .ped file... 0%2%4%6%9%11%13%15%18%20%22%24%27%29%31%34%36%38%40%43%45%47%49%52%54%56%59%61%63%65%68%70%72%74%77%79%81%84%86%88%90%93%95%97%100%.ped scan complete (for binary autoconversion).
    ## Performing single-pass .bed write (279 variants, 44 horses).
    ## 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%--file: plink-temporary.bed + plink-temporary.bim + plink-temporary.fam
    ## written.
    ## 279 variants loaded from .bim file.
    ## 44 horses (14 males, 14 females, 16 ambiguous) loaded from .fam.
    ## Ambiguous sex IDs written to plink.nosex .
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 44 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Warning: 32 het. haploid genotypes present (see plink.hh ); many commands treat
    ## these as missing.
    ## Total genotyping rate is 0.923265.
    ## 279 variants and 44 horses pass filters and QC.
    ## Note: No phenotypes present.
    ## Warning: Underscore(s) present in sample IDs.
    ## --recode vcf to plink.vcf ... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.

I also want to create a file that keeps track of read counts throughout
the pipeline, starting with raw read counts. I remove the sample
statistics file in case it already exists, since I will append counts to
it during each step.

``` bash
rm sampleStatistics.csv
for sample in seqs/*fastq.gz
do
    base=`basename -s "_L001_R1_001.fastq.gz" $sample`
    zcat $sample | echo $base $((`wc -l`/4)) >> sampleStatistics.csv
done
```

### TrimGalore

We must remove the first 40 base pairs from our sequences as this
contains sequences that are not of biological relevance. For this we use
trim_galore, and we relaxed some of the parameters to try and increase
the number of reads retained, but this does not have a large effect.
Below are the first 100 lines of output which summarize the effects of
running trim_galore on the first input file:

``` zsh
for sample in seqs/*R1_001.fastq.gz; do
    echo $sample
    base=`basename -s "_L001_R1_001.fastq.gz" $sample`
    echo $base
    trim_galore --length 15 --stringency 7 --clip_R1 40 \
        $sample -o intermediates
done
```

    ## seqs/B1_2014_S3_L001_R1_001.fastq.gz
    ## B1_2014_S3
    ## Multicore support not enabled. Proceeding with single-core trimming.
    ## Path to Cutadapt set as: 'cutadapt' (default)
    ## Cutadapt seems to be working fine (tested command 'cutadapt --version')
    ## Cutadapt version: 3.7
    ## single-core operation.
    ## No quality encoding type selected. Assuming that the data provided uses Sanger encoded Phred scores (default)
    ## 
    ## Output will be written into the directory: /home/stefan/Documents/Research/msc/analyses/NugenPilotAnalysis/intermediates/
    ## 
    ## 
    ## AUTO-DETECTING ADAPTER TYPE
    ## ===========================
    ## Attempting to auto-detect adapter type from the first 1 million sequences of the first file (>> seqs/B1_2014_S3_L001_R1_001.fastq.gz <<)
    ## 
    ## Found perfect matches for the following adapter sequences:
    ## Adapter type Count   Sequence    Sequences analysed  Percentage
    ## Illumina 1   AGATCGGAAGAGC   9814    0.01
    ## smallRNA 0   TGGAATTCTCGG    9814    0.00
    ## Nextera  0   CTGTCTCTTATA    9814    0.00
    ## Using Illumina adapter for trimming (count: 1). Second best hit was smallRNA (count: 0)
    ## 
    ## Writing report to '/home/stefan/Documents/Research/msc/analyses/NugenPilotAnalysis/intermediates/B1_2014_S3_L001_R1_001.fastq.gz_trimming_report.txt'
    ## 
    ## SUMMARISING RUN PARAMETERS
    ## ==========================
    ## Input filename: seqs/B1_2014_S3_L001_R1_001.fastq.gz
    ## Trimming mode: single-end
    ## Trim Galore version: 0.6.2
    ## Cutadapt version: 3.7
    ## Number of cores used for trimming: 1
    ## Quality Phred score cutoff: 20
    ## Quality encoding type selected: ASCII+33
    ## Adapter sequence: 'AGATCGGAAGAGC' (Illumina TruSeq, Sanger iPCR; auto-detected)
    ## Maximum trimming error rate: 0.1 (default)
    ## Minimum required adapter overlap (stringency): 7 bp
    ## Minimum required sequence length before a sequence gets removed: 15 bp
    ## All Read 1 sequences will be trimmed by 40 bp from their 5' end to avoid poor qualities or biases
    ## Output file(s) will be GZIP compressed
    ## 
    ## Cutadapt seems to be fairly up-to-date (version 3.7). Setting -j 1
    ## Writing final adapter and quality trimmed output to B1_2014_S3_L001_R1_001_trimmed.fq.gz
    ## 
    ## 
    ##   >>> Now performing quality (cutoff '-q 20') and adapter trimming in a single pass for the adapter sequence: 'AGATCGGAAGAGC' from file seqs/B1_2014_S3_L001_R1_001.fastq.gz <<< 
    ## This is cutadapt 3.7 with Python 3.9.9
    ## Command line parameters: -j 1 -e 0.1 -q 20 -O 7 -a AGATCGGAAGAGC seqs/B1_2014_S3_L001_R1_001.fastq.gz
    ## Processing reads on 1 core in single-end mode ...
    ## Finished in 0.23 s (24 µs/read; 2.54 M reads/minute).
    ## 
    ## === Summary ===
    ## 
    ## Total reads processed:                   9,814
    ## Reads with adapters:                         3 (0.0%)
    ## Reads written (passing filters):         9,814 (100.0%)
    ## 
    ## Total basepairs processed:     1,315,984 bp
    ## Quality-trimmed:                   3,854 bp (0.3%)
    ## Total written (filtered):      1,312,029 bp (99.7%)
    ## 
    ## === Adapter 1 ===
    ## 
    ## Sequence: AGATCGGAAGAGC; Type: regular 3'; Length: 13; Trimmed: 3 times
    ## 
    ## Minimum overlap: 7
    ## No. of allowed errors:
    ## 1-9 bp: 0; 10-13 bp: 1
    ## 
    ## Bases preceding removed adapters:
    ##   A: 0.0%
    ##   C: 33.3%
    ##   G: 33.3%
    ##   T: 33.3%
    ##   none/other: 0.0%
    ## 
    ## Overview of removed sequences
    ## length   count   expect  max.err error counts
    ## 20   1   0.0 1   0 1
    ## 27   1   0.0 1   0 1
    ## 54   1   0.0 1   1
    ## 
    ## RUN STATISTICS FOR INPUT FILE: seqs/B1_2014_S3_L001_R1_001.fastq.gz
    ## =============================================
    ## 9814 sequences processed in total
    ## Sequences removed because they became shorter than the length cutoff of 15 bp:   192 (2.0%)
    ## 
    ## seqs/B28_2014_S4_L001_R1_001.fastq.gz
    ## B28_2014_S4
    ## Multicore support not enabled. Proceeding with single-core trimming.
    ## Path to Cutadapt set as: 'cutadapt' (default)
    ## Cutadapt seems to be working fine (tested command 'cutadapt --version')
    ## Cutadapt version: 3.7
    ## single-core operation.
    ## No quality encoding type selected. Assuming that the data provided uses Sanger encoded Phred scores (default)
    ## 
    ## Output will be written into the directory: /home/stefan/Documents/Research/msc/analyses/NugenPilotAnalysis/intermediates/
    ## 
    ## 
    ## AUTO-DETECTING ADAPTER TYPE
    ...

Then append these read counts to our sample statistics file by first
writing them to a temporary text file.

``` bash
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
```

    ## rm: cannot remove 'tmp.txt': No such file or directory

### Single end read alignment

Next we will align the forward reads to the EquCab2 reference genome. We
have enabled very sensitive local alignment as it appears to retain a
similar amount of reads compared to the less sensitive parameters. Then
we move the .bam outputs to the intermediate file folder generated by
trim_galore. I start by specifying where the reference genome is
located. Specifically, I refer to the bowtie2 reference genome.

``` zsh
BOWTIE2_REFERENCE="/home/stefan/Documents/Research/msc/EquCab2/Sequence/Bowtie2Index/genome"

for sample in intermediates/*trimmed.fq.gz ; do
    echo $sample
    base=`basename -s "_trimmed.fq.gz" $sample`
    echo $base
    bowtie2 -p 2 --very-sensitive-local -U $sample \
        -x $BOWTIE2_REFERENCE | samtools sort -o $base".bam"
done

mv *.bam intermediates
```

    ## intermediates/B1_2014_S3_L001_R1_001_trimmed.fq.gz
    ## B1_2014_S3_L001_R1_001
    ## 9622 reads; of these:
    ##   9622 (100.00%) were unpaired; of these:
    ##     863 (8.97%) aligned 0 times
    ##     7522 (78.18%) aligned exactly 1 time
    ##     1237 (12.86%) aligned >1 times
    ## 91.03% overall alignment rate
    ## intermediates/B28_2014_S4_L001_R1_001_trimmed.fq.gz
    ## B28_2014_S4_L001_R1_001
    ## 7786 reads; of these:
    ##   7786 (100.00%) were unpaired; of these:
    ##     1285 (16.50%) aligned 0 times
    ##     5491 (70.52%) aligned exactly 1 time
    ##     1010 (12.97%) aligned >1 times
    ## 83.50% overall alignment rate
    ## intermediates/B33_2014_S5_L001_R1_001_trimmed.fq.gz
    ## B33_2014_S5_L001_R1_001
    ## 22243 reads; of these:
    ##   22243 (100.00%) were unpaired; of these:
    ##     2116 (9.51%) aligned 0 times
    ##     17450 (78.45%) aligned exactly 1 time
    ##     2677 (12.04%) aligned >1 times
    ## 90.49% overall alignment rate
    ## intermediates/B36_2014_S6_L001_R1_001_trimmed.fq.gz
    ## B36_2014_S6_L001_R1_001
    ## 14741 reads; of these:
    ##   14741 (100.00%) were unpaired; of these:
    ##     889 (6.03%) aligned 0 times
    ##     12090 (82.02%) aligned exactly 1 time
    ##     1762 (11.95%) aligned >1 times
    ## 93.97% overall alignment rate
    ## intermediates/B39_2014_S7_L001_R1_001_trimmed.fq.gz
    ## B39_2014_S7_L001_R1_001
    ## 153 reads; of these:
    ##   153 (100.00%) were unpaired; of these:
    ##     16 (10.46%) aligned 0 times
    ##     116 (75.82%) aligned exactly 1 time
    ##     21 (13.73%) aligned >1 times
    ## 89.54% overall alignment rate
    ## intermediates/E1_2014_S8_L001_R1_001_trimmed.fq.gz
    ## E1_2014_S8_L001_R1_001
    ## 24153 reads; of these:
    ##   24153 (100.00%) were unpaired; of these:
    ##     3126 (12.94%) aligned 0 times
    ##     18835 (77.98%) aligned exactly 1 time
    ##     2192 (9.08%) aligned >1 times
    ## 87.06% overall alignment rate
    ## intermediates/E12_2014_S13_L001_R1_001_trimmed.fq.gz
    ## E12_2014_S13_L001_R1_001
    ## 21887 reads; of these:
    ##   21887 (100.00%) were unpaired; of these:
    ##     2024 (9.25%) aligned 0 times
    ##     17426 (79.62%) aligned exactly 1 time
    ##     2437 (11.13%) aligned >1 times
    ## 90.75% overall alignment rate
    ## intermediates/E13_2014_S14_L001_R1_001_trimmed.fq.gz
    ## E13_2014_S14_L001_R1_001
    ## 14424 reads; of these:
    ##   14424 (100.00%) were unpaired; of these:
    ##     1268 (8.79%) aligned 0 times
    ##     11429 (79.24%) aligned exactly 1 time
    ##     1727 (11.97%) aligned >1 times
    ## 91.21% overall alignment rate
    ## intermediates/E14_2014_S15_L001_R1_001_trimmed.fq.gz
    ## E14_2014_S15_L001_R1_001
    ## 23646 reads; of these:
    ##   23646 (100.00%) were unpaired; of these:
    ##     1155 (4.88%) aligned 0 times
    ##     19827 (83.85%) aligned exactly 1 time
    ##     2664 (11.27%) aligned >1 times
    ## 95.12% overall alignment rate
    ## intermediates/E15_2014_S16_L001_R1_001_trimmed.fq.gz
    ## E15_2014_S16_L001_R1_001
    ## 14926 reads; of these:
    ##   14926 (100.00%) were unpaired; of these:
    ##     1039 (6.96%) aligned 0 times
    ##     12236 (81.98%) aligned exactly 1 time
    ##     1651 (11.06%) aligned >1 times
    ## 93.04% overall alignment rate
    ## intermediates/E17_2014_S17_L001_R1_001_trimmed.fq.gz
    ## E17_2014_S17_L001_R1_001
    ## 28041 reads; of these:
    ##   28041 (100.00%) were unpaired; of these:
    ##     1594 (5.68%) aligned 0 times
    ##     23032 (82.14%) aligned exactly 1 time
    ##     3415 (12.18%) aligned >1 times
    ## 94.32% overall alignment rate
    ## intermediates/E18_2014_S18_L001_R1_001_trimmed.fq.gz
    ## E18_2014_S18_L001_R1_001
    ## 12592 reads; of these:
    ##   12592 (100.00%) were unpaired; of these:
    ##     691 (5.49%) aligned 0 times
    ##     10364 (82.31%) aligned exactly 1 time
    ##     1537 (12.21%) aligned >1 times
    ## 94.51% overall alignment rate
    ## intermediates/E20_2014_S19_L001_R1_001_trimmed.fq.gz
    ## E20_2014_S19_L001_R1_001
    ## 13500 reads; of these:
    ##   13500 (100.00%) were unpaired; of these:
    ...

Now we get the read counts again.

``` bash
rm tmp.map.txt
# Number of reads that mapped to reference genome (EquCab2)
for sample in intermediates/*.bam; do
    base=`basename -s ".bam" $sample`
    echo "`samtools view -c -F 260 $sample` " >> tmp.map.txt
done

paste -d " " sampleStatistics.csv tmp.map.txt > tmp.txt
mv tmp.txt sampleStatistics.csv
rm tmp.map.txt
```

    ## rm: cannot remove 'tmp.map.txt': No such file or directory

Something that has caused a lot of headaches for me in the past has been
the necessity of adding read group information to bam files. Without
this step, information like the sample name is lost in the .bam file
which creates a lot of problems downstream, so we can use the GATK tool
“AddOrReplaceReadGroups” here. Again, we move the outputs to the
intermediate file folder.

``` zsh
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
```

    ## intermediates/B1_2014_S3_L001_R1_001.bam
    ## B1_2014_S3_L001_R1_001
    ## Using GATK jar /home/stefan/miniconda3/envs/variant_analysis/share/gatk4-4.1.9.0-0/gatk-package-4.1.9.0-local.jar
    ## Running:
    ##     java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /home/stefan/miniconda3/envs/variant_analysis/share/gatk4-4.1.9.0-0/gatk-package-4.1.9.0-local.jar AddOrReplaceReadGroups I=intermediates/B1_2014_S3_L001_R1_001.bam O=B1_2014_S3_L001_R1_001_rg.bam RGID=1 RGLB=lib1 RGPL=nugen RGPU=unit1 RGSM=intermediates/B1_2014_S3_L001_R1_001.bam
    ## INFO 2022-03-02 14:21:16 AddOrReplaceReadGroups  
    ## 
    ## ********** NOTE: Picard's command line syntax is changing.
    ## **********
    ## ********** For more information, please see:
    ## ********** https://github.com/broadinstitute/picard/wiki/Command-Line-Syntax-Transition-For-Users-(Pre-Transition)
    ## **********
    ## ********** The command line looks like this in the new syntax:
    ## **********
    ## **********    AddOrReplaceReadGroups -I intermediates/B1_2014_S3_L001_R1_001.bam -O B1_2014_S3_L001_R1_001_rg.bam -RGID 1 -RGLB lib1 -RGPL nugen -RGPU unit1 -RGSM intermediates/B1_2014_S3_L001_R1_001.bam
    ## **********
    ## 
    ## 
    ## 14:21:16.965 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/home/stefan/miniconda3/envs/variant_analysis/share/gatk4-4.1.9.0-0/gatk-package-4.1.9.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    ## [Wed Mar 02 14:21:16 MST 2022] AddOrReplaceReadGroups INPUT=intermediates/B1_2014_S3_L001_R1_001.bam OUTPUT=B1_2014_S3_L001_R1_001_rg.bam RGID=1 RGLB=lib1 RGPL=nugen RGPU=unit1 RGSM=intermediates/B1_2014_S3_L001_R1_001.bam    VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=2 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false GA4GH_CLIENT_SECRETS=client_secrets.json USE_JDK_DEFLATER=false USE_JDK_INFLATER=false
    ## Mar 02, 2022 2:21:17 PM shaded.cloud_nio.com.google.auth.oauth2.ComputeEngineCredentials runningOnComputeEngine
    ## INFO: Failed to detect whether we are running on Google Compute Engine.
    ## [Wed Mar 02 14:21:17 MST 2022] Executing as stefan@stefan-81x2 on Linux 5.15.25-1-MANJARO amd64; OpenJDK 64-Bit Server VM 1.8.0_312-b07; Deflater: Intel; Inflater: Intel; Provider GCS is available; Picard version: 4.1.9.0
    ## INFO 2022-03-02 14:21:17 AddOrReplaceReadGroups  Created read-group ID=1 PL=nugen LB=lib1 SM=intermediates/B1_2014_S3_L001_R1_001.bam
    ## 
    ## [Wed Mar 02 14:21:17 MST 2022] picard.sam.AddOrReplaceReadGroups done. Elapsed time: 0.01 minutes.
    ## Runtime.totalMemory()=342884352
    ## Tool returned:
    ## 0
    ## intermediates/B28_2014_S4_L001_R1_001.bam
    ## B28_2014_S4_L001_R1_001
    ## Using GATK jar /home/stefan/miniconda3/envs/variant_analysis/share/gatk4-4.1.9.0-0/gatk-package-4.1.9.0-local.jar
    ## Running:
    ##     java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /home/stefan/miniconda3/envs/variant_analysis/share/gatk4-4.1.9.0-0/gatk-package-4.1.9.0-local.jar AddOrReplaceReadGroups I=intermediates/B28_2014_S4_L001_R1_001.bam O=B28_2014_S4_L001_R1_001_rg.bam RGID=1 RGLB=lib1 RGPL=nugen RGPU=unit1 RGSM=intermediates/B28_2014_S4_L001_R1_001.bam
    ## INFO 2022-03-02 14:21:19 AddOrReplaceReadGroups  
    ## 
    ## ********** NOTE: Picard's command line syntax is changing.
    ## **********
    ## ********** For more information, please see:
    ## ********** https://github.com/broadinstitute/picard/wiki/Command-Line-Syntax-Transition-For-Users-(Pre-Transition)
    ## **********
    ## ********** The command line looks like this in the new syntax:
    ## **********
    ## **********    AddOrReplaceReadGroups -I intermediates/B28_2014_S4_L001_R1_001.bam -O B28_2014_S4_L001_R1_001_rg.bam -RGID 1 -RGLB lib1 -RGPL nugen -RGPU unit1 -RGSM intermediates/B28_2014_S4_L001_R1_001.bam
    ## **********
    ## 
    ## 
    ## 14:21:19.270 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/home/stefan/miniconda3/envs/variant_analysis/share/gatk4-4.1.9.0-0/gatk-package-4.1.9.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    ## [Wed Mar 02 14:21:19 MST 2022] AddOrReplaceReadGroups INPUT=intermediates/B28_2014_S4_L001_R1_001.bam OUTPUT=B28_2014_S4_L001_R1_001_rg.bam RGID=1 RGLB=lib1 RGPL=nugen RGPU=unit1 RGSM=intermediates/B28_2014_S4_L001_R1_001.bam    VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=2 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false GA4GH_CLIENT_SECRETS=client_secrets.json USE_JDK_DEFLATER=false USE_JDK_INFLATER=false
    ## Mar 02, 2022 2:21:19 PM shaded.cloud_nio.com.google.auth.oauth2.ComputeEngineCredentials runningOnComputeEngine
    ## INFO: Failed to detect whether we are running on Google Compute Engine.
    ## [Wed Mar 02 14:21:19 MST 2022] Executing as stefan@stefan-81x2 on Linux 5.15.25-1-MANJARO amd64; OpenJDK 64-Bit Server VM 1.8.0_312-b07; Deflater: Intel; Inflater: Intel; Provider GCS is available; Picard version: 4.1.9.0
    ## INFO 2022-03-02 14:21:19 AddOrReplaceReadGroups  Created read-group ID=1 PL=nugen LB=lib1 SM=intermediates/B28_2014_S4_L001_R1_001.bam
    ## 
    ## [Wed Mar 02 14:21:19 MST 2022] picard.sam.AddOrReplaceReadGroups done. Elapsed time: 0.01 minutes.
    ## Runtime.totalMemory()=347078656
    ## Tool returned:
    ## 0
    ## intermediates/B33_2014_S5_L001_R1_001.bam
    ## B33_2014_S5_L001_R1_001
    ## Using GATK jar /home/stefan/miniconda3/envs/variant_analysis/share/gatk4-4.1.9.0-0/gatk-package-4.1.9.0-local.jar
    ## Running:
    ##     java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /home/stefan/miniconda3/envs/variant_analysis/share/gatk4-4.1.9.0-0/gatk-package-4.1.9.0-local.jar AddOrReplaceReadGroups I=intermediates/B33_2014_S5_L001_R1_001.bam O=B33_2014_S5_L001_R1_001_rg.bam RGID=1 RGLB=lib1 RGPL=nugen RGPU=unit1 RGSM=intermediates/B33_2014_S5_L001_R1_001.bam
    ## INFO 2022-03-02 14:21:21 AddOrReplaceReadGroups  
    ## 
    ## ********** NOTE: Picard's command line syntax is changing.
    ## **********
    ## ********** For more information, please see:
    ## ********** https://github.com/broadinstitute/picard/wiki/Command-Line-Syntax-Transition-For-Users-(Pre-Transition)
    ## **********
    ## ********** The command line looks like this in the new syntax:
    ## **********
    ## **********    AddOrReplaceReadGroups -I intermediates/B33_2014_S5_L001_R1_001.bam -O B33_2014_S5_L001_R1_001_rg.bam -RGID 1 -RGLB lib1 -RGPL nugen -RGPU unit1 -RGSM intermediates/B33_2014_S5_L001_R1_001.bam
    ## **********
    ## 
    ## 
    ## 14:21:21.676 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/home/stefan/miniconda3/envs/variant_analysis/share/gatk4-4.1.9.0-0/gatk-package-4.1.9.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    ## [Wed Mar 02 14:21:21 MST 2022] AddOrReplaceReadGroups INPUT=intermediates/B33_2014_S5_L001_R1_001.bam OUTPUT=B33_2014_S5_L001_R1_001_rg.bam RGID=1 RGLB=lib1 RGPL=nugen RGPU=unit1 RGSM=intermediates/B33_2014_S5_L001_R1_001.bam    VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=2 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false GA4GH_CLIENT_SECRETS=client_secrets.json USE_JDK_DEFLATER=false USE_JDK_INFLATER=false
    ## Mar 02, 2022 2:21:22 PM shaded.cloud_nio.com.google.auth.oauth2.ComputeEngineCredentials runningOnComputeEngine
    ## INFO: Failed to detect whether we are running on Google Compute Engine.
    ## [Wed Mar 02 14:21:22 MST 2022] Executing as stefan@stefan-81x2 on Linux 5.15.25-1-MANJARO amd64; OpenJDK 64-Bit Server VM 1.8.0_312-b07; Deflater: Intel; Inflater: Intel; Provider GCS is available; Picard version: 4.1.9.0
    ## INFO 2022-03-02 14:21:22 AddOrReplaceReadGroups  Created read-group ID=1 PL=nugen LB=lib1 SM=intermediates/B33_2014_S5_L001_R1_001.bam
    ## 
    ## [Wed Mar 02 14:21:22 MST 2022] picard.sam.AddOrReplaceReadGroups done. Elapsed time: 0.01 minutes.
    ## Runtime.totalMemory()=343408640
    ## Tool returned:
    ## 0
    ## intermediates/B36_2014_S6_L001_R1_001.bam
    ## B36_2014_S6_L001_R1_001
    ## Using GATK jar /home/stefan/miniconda3/envs/variant_analysis/share/gatk4-4.1.9.0-0/gatk-package-4.1.9.0-local.jar
    ## Running:
    ##     java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /home/stefan/miniconda3/envs/variant_analysis/share/gatk4-4.1.9.0-0/gatk-package-4.1.9.0-local.jar AddOrReplaceReadGroups I=intermediates/B36_2014_S6_L001_R1_001.bam O=B36_2014_S6_L001_R1_001_rg.bam RGID=1 RGLB=lib1 RGPL=nugen RGPU=unit1 RGSM=intermediates/B36_2014_S6_L001_R1_001.bam
    ## INFO 2022-03-02 14:21:24 AddOrReplaceReadGroups  
    ## 
    ## ********** NOTE: Picard's command line syntax is changing.
    ## **********
    ## ********** For more information, please see:
    ## ********** https://github.com/broadinstitute/picard/wiki/Command-Line-Syntax-Transition-For-Users-(Pre-Transition)
    ## **********
    ## ********** The command line looks like this in the new syntax:
    ...

### Call variants using mpileup

Now we can call variants, but since we have a set of variants that act
as our ‘reference data set’, we can create a list to input to variant
calling so that mpileup only outputs these particular sites. To do that,
we have to compress and index the vcf file we generated at the start of
the workflow. Then we create the targets file, and I like to have the
file decompressed so I unzip it again.

``` zsh
bgzip ref.vcf
tabix -p vcf ref.vcf.gz

# Create targets file from reference variants
bcftools query -f '%CHROM\t%POS\t%REF,%ALT\n' ref.vcf.gz | \
    bgzip -c > targets.tsv.gz && tabix -s1 -b2 targets.tsv.gz

gunzip ref.vcf.gz
```

Now we can call the variants:

``` zsh
FASTA_REFERENCE="/home/stefan/Documents/Research/msc/EquCab2/Sequence/WholeGenomeFasta/genome.fa"
samtools merge merged.bam intermediates/*_rg.bam
bcftools mpileup -d1000000 --skip-indels -Ov merged.bam \
    --fasta-ref $FASTA_REFERENCE -a FORMAT/DP \
    --targets-file targets.tsv.gz -o merged.mpileup 
```

    ## [mpileup] 50 samples in 1 input files
    ## [mpileup] maximum number of reads per input file set to -d 1000000

For some quality control, we can remove null alleles by using the -v
parameter, and then

``` zsh
bcftools call merged.mpileup -v -f GQ -t targets.tsv.gz -O v -m > merged.vcf

bcftools filter -S . -e "FMT/GQ<20" merged.vcf -o merged.filt.vcf
```

    ## Note: none of --samples-file, --ploidy or --ploidy-file given, assuming all sites are diploid

### Create genotype and depth count tables

Now we can use the utility from GATK to create nicely formatted tables
detailing the genotypes and read depth for each genotype. First I index
the vcf files again (for use with GATK), sort all output VCFs and then
do some quick text formatting with sed.

``` zsh
bgzip merged.mpileup ; tabix -p vcf merged.mpileup.gz ; gunzip merged.mpileup
bgzip merged.vcf ; tabix -p vcf merged.vcf.gz ; gunzip merged.vcf
bgzip merged.filt.vcf ; tabix -p vcf merged.filt.vcf.gz ; gunzip merged.filt.vcf

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
```

And then I get all of the desired output tables.

``` zsh
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
```

    ## Using GATK jar /home/stefan/miniconda3/envs/variant_analysis/share/gatk4-4.1.9.0-0/gatk-package-4.1.9.0-local.jar
    ## Running:
    ##     java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /home/stefan/miniconda3/envs/variant_analysis/share/gatk4-4.1.9.0-0/gatk-package-4.1.9.0-local.jar VariantsToTable -V merged.sort.mpileup -F CHROM -F POS -F REF -F ALT -GF DP -O merged.depths.mpileup.table
    ## 14:24:07.093 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/home/stefan/miniconda3/envs/variant_analysis/share/gatk4-4.1.9.0-0/gatk-package-4.1.9.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    ## Mar 02, 2022 2:24:07 PM shaded.cloud_nio.com.google.auth.oauth2.ComputeEngineCredentials runningOnComputeEngine
    ## INFO: Failed to detect whether we are running on Google Compute Engine.
    ## 14:24:07.244 INFO  VariantsToTable - ------------------------------------------------------------
    ## 14:24:07.244 INFO  VariantsToTable - The Genome Analysis Toolkit (GATK) v4.1.9.0
    ## 14:24:07.244 INFO  VariantsToTable - For support and documentation go to https://software.broadinstitute.org/gatk/
    ## 14:24:07.245 INFO  VariantsToTable - Executing as stefan@stefan-81x2 on Linux v5.15.25-1-MANJARO amd64
    ## 14:24:07.245 INFO  VariantsToTable - Java runtime: OpenJDK 64-Bit Server VM v1.8.0_312-b07
    ## 14:24:07.245 INFO  VariantsToTable - Start Date/Time: March 2, 2022 2:24:07 MST PM
    ## 14:24:07.245 INFO  VariantsToTable - ------------------------------------------------------------
    ## 14:24:07.245 INFO  VariantsToTable - ------------------------------------------------------------
    ## 14:24:07.245 INFO  VariantsToTable - HTSJDK Version: 2.23.0
    ## 14:24:07.245 INFO  VariantsToTable - Picard Version: 2.23.3
    ## 14:24:07.245 INFO  VariantsToTable - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    ## 14:24:07.245 INFO  VariantsToTable - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    ## 14:24:07.245 INFO  VariantsToTable - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    ## 14:24:07.245 INFO  VariantsToTable - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    ## 14:24:07.245 INFO  VariantsToTable - Deflater: IntelDeflater
    ## 14:24:07.246 INFO  VariantsToTable - Inflater: IntelInflater
    ## 14:24:07.246 INFO  VariantsToTable - GCS max retries/reopens: 20
    ## 14:24:07.246 INFO  VariantsToTable - Requester pays: disabled
    ## 14:24:07.246 INFO  VariantsToTable - Initializing engine
    ## 14:24:07.494 INFO  FeatureManager - Using codec VCFCodec to read file file:///home/stefan/Documents/Research/msc/analyses/NugenPilotAnalysis/merged.sort.mpileup
    ## 14:24:07.510 INFO  VariantsToTable - Done initializing engine
    ## 14:24:07.515 WARN  VariantsToTable - Allele-specific fields will only be split if splitting multi-allelic variants is specified (`--split-multi-allelic` or `-SMA`
    ## 14:24:07.516 INFO  ProgressMeter - Starting traversal
    ## 14:24:07.516 INFO  ProgressMeter -        Current Locus  Elapsed Minutes    Variants Processed  Variants/Minute
    ## 14:24:07.603 INFO  ProgressMeter -             unmapped              0.0                   279         194651.2
    ## 14:24:07.604 INFO  ProgressMeter - Traversal complete. Processed 279 total variants in 0.0 minutes.
    ## 14:24:07.604 INFO  VariantsToTable - Shutting down engine
    ## [March 2, 2022 2:24:07 MST PM] org.broadinstitute.hellbender.tools.walkers.variantutils.VariantsToTable done. Elapsed time: 0.01 minutes.
    ## Runtime.totalMemory()=453509120
    ## Using GATK jar /home/stefan/miniconda3/envs/variant_analysis/share/gatk4-4.1.9.0-0/gatk-package-4.1.9.0-local.jar
    ## Running:
    ##     java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /home/stefan/miniconda3/envs/variant_analysis/share/gatk4-4.1.9.0-0/gatk-package-4.1.9.0-local.jar VariantsToTable -V merged.sort.vcf -F CHROM -F POS -F REF -F ALT -GF DP -O merged.depths.table
    ## 14:24:09.806 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/home/stefan/miniconda3/envs/variant_analysis/share/gatk4-4.1.9.0-0/gatk-package-4.1.9.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    ## Mar 02, 2022 2:24:09 PM shaded.cloud_nio.com.google.auth.oauth2.ComputeEngineCredentials runningOnComputeEngine
    ## INFO: Failed to detect whether we are running on Google Compute Engine.
    ## 14:24:09.955 INFO  VariantsToTable - ------------------------------------------------------------
    ## 14:24:09.956 INFO  VariantsToTable - The Genome Analysis Toolkit (GATK) v4.1.9.0
    ## 14:24:09.956 INFO  VariantsToTable - For support and documentation go to https://software.broadinstitute.org/gatk/
    ## 14:24:09.956 INFO  VariantsToTable - Executing as stefan@stefan-81x2 on Linux v5.15.25-1-MANJARO amd64
    ## 14:24:09.956 INFO  VariantsToTable - Java runtime: OpenJDK 64-Bit Server VM v1.8.0_312-b07
    ## 14:24:09.956 INFO  VariantsToTable - Start Date/Time: March 2, 2022 2:24:09 MST PM
    ## 14:24:09.956 INFO  VariantsToTable - ------------------------------------------------------------
    ## 14:24:09.956 INFO  VariantsToTable - ------------------------------------------------------------
    ## 14:24:09.957 INFO  VariantsToTable - HTSJDK Version: 2.23.0
    ## 14:24:09.957 INFO  VariantsToTable - Picard Version: 2.23.3
    ## 14:24:09.957 INFO  VariantsToTable - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    ## 14:24:09.957 INFO  VariantsToTable - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    ## 14:24:09.957 INFO  VariantsToTable - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    ## 14:24:09.957 INFO  VariantsToTable - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    ## 14:24:09.957 INFO  VariantsToTable - Deflater: IntelDeflater
    ## 14:24:09.957 INFO  VariantsToTable - Inflater: IntelInflater
    ## 14:24:09.957 INFO  VariantsToTable - GCS max retries/reopens: 20
    ## 14:24:09.957 INFO  VariantsToTable - Requester pays: disabled
    ## 14:24:09.957 INFO  VariantsToTable - Initializing engine
    ## 14:24:10.195 INFO  FeatureManager - Using codec VCFCodec to read file file:///home/stefan/Documents/Research/msc/analyses/NugenPilotAnalysis/merged.sort.vcf
    ## 14:24:10.212 INFO  VariantsToTable - Done initializing engine
    ## 14:24:10.217 WARN  VariantsToTable - Allele-specific fields will only be split if splitting multi-allelic variants is specified (`--split-multi-allelic` or `-SMA`
    ## 14:24:10.218 INFO  ProgressMeter - Starting traversal
    ## 14:24:10.218 INFO  ProgressMeter -        Current Locus  Elapsed Minutes    Variants Processed  Variants/Minute
    ## 14:24:10.302 INFO  ProgressMeter -             unmapped              0.0                   269         192142.9
    ## 14:24:10.303 INFO  ProgressMeter - Traversal complete. Processed 269 total variants in 0.0 minutes.
    ## 14:24:10.303 INFO  VariantsToTable - Shutting down engine
    ## [March 2, 2022 2:24:10 MST PM] org.broadinstitute.hellbender.tools.walkers.variantutils.VariantsToTable done. Elapsed time: 0.01 minutes.
    ## Runtime.totalMemory()=471334912
    ## Using GATK jar /home/stefan/miniconda3/envs/variant_analysis/share/gatk4-4.1.9.0-0/gatk-package-4.1.9.0-local.jar
    ## Running:
    ##     java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /home/stefan/miniconda3/envs/variant_analysis/share/gatk4-4.1.9.0-0/gatk-package-4.1.9.0-local.jar VariantsToTable -V merged.filt.sort.vcf -F CHROM -F POS -F REF -F ALT -GF DP -O merged.filt.depths.table
    ## 14:24:12.496 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/home/stefan/miniconda3/envs/variant_analysis/share/gatk4-4.1.9.0-0/gatk-package-4.1.9.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    ## Mar 02, 2022 2:24:12 PM shaded.cloud_nio.com.google.auth.oauth2.ComputeEngineCredentials runningOnComputeEngine
    ## INFO: Failed to detect whether we are running on Google Compute Engine.
    ## 14:24:12.651 INFO  VariantsToTable - ------------------------------------------------------------
    ## 14:24:12.651 INFO  VariantsToTable - The Genome Analysis Toolkit (GATK) v4.1.9.0
    ## 14:24:12.651 INFO  VariantsToTable - For support and documentation go to https://software.broadinstitute.org/gatk/
    ## 14:24:12.651 INFO  VariantsToTable - Executing as stefan@stefan-81x2 on Linux v5.15.25-1-MANJARO amd64
    ## 14:24:12.651 INFO  VariantsToTable - Java runtime: OpenJDK 64-Bit Server VM v1.8.0_312-b07
    ## 14:24:12.652 INFO  VariantsToTable - Start Date/Time: March 2, 2022 2:24:12 MST PM
    ## 14:24:12.652 INFO  VariantsToTable - ------------------------------------------------------------
    ## 14:24:12.652 INFO  VariantsToTable - ------------------------------------------------------------
    ## 14:24:12.652 INFO  VariantsToTable - HTSJDK Version: 2.23.0
    ## 14:24:12.652 INFO  VariantsToTable - Picard Version: 2.23.3
    ## 14:24:12.652 INFO  VariantsToTable - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    ## 14:24:12.652 INFO  VariantsToTable - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    ## 14:24:12.652 INFO  VariantsToTable - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    ## 14:24:12.652 INFO  VariantsToTable - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    ## 14:24:12.652 INFO  VariantsToTable - Deflater: IntelDeflater
    ## 14:24:12.652 INFO  VariantsToTable - Inflater: IntelInflater
    ## 14:24:12.653 INFO  VariantsToTable - GCS max retries/reopens: 20
    ## 14:24:12.653 INFO  VariantsToTable - Requester pays: disabled
    ## 14:24:12.653 INFO  VariantsToTable - Initializing engine
    ## 14:24:12.886 INFO  FeatureManager - Using codec VCFCodec to read file file:///home/stefan/Documents/Research/msc/analyses/NugenPilotAnalysis/merged.filt.sort.vcf
    ## 14:24:12.903 INFO  VariantsToTable - Done initializing engine
    ## 14:24:12.908 WARN  VariantsToTable - Allele-specific fields will only be split if splitting multi-allelic variants is specified (`--split-multi-allelic` or `-SMA`
    ## 14:24:12.909 INFO  ProgressMeter - Starting traversal
    ## 14:24:12.909 INFO  ProgressMeter -        Current Locus  Elapsed Minutes    Variants Processed  Variants/Minute
    ...

### Final formatting and clean-up

Here I run an Rscript I made to clean-up some of the formatting (things
like column names and so on) in the final output files.

``` zsh
cd extra_scripts ; Rscript tableCleaningAndFormatting.R ; cd ..
```

    ## 
    ## Attaching package: ‘dplyr’
    ## 
    ## The following objects are masked from ‘package:stats’:
    ## 
    ##     filter, lag
    ## 
    ## The following objects are masked from ‘package:base’:
    ## 
    ##     intersect, setdiff, setequal, union
    ## 
    ## [1] "SI_hair_152" "SI_hair_471" "SI_hair_516" "SI_hair_516"
    ## Warning message:
    ## In merge.data.frame(nug.long, nug.dep.long, by = c("Var1", "Var2")) :
    ##   column names ‘value.x’, ‘value.y’ are duplicated in the result
    ## Warning message:
    ## In merge.data.frame(nug.long, nug.gen.filt.long, by = c("Var1",  :
    ##   column names ‘value.x’, ‘value.y’ are duplicated in the result
    ## [1] 16345.1
    ## [1] 10113.22
    ## [1] 38817
    ## [1] 154
    ## [1] 97.99368
    ## [1] 0.4413755
    ## [1] 88.57712
    ## [1] 4.62811
    ## [1] 12462.04
    ## [1] 7938.3
    ## [1] 12273.48
    ## [1] 7811.272
    ## [1] 12273.48
    ## [1] 7811.272
    ## [1] 220.2708
    ## [1] 77.61929
    ##   Sample          V1    V2    V3    V4    V5    V6    V7  V8 concordance
    ## 1     B1  B1_2014_S3  9814  9622  8759  7258  7167  7167 248    98.18548
    ## 2    B28 B28_2014_S4  7924  7786  6501  5328  5224  5224 232    98.88393
    ## 3    B33 B33_2014_S5 22763 22243 20127 16872 16534 16534 260    98.26923
    ## 4    B36 B36_2014_S6 15030 14741 13852 11658 11491 11491 251    98.77551
    ## 5    B39 B39_2014_S7   154   153   137   107   105   105   5   100.00000
    ## 6     E1  E1_2014_S8 24727 24153 21027 18193 17845 17845 257    98.39357

Here is some optional clean up to move or remove files we no longer
need.

``` zsh
mv *table data
mv sampleStatistics.csv data

rm -rf intermediates
rm targets*
rm ref.vcf*
rm merged.*
```

And lastly, deactivate the conda environment. This workflow will have
created the genotype + depth tables (now in the “data” folder) necessary
to run the post-processing analysis for this manuscript.

``` zsh
conda deactivate
```
