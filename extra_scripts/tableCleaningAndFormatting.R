library(dplyr)
library(reshape2)
library(gtools)

meta <- read.table("../data/info_about_samples_in_test_run.csv", sep = "\t", header = T)
ref.gen <- read.table("../data/ref.genotypes.table", header = T, sep = "\t")

# Remove unnecessary strings in column names
colnames(ref.gen) <- colnames(ref.gen) %>% 
  gsub("\\.GT", "", .) %>%
  gsub("Sable_", "", .) %>%
  gsub("\\.", "-", .)
#colnames(ref.gen)[which(!colnames(ref.gen) %in% meta$SNP_chip_sample)]

# Map sample names
# We have some repeat samples that will map to multiple individuals
# from the reference genotypes to the Nugen genotypes
# Repeat samples (reference): SI_hair_152, SI_hair_471, SI_hair_516 
meta$SNP.w.rep <- meta$SNP_chip_sample
#which(duplicated(meta$SNP.w.rep)) # Indices of duplicate
meta$SNP.w.rep[which(duplicated(meta$SNP.w.rep))] # Duplicates
meta$SNP.w.rep[which(duplicated(meta$SNP.w.rep))] <- 
  c("SI_hair_152b", "SI_hair_471b", "SI_hair_516b", "SI_hair_516c")
ref.gen$SI_hair_152b <- ref.gen$SI_hair_152
ref.gen$SI_hair_471b <- ref.gen$SI_hair_471
ref.gen$SI_hair_516b <- ref.gen$SI_hair_516
ref.gen$SI_hair_516c <- ref.gen$SI_hair_516

# These are the sample names on Nugen
meta$sample <- paste0(meta$fecal_sample_nugen, "_S", meta$SAMPLE_ID)

# Map the columns (change reference genotype names)
newColNames <- c("CHROM", "POS", "REF", "ALT", 
                 meta$sample[match(names(ref.gen[5:ncol(ref.gen)]), meta$SNP.w.rep)])
names(ref.gen) <- newColNames
rm(newColNames)

# Make sure that sample names are same order as other tables
#ref.gen[order(ref.gen[-c(1:4)]), -c(1:4)]
ref.gen.only <- ref.gen[-c(1:4)]
ref.gen.only.sort <- ref.gen.only[, order(colnames(ref.gen.only))]
ref.gen <- cbind(ref.gen[c(1:4)], ref.gen.only.sort)

rm(ref.gen.only) ; rm(ref.gen.only.sort)

write.table(ref.gen, "../data/ref.geno.formatted.table", quote = F, sep = "\t",
            col.names = T, row.names = F)

sampleStats <- read.table("../data/sampleStatistics.csv", sep = " ")
ref.gen <- read.table("../data/ref.geno.formatted.table", header = T, sep = "\t")
mpileup.dep <- read.table("../data/merged.depths.mpileup.table", header = T, sep = "\t")
nug.dep <- read.table("../data/merged.depths.table", header = T, sep = "\t")
nug.gen.filt <- read.table("../data/merged.filt.genotypes.table", header = T, sep = "\t")
nug.dep.filt <- read.table("../data/merged.filt.depths.table", header = T, sep = "\t")


cleanGATKTable <- function(table) {
  rownames(table) <- paste0(table$CHROM, ".", table$POS)
  table[, c("CHROM", "POS", "REF", "ALT")] <- NULL
  colnames(table) <- gsub("\\..*", "", colnames(table))
  colnames(table) <- gsub("_.*", "", colnames(table))
  colnames(table) <- gsub(".DP", "", colnames(table))
  colnames(table) <- gsub(".GT", "", colnames(table))
  return(as.matrix(table))
}

ref.long <- melt(cleanGATKTable(ref.gen))
mpileup.dep.long <- melt(cleanGATKTable(mpileup.dep))
nug.dep.long <- melt(cleanGATKTable(nug.dep))
nug.gen.filt.long <- melt(cleanGATKTable(nug.gen.filt))
nug.dep.filt.long <- melt(cleanGATKTable(nug.dep.filt))

nug.long <- merge(ref.long, mpileup.dep.long, by = c("Var1", "Var2"))
nug.long <- merge(nug.long, nug.dep.long, by = c("Var1", "Var2"))
nug.long <- merge(nug.long, nug.dep.long, by = c("Var1", "Var2"))
nug.long <- merge(nug.long, nug.gen.filt.long, by=c("Var1", "Var2"))

names(nug.long) <- c("Chrom.pos", "Sample", "Ref", "mpileup.depth",
                     "call.depth", "filt.depth", "Nug")

both.nug <- nug.long %>% 
  filter(Ref != "./.") %>%
  filter(Nug != "./.")

both.nug$ref.a <- sapply(strsplit(both.nug$Ref, split = "/"), "[[", 1)
both.nug$ref.b <- sapply(strsplit(both.nug$Ref, split = "/"), "[[", 2)
both.nug$nug.a <- sapply(strsplit(both.nug$Nug, split = "/"), "[[", 1)
both.nug$nug.b <- sapply(strsplit(both.nug$Nug, split = "/"), "[[", 2)

both.nug$match <- "0"
for (i in 1:nrow(both.nug)) {
  #print(i)
  both.nug[i, "match"] <- length(which(both.nug[i, c("ref.a", "ref.b")] %in% 
                                         both.nug[i, c("nug.a", "nug.b")]))
}

sampleConcordanceNug <- both.nug %>% 
  group_by(Sample) %>% 
  summarize(concordance = mean(as.integer(match)/2*100))

# Remove unidentified sample
sampleStats <- sampleStats[-49, ]

# Range of raw reads
mean(sampleStats$V2) ; sd(sampleStats$V2)
max(sampleStats$V2) ; min(sampleStats$V2)

# What proportion of reads were left after trimgalore?
mean(sampleStats$V3/sampleStats$V2*100)
sd(sampleStats$V3/sampleStats$V2*100)

# What proportion of reads were aligned to EquCab2?
mean(sampleStats$V4/sampleStats$V2*100)
sd(sampleStats$V4/sampleStats$V2*100)

# Number of reads after generating mpileup file
mpileup.dep$Undetermined_S0.DP <- NULL
mpileup.reads <- data.frame(colSums(mpileup.dep[-c(1:4)]))
#names(mpileup.reads) <- c("V5")
sampleStats$V5 <- mpileup.reads[sort(rownames(mpileup.reads)), ]
#mpileup.reads
mean(sampleStats$V5) ; sd(sampleStats$V5)

# Number of reads after using bcftools call 
nug.dep$Undetermined_S0.DP <- NULL
call.reads <- data.frame(colSums(nug.dep[-c(1:4)]))
#names(call.reads) <- c("V6")
sampleStats$V6 <- call.reads[sort(rownames(call.reads)), ]
#call.reads
mean(sampleStats$V6) ; sd(sampleStats$V6)

# Number of reads after filtering variants
nug.dep.filt$Undetermined_S0.DP <- NULL
filt.reads <- data.frame(colSums(nug.dep.filt[-c(1:4)]))
#names(filt.reads) <- c("V7")
sampleStats$V7 <- filt.reads[sort(rownames(filt.reads)), ]
#filt.reads
mean(sampleStats$V7) ; sd(sampleStats$V7)

# What is the proportion of the 279 targets with a genotype? mpileup
nug.gen <- read.table("../data/merged.filt.genotypes.table", header = T, sep = "\t")
nug.gen$Undetermined_S0.GT <- NULL
nug.gen.counts <- data.frame(colSums(nug.gen[-c(1:4)] != "./."))
#names(nug.gen.counts) <- c("V8")
sampleStats$V8 <- nug.gen.counts[sort(rownames(nug.gen.counts)), ]
mean(sampleStats$V8) ; sd(sampleStats$V8)

# What is the percent concordance? Match by sample names
sampleStats$Sample <- gsub("_.*", "", sampleStats$V1)
sampleStats <- merge(sampleStats, sampleConcordanceNug, by="Sample")
head(sampleStats)
# Significant digits for concordance:
sampleStats$concordance <- round(sampleStats$concordance, digits = 1)

# Format table and write to file
names(sampleStats) <- c("Sample", "fullSampleName", "raw", "trimmed", 
                        "aligned", "mpileup.reads", "call.reads",
                        "filt.reads", "num.geno", "concordance")

write.table(sampleStats, "../data/sampleStatisticsFull.csv", sep = ",",
            quote = F, row.names = F)