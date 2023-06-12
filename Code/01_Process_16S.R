# Download and analyze NEON 16s soil microbe marker gene sequence data from GRSM 
library(plyr) #Tidyverse
library(dplyr) #Tidyverse
library(neonUtilities) #Required for neonMicrobe
library(ShortRead) #Required for amplicon analysis
library(Biostrings) #Required for amplicon analysis
library(dada2) #Required for amplicon analysis
library(ggplot2) #For EDA and Figures
library(phangorn) #For generation of Phylogenetic Tree 

# Load neonMicrobe and dependent packages to run 16S amplicon analysis
install.packages("devtools")
devtools::install_github("claraqin/neonMicrobe")

.cran_packages <-  c("ggplot2", "gridExtra")
.bioc_packages <- c("dada2", "msa", "phyloseq")
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

library(neonMicrobe)

# Set working directory
setBaseDirectory("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Duke/GRSM_Biome")

#Re-load metadata to ensure no changes since 00_Download_Raw_Data script
meta_16s <- downloadSequenceMetadata(startYrMo = "2014-07", endYrMo = "2022-07",
                                     sites = c("GRSM"),
                                     targetGene = "16S")

meta_16s_qc <- qcMetadata(meta_16s, pairedReads = "Y")

downloadRawSequenceData(meta_16s_qc)

# Process amplicon sequences: quality control, trimming, alignment

meta <- qcMetadata(meta_16s, pairedReads = "Y", rmFlagged = "Y")
meta$date_ymd <- as.Date(format(as.Date(meta$collectDate, format="%Y-%m-%d %H:%M:%S"), "%Y-%m-%d"))

# Define individual amplicon runs
unique_runs <- unique(meta$sequencerRunID)
fl_nm <- meta$rawDataFileName
head(fl_nm)

# Trim and align amplicon sequences
trim_trackReads <- trimPrimers16S(
  fl_nm, in_subdir = "raw", out_subdir = "1_trimmed",
  meta = meta, 
  multithread = TRUE 
)

r1_files <- file.path(NEONMICROBE_DIR_SEQUENCE(), "16S", grep("_R1", fl_nm, value=TRUE))
plotEEProfile(r1_files, aggregate = TRUE) 

r2_files <- file.path(NEONMICROBE_DIR_SEQUENCE(), "16S", grep("_R2", fl_nm, value=TRUE))
plotEEProfile(r2_files, aggregate = TRUE) 

filter_trackReads <- qualityFilter16S(
  fl_nm, in_subdir = "1_trimmed", out_subdir = "2_filtered",
  meta = meta, truncLen = c(250, 230), maxEE = c(2, 8),
  multithread = TRUE
)

# Infer sample composition with DADA
for(i in 1:length(unique_runs)) {
  meta_thisrun <- meta[which(meta$sequencerRunID==unique_runs[i]),]
  fl_nm_thisrun <- meta_thisrun$rawDataFileName
  dada_out <- runDada16S(
    fl_nm_thisrun, in_subdir = "2_filtered", meta = meta,
    out_seqtab = file.path(NEONMICROBE_DIR_OUTPUTS(), "grsm",
                           paste0("grsm_asv_", unique_runs[i], ".Rds")),
    out_track = file.path(NEONMICROBE_DIR_OUTPUTS(), "grsm",
                          paste0("grsm_track_", unique_runs[i], ".csv")),
    verbose = FALSE,
    multithread = TRUE
  )
}

track <- combineReadTrackingTables16S(trim_trackReads, 
                                      filter_trackReads, 
                                      dada_out$track)

#Write out processed amplicon sequences as checkpoint dataset
write.csv(
  track, 
  file.path(NEONMICROBE_DIR_TRACKREADS(), "16S",
            paste0("track_16s_", sub(" ", "_", gsub(":", "", Sys.time())), ".csv"))
)
rm(trim_trackReads)
rm(filter_trackReads)

seqtab_joined <- mergeSequenceTables(
  tables = file.path(NEONMICROBE_DIR_OUTPUTS(), "grsm",
                     paste0("grsm_asv_", unique_runs, ".Rds"))
)

# Collapses NEON 16S ASV table

seqtab_collapse_filename <- file.path(NEONMICROBE_DIR_OUTPUTS(), "grsm",
                                      "NEON_16S_seqtab_nochim_grsm_COLLAPSED.Rds")
t0 <- Sys.time()
seqtab_collapse <- collapseNoMismatch(seqtab_joined)
saveRDS(seqtab_collapse, seqtab_collapse_filename) #Save Output as checkpoint
t1 <- Sys.time()
t1 - t0 # Check time. 4.37 hrs on Mac w/ M1 chip

# Read in aligned, high quality amplicon table, new checkpoint so that prior
# steps do not need to be run if one wants to alter the parameters around 
# amplicon 
seqtab_collapse <- readRDS(seqtab_collapse_filename)

keep_taxa <- which(colSums(seqtab_collapse) > 10)
keep_samples <- which(rowSums(seqtab_collapse) > 1000)
seqtab_collapse <- seqtab_collapse[,keep_taxa]
seqtab_collapse<- seqtab_collapse[keep_samples,]

# Check size of dataset
print(dim(seqtab_collapse))

# Load Silva to assign taxonomy to ASVs (https://www.arb-silva.de)
silva.url <- "https://zenodo.org/record/3986799/files/silva_nr99_v138_train_set.fa.gz"
download.file(silva.url, NEONMICROBE_DIR_TAXREF(), basename(silva.url))

tax <- assignTaxonomy(
  seqtab_collapse, file.path("/Users/christopher/Library/Mobile Documents/com~apple~CloudDocs/Documents/Duke/GRSM_Biome/data/tax_ref/silva_nr99_v138.1_train_set.fa"), 
  verbose = TRUE, multithread = TRUE
)

saveRDS(tax, file.path(NEONMICROBE_DIR_OUTPUTS(), "grsm.Rds")) #Save Taxonomy Database


# Load in Environmental & Soils Data to add to sequence data for analysis
# Environmental/Soils data collects simultaneously to soil core data that 
# contains information on microbial composition and abundances (ASVs). 
# Aligned with quality controlled metadata
soils <- downloadSoilData(startYrMo = "2014-07", endYrMo = "2022-07",sites = c("GRSM"))
select_columns <- c("domainID", "siteID", "plotID", "plotType", "nlcdClass", "collectDate", 
                    "sampleTiming", "standingWaterDepth", "sampleID", 
                    "horizon", "soilTemp", "litterDepth", "sampleTopDepth", "sampleBottomDepth", 
                    "geneticSampleID", "soilInCaClpH") #environmetnal data of interest
soils <- soils[, select_columns]
sampledata <- merge(meta, soils, all.x=T, by = "geneticSampleID")

remove <- which(duplicated(sampledata$geneticSampleID) == T) #163 Duplicated
sampledata <- sampledata[-remove,]

# Add rownames to sample data table 
rownames(sampledata) <- sampledata$dnaSampleID

seqtab_orig <- seqtab_collapse #define new phyloseq object, prevent overwrite
taxa_orig <- tax #define new phyloseq object, prevent overwrite 

# Remove low-abundance taxa (optional)
keep <- which(colSums(seqtab_orig) > 20)
seqtab <- seqtab_orig[,keep]
print(dim(seqtab))
#> [1]    86 17038

# Remove low-quality samples (optional)
keep <- which(rowSums(seqtab) > 5000)
seqtab <- seqtab[keep,]
print(dim(seqtab))

# Subset to samples present in both dataframes (ASV + Environment)
common_samples <- intersect(rownames(sampledata), rownames(seqtab))
sampledata <- sampledata[common_samples,]
seqtab <- seqtab[common_samples,]

# Check rowname agreement for phyloseq
identical(rownames(seqtab), rownames(sampledata))

# Create phylogenetic trees
seqs <- getSequences(seqtab)
names(seqs) <- seqs # This propagates to the tip labels of the tree
mult <- msa(seqs, method="ClustalW", type="dna", order="input")

# Align sequences, trimmed to approx. 230 bp for tree node generation
phang.align <- as.phyDat(mult, type="DNA", names=getSequence(seqtab))
dm <- dist.ml(phang.align, verbose = TRUE)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

# negative edges length changed to 0!
# run multiple gaussian bootstraps to generate final phylogenetic tree
# Warning! Process takes approx. 80 hrs on Mac M1 chip
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)

# Combine into phyloseq object
ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE),
               sample_data(sampledata),
               tax_table(taxa_orig),
               phy_tree(fitGTR$tree))

# store the DNA sequences of our ASVs in the refseq slot of the phyloseq object,
# and then rename our taxa to a short string
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

saveRDS(ps, file.path(NEONMICROBE_DIR_OUTPUTS(), "NEON_16S_phyloseq_grsm.Rds"))

# This is now the working Phyloseq object for all analyses!