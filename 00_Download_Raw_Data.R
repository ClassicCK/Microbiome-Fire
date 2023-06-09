#Load Required Packages
library(plyr) #Tidyverse
library(dplyr) #Tidyverse
library(neonUtilities) #Required for neonMicrobe
library(neonMicrobe) #Functions to work with NEON Data
library(ShortRead) #Functions to process Sequence Data
library(Biostrings) #Functions for Sequence Data
library(dada2) #Functions for Sequence Data


# Set Base Directory to Download Data
setBaseDirectory("C:/Users/clk38/iCloudDrive/Documents/Duke/GRSM_Biome")
makeDataDirectories() #Creates file structure to store Metadata/Amplicon Data

#Download Bacterial (16s) and Fungal (ITS) Amplicon Sequence Metadata
meta_16s <- downloadSequenceMetadata(startYrMo = "2015-07", endYrMo = "2022-07",
                                     sites = c("GRSM"),
                                     targetGene = "16S")

meta_ITS <- downloadSequenceMetadata(startYrMo = "2015-07", endYrMo = "2022-07",
                                     sites = c("GRSM"),
                                     targetGene = "ITS")

#Perform Quality Control on Metadata (see neonMicrobe Vignette for Methodology)
meta_16s_qc <- qcMetadata(meta_16S, pairedReads = "Y")
meta_ITS_qc <- qcMetadata(meta_ITS, pairedReads = "Y")

#Collapse into Table
with(meta_16s_qc, table(siteID, sequencerRunID))
with(meta_ITS_qc, table(siteID, sequencerRunID))

#Use Metadata to pull sequence data from online NEON repository
downloadRawSequenceData(meta_16s_qc)
downloadRawSequenceData(meta_ITS_qc)