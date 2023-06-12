# Load in required pacakges
library(magrittr) # Pipe operator package
library(neonMicrobe) # For work with NEON Aplicon sequence data
library(microViz) # RDA, PCA, PCoA, and other oridnations for sequence data
library(mattsUtils) # Subset and Filter phyloseq object
library(BiocManager) # Install Bioconductor Packages
library(phyloseq) # Manage amplicon sequence data in structured environment
library(plyr) # Tidyverse
library(dplyr) # Tidyverse
library(tidyr) # Tidyverse
library(ggplot2) # Data visualization
library(cowplot) # Data visualization
library(ggthemes) # Data visualization
library(RColorBrewer) # Data visualization
library(ggh4x) # Data visualization
library(pals)# Data visualization
library(indicspecies) # Indicator Species Analysis Package
library(ggtree)# Data visualization
library(xtable)# Save table as LateX

# Load in previosuly generated phyloseq object (01_Process_16S.R)
list.files(.libPaths()[1])  %>%
  purrr::walk(~library(.x, character.only = TRUE))

ps <- readRDS(file.path(NEONMICROBE_DIR_OUTPUTS(), "NEON_16S_phyloseq_grsm.Rds"))
ps <- taxa_prune(ps, "Archaea", classification="Kingdom") #Remove Archaea
ps <- subset_taxa(ps, is.na(Phylum) == F) # Remove ASVs with no Silva Phlya
ps <- tax_fix(ps) # quality control check on phyloseq object to ensure all data (enviornment, phylogeny, taxonomy, ASVs) aligns
ps <- phyloseq_validate(ps, remove_undetected = TRUE) # Removes remaining unaligned ASVs

# Subset data to relevant environmental variables, and assign burn status 
# to the environmental data. Burn status determined from reported burn
# information provided by the National Park Service and aligned to the NEON
# plots in Figure 1 (generated in 'mapping.r' script)
keep <- c("plotID.x", "date_ymd", "nlcdClass", "dnaSampleID",
          "sampleTiming", "sampleID", "horizon", "soilTemp", "soilInCaClpH")
keep.sites <- c("GRSM_002", "GRSM_001","GRSM_007", "GRSM_016", "GRSM_055", 
                "GRSM_058", "GRSM_059", "GRSM_060", "GRSM_003") #Burned/Unburned Sites

ps <- subset_samples(ps, ps@sam_data$plotID.x %in% keep.sites)

ps@sam_data <- ps@sam_data[,keep]

# Assign Burn Status
for(i in 1:nrow(ps@sam_data)){
  if(ps@sam_data$plotID.x[i] == "GRSM_001"){
    ps@sam_data$fire[i] <- "Unburned"
  }
  if(ps@sam_data$plotID.x[i] == "GRSM_002"){
    ps@sam_data$fire[i] <- "Unburned"
  }
  if(ps@sam_data$plotID.x[i] == "GRSM_003"){
    ps@sam_data$fire[i] <- "Burned"
  }
  if(ps@sam_data$plotID.x[i] == "GRSM_007"){
    ps@sam_data$fire[i] <- "Burned"
  }
  if(ps@sam_data$plotID.x[i] == "GRSM_016"){
    ps@sam_data$fire[i] <- "Unburned"
  }
  if(ps@sam_data$plotID.x[i] == "GRSM_055"){
    ps@sam_data$fire[i] <- "Burned"
  }
  if(ps@sam_data$plotID.x[i] == "GRSM_058"){
    ps@sam_data$fire[i] <- "Burned"
  }
  if(ps@sam_data$plotID.x[i] == "GRSM_059"){
    ps@sam_data$fire[i] <- "Burned"
  }
  if(ps@sam_data$plotID.x[i]== "GRSM_060"){
    ps@sam_data$fire[i] <- "Burned"
  }
}

# Generate a 'year' column for downstream plotting and analysis, and set factor
# level as Pre- (2016) and Post- (2017) Chimney Tops 2 Fire
ps@sam_data$year <- stringr::str_extract(ps@sam_data$date_ymd, "^.{4}")

for(i in 1:nrow(ps@sam_data)){
  if(ps@sam_data$year[i] == "2016"){
    ps@sam_data$Time[i] <- "Pre-"
  }
  if(ps@sam_data$year[i] == "2017"){
    ps@sam_data$Time[i] <- "Post-"
  }
}

# Generate Alpha and Beta Diversity estimates from phyloseq object
ps_rare <- rarefy_even_depth(ps, sample.size=10000, rngseed=101010) #Value from Rocca et. al. (2022) Paper
alpha_div <- cbind("dnaSampleID" = get_variable(ps_rare, "dnaSampleID"), # Generate Alpha- and Shannon Diversity Estimates
                   estimate_richness(ps_rare, measures = c("Observed", "Shannon")))
alpha_div <- merge(sampledata, alpha_div, by="dnaSampleID", all.y=TRUE) #Merge with Environmental/Soils Data 
alpha_div$fire <- NA #Set up new column to generate fire burn status

#Add in Fire Data
for(i in 1:nrow(alpha_div)){
  if(alpha_div$plotID.x[i] == "GRSM_001"){
    alpha_div$fire[i] <- "Unburned"
  }
  if(alpha_div$plotID.x[i] == "GRSM_002"){
    alpha_div$fire[i] <- "Unburned"
  }
  if(alpha_div$plotID.x[i] == "GRSM_003"){
    alpha_div$fire[i] <- "Burned"
  }
  if(alpha_div$plotID.x[i] == "GRSM_007"){
    alpha_div$fire[i] <- "Burned"
  }
  if(alpha_div$plotID.x[i] == "GRSM_016"){
    alpha_div$fire[i] <- "Unburned"
  }
  if(alpha_div$plotID.x[i] == "GRSM_055"){
    alpha_div$fire[i] <- "Burned"
  }
  if(alpha_div$plotID.x[i] == "GRSM_058"){
    alpha_div$fire[i] <- "Burned"
  }
  if(alpha_div$plotID.x[i] == "GRSM_059"){
    alpha_div$fire[i] <- "Burned"
  }
  if(alpha_div$plotID.x[i]== "GRSM_060"){
    alpha_div$fire[i] <- "Burned"
  }
}

# set up generation of year and Pre-/Post- Burn variables
alpha_div$year <- stringr::str_extract(alpha_div$date_ymd, "^.{4}")
alpha_div$yearmonth <- stringr::str_extract(alpha_div$date_ymd, "^.{7}")

alpha_div$Time <- NA
for(i in 1:nrow(alpha_div)){
  if(alpha_div$year[i] == "2016"){
    alpha_div$Time[i] <- "Pre-"
  }
  if(alpha_div$year[i] == "2017"){
    alpha_div$Time[i] <- "Post-"
  }
}

alpha_div$fire <- factor(alpha_div$fire, levels = c("Unburned", "Burned"))
alpha_div$Time <- factor(alpha_div$Time, levels = c("Pre-", "Post-"))

# Statistical Tests
alpha.burned <- alpha_div[alpha_div$fire == "Burned",] #Burned Plots
alpha.unburned <- alpha_div[alpha_div$fire == "Unburned",] #Unburned Plots

# Wilcox Test between Pre- and Post- Burn in Burned Plots
pre <- alpha.burned[alpha.burned$Time == "Pre-",c("Observed", "Shannon")]
post <- alpha.burned[alpha.burned$Time == "Post-",c("Observed", "Shannon")]

burned.obs <- wilcox.test(x = pre$Observed, y = post$Observed,
                          alternative = c("two.sided"),
                          mu = 0, paired = FALSE, exact = FALSE, correct = TRUE,
                          conf.int = FALSE, conf.level = 0.95)

burned.shannon <- wilcox.test(x = pre$Shannon, y = post$Shannon,
                              alternative = c("two.sided"),
                              mu = 0, paired = FALSE, exact = NULL, correct = TRUE,
                              conf.int = FALSE, conf.level = 0.95)
print(burned.obs)
summary(burned.obs)

print(burned.shannon)
summary(burned.shannon)

# Wilcox Test between Pre- and Post- Burn in Burned Plots
pre <- alpha.unburned[alpha.unburned$Time == "Pre-",c("Observed", "Shannon")]
post <- alpha.unburned[alpha.unburned$Time == "Post-",c("Observed", "Shannon")]

unburned.obs <- wilcox.test(x = pre$Observed, y = post$Observed,
                            alternative = c("two.sided"),
                            mu = 0, paired = FALSE, exact = NULL, correct = TRUE,
                            conf.int = FALSE, conf.level = 0.95)

unburned.shannon <- wilcox.test(x = pre$Shannon, y = post$Shannon,
                                alternative = c("two.sided"),
                                mu = 0, paired = FALSE, exact = NULL, correct = TRUE,
                                conf.int = FALSE, conf.level = 0.95)

print(unburned.obs)
summary(unburned.obs)

print(unburned.shannon)
summary(unburned.shannon)

# Wilcox Test between Burned and Unburned Plots Pre- Burn (make sure there isn't already a significant difference before fire)
un <- alpha.unburned[alpha.unburned$Time == "Pre-",c("Observed", "Shannon")] #Unburned plots in 2016
b <- alpha.burned[alpha.burned$Time == "Pre-",c("Observed", "Shannon")] #Burned plots in 2016

pre.test.obs <- wilcox.test(x = un$Observed, y = b$Observed,
                            alternative = c("two.sided"),
                            mu = 0, paired = FALSE, exact = FALSE, correct = TRUE,
                            conf.int = FALSE, conf.level = 0.95)

pre.test.shannon <- wilcox.test(x = un$Shannon, y = b$Shannon,
                                alternative = c("two.sided"),
                                mu = 0, paired = FALSE, exact = FALSE, correct = TRUE,
                                conf.int = FALSE, conf.level = 0.95)

print(pre.test.obs)
summary(pre.test.obs)

print(pre.test.shannon )
summary(pre.test.shannon)

# PERMANOVA???????


# Figure 2. Boxplot of Pre- and Post- Burn in Burned and Unburned plots for 
# total (A) ASV richness (alpha diversity) and (B) Shannon Diversity. Panels 
# (C) and (D) show recovery of alpha and shannon diversity over time in 
# burned plots only.

Fig2a <- ggplot(alpha_div, aes(x=fire, y=Observed, fill = Time)) +
  geom_boxplot(na.rm = T) +
  ylab("Observed ASV Richness") +
  xlab("Fire Status") +
  scale_fill_manual(values = c("#002244","#FB4F14")) +
  theme_cowplot(12) +
  labs(color = expression("Time")) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position="none") +
  stat_compare_means(mapping = NULL, comparisons = NULL, hide.ns = TRUE,
                     label = "p.signif",  label.x = NULL, label.y = NULL)

Fig2b <- ggplot(alpha_div, aes(x=fire, y=Shannon, fill = Time)) +
  geom_boxplot(na.rm = T) +
  ylab("Observed Shannon Diversity") +
  xlab("Fire Status") +
  scale_fill_manual(values = c("#002244","#FB4F14")) +
  theme_cowplot(12) +
  labs(color = expression("Time")) +
  ylim(5, 7)  +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position="none") +
  stat_compare_means(mapping = NULL, comparisons = NULL, hide.ns = TRUE,
                     label = "p.signif",  label.x = NULL, label.y = NULL)

# Pull out just Burned and Post-Burn to generate time/recovery component
alpha_burned <- subset(alpha_div, alpha_div$fire == "Burned")
alpha_burned <- subset(alpha_burned, alpha_burned$Time == "Post-")
alpha_burned$sampleTiming <- factor(alpha_burned$sampleTiming, levels = c("winterSpringTransition", "peakGreenness",  "fallWinterTransition"))

#Sort into pre and post, combine all pre into one observation. 
#Leave post as seperate observations that are +3, +6, +9 months post burn
for(i in 1:nrow(alpha_div)){
  if(alpha_div$year[i] == 2016){
    alpha_div$timepoint[i] <- "2016"  
  }
  if(alpha_div$year[i] == 2017){
    if(alpha_div$sampleTiming[i] == "winterSpringTransition"){
      alpha_div$timepoint[i] <- "+3 Mon."
    }
    if(alpha_div$sampleTiming[i] == "peakGreenness"){
      alpha_div$timepoint[i] <- "+6 Mon."
    }
    if(alpha_div$sampleTiming[i] == "fallWinterTransition"){
      alpha_div$timepoint[i] <- "+9 Mon."
    }
  }
}

alpha_div$timepoint <- factor(alpha_div$timepoint, levels = c("2016", 
                                                              "+3 Mon.",
                                                              "+6 Mon.",
                                                              "+9 Mon."))

#Generate time recovery plots for alpha and shannon diversity
Fig2C <- ggplot(alpha_div, aes(x=timepoint, y=Observed, fill=timepoint)) +
  geom_boxplot(na.rm = T) +
  ylab("Observed ASV Richness") +
  xlab("Sample Timing") + 
  theme_cowplot(12) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12), 
        legend.position="none") +
  scale_fill_manual(values = c("#002244","#FB4F14", "#FB4F14", "#FB4F14"))

Fig2D <- ggplot(alpha_div, aes(x=timepoint, y=Shannon, fill=timepoint)) +
  geom_boxplot(na.rm = T) +
  ylab("Observed Shannon Diversity") +
  xlab("Sample Timing") + 
  theme_cowplot(12) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12), 
        legend.position="none") +
  scale_fill_manual(values = c("#002244","#FB4F14", "#FB4F14", "#FB4F14"))

legend <- get_legend(
  Fig2a + guides(color = guide_legend(nrow = 1)) + theme(legend.position = "bottom", 
                                                         legend.title = element_text(size = 12), 
                                                         legend.text = element_text(size = 12)))
Figure_Tmp <- plot_grid(Fig2a,Fig2b,Fig2C,Fig2D, labels = c("A", "B", "C", "D"), nrow = 2)
Figure2 <- plot_grid(Figure_Tmp, legend, ncol = 1, rel_heights  = c(1, .1), align = "v", axis = "lr") 

ggsave(filename= paste0(getwd(), "/Figures/Figure2.pdf"), 
       plot = Figure2, width = 180, height = 180, units=c("mm"), dpi=600) 

# Figure 3. 

# Define shapes for plot points
good.shapes = c(1:25,33:127)

# Figure 3A
new.ps <- ps %>%
  tax_glom(taxrank = "Genus")%>% #Filter to Genus level for plotting, reduce ASV to lowest Silva taxonomic unit 
  psmelt()%>%
  arrange(Genus)

mycolors<-as.vector(brewer.paired(32)) #Set Color Scheme

new.ps$fire <- factor(new.ps$fire, levels = c("Unburned", "Burned"))
new.ps$Time <- factor(new.ps$Time, levels = c("Pre-", "Post-"))

Fig3A <- ggplot(new.ps,aes(x=Time, y=Abundance, fill=Phylum))+
  geom_bar(stat = "identity", position = "stack") + 
  facet_nested(.~fire+plotID.x,scales="free") +
  theme_cowplot(12) +
  scale_fill_manual(values=(mycolors), name="Phylum")+
  xlab("Samples") + ylab("Relative Abundance") +
  theme(axis.text=element_text(size=8),
        axis.title=element_blank(),
        axis.text.x =element_text(size=8),
        axis.ticks = element_blank(),
        legend.text=element_text(size=8),
        legend.position = "bottom", 
        strip.background = element_rect(fill="white", colour = NA),
        strip.text = element_text(color="black",size=8, face="bold")) +
  guides(fill=guide_legend(keywidth = 0.5,keyheight = 0.5, ncol = 6))

ggsave(filename= paste0(getwd(), "/Figures/Figure3A.pdf"), 
       plot = Fig3A, width = 180, height = 180, units=c("mm"), dpi=600)

# Perform PCoA and RDA analysis to generate Figure 3B and 3C. Print tables
# for supplementary information

# Figure 3B. MDA plot of Bray-Curtis Distance
# Color by Pre/Post Burn, Shape by Plot ID
ps.bray <-  ps %>%
  tax_filter(min_prevalence = 0.1, tax_level = "Genus") %>% #Filter for at least 10 observations per Genus (Rocca 2022)
  dist_calc(dist = "bray") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(color = "Time", shape = "plotID.x", size = 2, alpha = 1,
           show.legend = FALSE) +
  scale_colour_manual(values = c("#FB4F14","#002244", "")) +
  scale_shape_manual(values=good.shapes)

# Generate Secondary repeat figure to pull legend from
ps.bray.legend <-  ps %>%
  tax_filter(min_prevalence = 0.1, tax_level = "Genus") %>% #Filter for at least 10 observations per Genus (Rocca 2022)
  dist_calc(dist = "bray") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(color = "Time", shape = "plotID.x", size = 2, alpha = 1,
           show.legend = TRUE) +
  scale_colour_manual(values = c("#FB4F14","#002244", "")) +
  scale_shape_manual(values=good.shapes)


Fig3B <- ps.bray +
  theme(axis.text=element_text(size=8),
        axis.title=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        legend.text=element_text(size=8)) +
  theme_cowplot(12) +
  labs(color='Fire Status', shape = "Plot") + 
  ggside::geom_xsidedensity(aes(fill = Time), alpha = 0.5, show.legend = FALSE) +
  ggside::geom_ysidedensity(aes(fill = Time), alpha = 0.5, show.legend = FALSE) +
  ggside::theme_ggside_void()

ggsave(filename= paste0(getwd(), "/Figures/Figure3B.pdf"), 
       plot = Fig3B, width = 180, height = 180, units=c("mm"), dpi=600)

# Get PCoA MDS Figure
ps.pcoa <- ps %>%
  tax_filter(min_prevalence = 0.1, tax_level = "Genus") %>%
  dist_calc(dist = "bray") %>%
  ord_calc(method = "PCoA") %>%
  ord_get() %>%
  phyloseq::plot_scree() +
  theme(axis.text=element_text(size=8),
        axis.title=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        legend.text=element_text(size=8)) +
  theme_cowplot(12) +
  xlab("Axis") +
  ylab("Eigenvalue")

ggsave(filename= paste0(getwd(), "/Figures/FigureS2.pdf"), 
       plot = ps.pcoa, width = 180, height = 90, units=c("mm"), dpi=600)


# Figure 3C. RDA Plot (Redundancy Analysis)
ps.weighted <- ps %>% #Weight sites by scaling and centering
  ps_mutate( #Redefine Factor Levels for Eigenvector plotting
    pH = c(scale(soilInCaClpH, center = TRUE, scale = TRUE)),
    Fire = if_else(fire == "Burned" & Time == "Post-", 1, 0),
    Spring = if_else(sampleTiming == "winterSpringTransition" & Time == "Post-", 1, 0),
    Summer = if_else(sampleTiming == "peakGreenness" & Time == "Post-", 1, 0),
    Fall = if_else(sampleTiming == "fallWinterTransition" & Time == "Post-", 1, 0),
    GRSM_001 = if_else(plotID.x == "GRSM_001", 1, 0),
    GRSM_002 = if_else(plotID.x == "GRSM_002", 1, 0),
    GRSM_003 = if_else(plotID.x == "GRSM_003", 1, 0),
    GRSM_007 = if_else(plotID.x == "GRSM_007", 1, 0),
    GRSM_016 = if_else(plotID.x == "GRSM_016", 1, 0),
    GRSM_055 = if_else(plotID.x == "GRSM_055", 1, 0),
    GRSM_058 = if_else(plotID.x == "GRSM_058", 1, 0),
    GRSM_059 = if_else(plotID.x == "GRSM_059", 1, 0),
    GRSM_060 = if_else(plotID.x == "GRSM_060", 1, 0))

constr.ord <- ps.weighted %>% #Run constrained ordination, holding site ID as a condition
  tax_filter(min_prevalence = 0.1, tax_level = "Genus") %>% #Filter for at least 10 observations per Genus (Rocca 2022)
  tax_agg("Genus") %>%
  tax_transform("hellinger") %>%
  ord_calc(method = "RDA", conditions = c("GRSM_001", "GRSM_002", "GRSM_003", 
                                          "GRSM_007", "GRSM_055", "GRSM_058", 
                                          "GRSM_059","GRSM_060"), 
           constraints = c("pH", "Fire", "Spring", "Summer", "Fall")) %>% 
  ord_plot(color = "Time", shape = "plotID.x", size = 2, alpha = 1, show.legend = FALSE) +
  scale_colour_manual(values = c("#FB4F14","#002244", "")) +
  scale_shape_manual(values=good.shapes)

# Plot RDA
Fig3C <- constr.ord +
  theme(axis.text=element_text(size=8),
        axis.title=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        legend.text=element_text(size=8)) +
  theme_cowplot(12) +
  labs(color='Time', shape = "Plot") + 
  ggside::geom_xsidedensity(aes(fill = Time), alpha = 0.5, show.legend = FALSE) +
  ggside::geom_ysidedensity(aes(fill = Time), alpha = 0.5, show.legend = FALSE) +
  ggside::theme_ggside_void()

ggsave(filename= paste0(getwd(), "/Figures/Figure3A.pdf"), 
       plot = Fig3C, width = 180, height = 180, units=c("mm"), dpi=600)

# Get RDA Table
ps.rda <- ps.weighted %>% #Run constrained ordination, holding site ID as a condition
  tax_filter(min_prevalence = 0.1, tax_level = "Genus") %>% #Filter for at least 10 observations per Genus (Rocca 2022)
  tax_agg("Genus") %>%
  tax_transform("hellinger") %>%
  ord_calc(method = "RDA", conditions = c("GRSM_001", "GRSM_002", "GRSM_003", 
                                          "GRSM_007", "GRSM_055", "GRSM_058", 
                                          "GRSM_059","GRSM_060"), 
           constraints = c("pH", "Fire", "Spring", "Summer", "Fall"))

print(ord_get(ps.rda))

# Put together Figure 3
legend <- get_legend(ps.bray.legend + guides(color = guide_legend(nrow = 1)) + 
                       theme(legend.position = "bottom", legend.title = element_text(size = 8), legend.text = element_text(size = 8)))
Figure.Tmp <- plot_grid(Fig3B, Fig3C, nrow = 1, align = "v", axis = "tb", labels = c("B", "C"))
Figure.Tmp.2 <- plot_grid(Figure.Tmp, legend, ncol= 1, align = "h", rel_heights  = c(1,.1), axis = "lr")
Figure3 <- plot_grid(Fig3A, Figure.Tmp.2, ncol = 1, align = "v", axis = "tb", labels = c("A", "")) 

ggsave(filename= paste0(getwd(),"/Figures/Figure3.pdf"), 
       plot = Figure3, width = 180, height = 240, units=c("mm"), dpi=600)

# Figure 4 (Phylogenetic Tree Color Coded by Indicator Species Analysis)
# If the site classification vector (environmental data) is obtained independently of species data, 
# the significance of statistical tests carried out on the indicator species will be meaningful. 
# For example, one could classify the sites using environmental data before indicator species analysis. 
# An example is found in Borcard et al. 2011

# Specific Taxa Level Indicator Species Analysis
ps.taxa <- tax_glom(ps, taxrank = "Genus") #Choose Taxa Level for Indicator "Species". Here, we are analyzing at the genus level

# Need to separate into matrices the environmental data and the species data, 
# and for our indicator species analysis, we are interested in post-burn 
# plots that had fire. Thus, comparing species that are indicators of whether 
# a burn did occur or did not occur 

ps.matrix <- veganotu(ps.taxa) #Species Matrix
ps.env <- veganenv(ps.taxa) #Environmental Matrix 

ps.ven.full <- ps.env[,c("fire", "Time")] #Both time points for env matrix
reduced <- which(ps.ven.full[,"fire"] == "Burned") #Only Burned Plots
ps.ven.reduced <- ps.ven.full[reduced,] #Burned Plots, both time points
ps.matrix.reduced <- ps.matrix[reduced,] #Burned Plots, both time points
ps.ven.reduced <- as.data.frame(ps.ven.reduced) 
ps.ven.reduced$Time <- as.factor(ps.ven.reduced$Time)

ps.ven.reduced <- ps.ven.reduced[,-1] #Eliminate column that is only rowID

#Indic Species Analysis for Genus Level
set.seed(811)
indval <- multipatt(ps.matrix.reduced, ps.ven.reduced, 
                    fun = "r.g", control = how(nperm=811)) #set seed and permutations

print(summary(indval))

#Assign Indicator Species to proper factor (indicative of burn or no burn)
#Response to Fire
positive.response <- c('ASV589',
                       'ASV1920',
                       'ASV382',
                       'ASV2222',
                       'ASV86',
                       'ASV405',
                       'ASV5508',
                       'ASV2332',
                       'ASV2461',
                       'ASV1439',
                       'ASV4056',
                       'ASV2775',
                       'ASV2064',
                       'ASV2108',
                       'ASV120',
                       'ASV2331',
                       'ASV2021',
                       'ASV673',
                       'ASV5312',
                       'ASV3257',
                       'ASV1060',
                       'ASV975',
                       'ASV1224',
                       'ASV1327',
                       'ASV4271',
                       'ASV749',
                       'ASV1190',
                       'ASV330',
                       'ASV341',
                       'ASV903',
                       'ASV1577',
                       'ASV1182',
                       'ASV508',
                       'ASV514',
                       'ASV736',
                       'ASV694',
                       'ASV860',
                       'ASV1582',
                       'ASV1742',
                       'ASV3755',
                       'ASV1400',
                       'ASV953')

negative.response <- c('ASV813',
                       'ASV721',
                       'ASV155',
                       'ASV490',
                       'ASV25',
                       'ASV5',
                       'ASV1017',
                       'ASV15',
                       'ASV19',
                       'ASV365',
                       'ASV1398',
                       'ASV4542',
                       'ASV536',
                       'ASV112',
                       'ASV498',
                       'ASV181',
                       'ASV388',
                       'ASV3023',
                       'ASV199',
                       'ASV1276',
                       'ASV504',
                       'ASV298',
                       'ASV804',
                       'ASV3489',
                       'ASV242',
                       'ASV3843',
                       'ASV1442',
                       'ASV379',
                       'ASV1828',
                       'ASV153',
                       'ASV1884',
                       'ASV83',
                       'ASV2055',
                       'ASV569',
                       'ASV99',
                       'ASV4282',
                       'ASV1598',
                       'ASV550',
                       'ASV1889',
                       'ASV340',
                       'ASV202',
                       'ASV769',
                       'ASV57',
                       'ASV12',
                       'ASV4552',
                       'ASV3730',
                       'ASV547',
                       'ASV5388',
                       'ASV594',
                       'ASV251',
                       'ASV4548',
                       'ASV1232',
                       'ASV2235',
                       'ASV76',
                       'ASV1057',
                       'ASV425',
                       'ASV2',
                       'ASV371',
                       'ASV901',
                       'ASV2232',
                       'ASV968',
                       'ASV3686',
                       'ASV1326',
                       'ASV4094',
                       'ASV1306',
                       'ASV39',
                       'ASV4372',
                       'ASV17',
                       'ASV4665',
                       'ASV374',
                       'ASV191',
                       'ASV4',
                       'ASV1384',
                       'ASV132',
                       'ASV3434',
                       'ASV429',
                       'ASV1288',
                       'ASV1142',
                       'ASV629',
                       'ASV5559',
                       'ASV2608',
                       'ASV926',
                       'ASV907',
                       'ASV2740',
                       'ASV1971')

positive.response <- as.data.frame(positive.response)
negative.response <- as.data.frame(negative.response)

# rename factor levels and species IDs
negative.response$response <- "Burn (-)"
positive.response$response <- "Burn (+)"

colnames(negative.response)[1] <- "ASV"
colnames(positive.response)[1] <- "ASV"

responses <- rbind(positive.response, negative.response)

#Bind Indicator Species Analysis Data to Phyloseq Object
ps.genus <- ps %>%
  tax_glom(taxrank = "Genus") #Aggregate to same level as Indicator Species Analysis

tree.data <- tax_table(ps.genus)
rownames(responses) <- responses$ASV

tree.data <- merge(tree.data, responses, by = 0, all.x=TRUE) #Merge Indicator Responses
colnames(tree.data)[1] <- "ASV"
tree.data <- tree.data[,-8]
print(colnames(tree.data))
rownames(tree.data) <- tree.data$ASV 

#Assign neutral responses to genera that are non-indicative/significant
for(i in 1:nrow(tree.data)){ 
  if(is.na(tree.data$response[i] == TRUE)){
    tree.data$response[i] <- "Burn (0)"
  }
  else{}
}

# Generate Phylogenetic Tree Color Coded by Clade
tree.data$SuperPhylum <- "NA" # Add Clade (Super Phylum) to Tree Matrix
tr <- ps.taxa@phy_tree # Generate Seperate tree for linking in Figure 4

# Add in Clade (Super Phylum) groupings to trees. Taken from MOL definitions (https://earthmicrobiome.org)
for(i in 1:nrow(tree.data)){
  if(tree.data$Phylum[i] %in% c("Proteobacteria", "Myxococcota", "Acidobacteriota", "Bdellovibrionota", "RCP2-54", "Desulfobacterota", "Methylomirabilota", "NB1-j", 
                                "SAR324 clade(Marine group B)", "Dependentiae",
                                "Nitrospirota", "Entotheonellaeota", "MBNT15")){
    tree.data$SuperPhylum[i] <- "PANNAM"
  }
  if(tree.data$Phylum[i] %in% c("Actinobacteriota", "Firmicutes", "Armatimonadota", 
                                "WPS-2")){
    tree.data$SuperPhylum[i] <- "Terrabacteria"
  }
  if(tree.data$Phylum[i] %in% c("Planctomycetota", "Verrucomicrobiota")){
    tree.data$SuperPhylum[i] <- "PVC"
  }
  if(tree.data$Phylum[i] %in% c("Bacteroidota", "Latescibacterota", "FCPU426",
                                "Gemmatimonadota", "Fibrobacterota")){
    tree.data$SuperPhylum[i] <- "FCB"
  }
  if(tree.data$Phylum[i] %in% c("Patescibacteria", "GAL15", "Bacteria Kingdom", "FCPU426")){
    tree.data$SuperPhylum[i] <- "Bacteria Candidate Group"
  }
  if(tree.data$Phylum[i] %in% c("Chloroflexi")){
    tree.data$SuperPhylum[i] <- "CCD"
  }
  if(tree.data$Phylum[i] %in% c("Cyanobacteria")){
    tree.data$SuperPhylum[i] <- "CMS"
  }
  if(tree.data$Phylum[i] %in% c("Elusimicrobiota", "Spirochaetota")){
    tree.data$SuperPhylum[i] <- "Gracilicutes"
  }
  if(tree.data$Phylum[i] %in% c("Deinococcota")){
    tree.data$SuperPhylum[i] <- "DST"
  }
}

#Add SuperPhylum Column to phyloseq object
ps.taxa <- tax_names2rank(ps.taxa, colname = "SuperPhylum")

for(i in 1:nrow(ps.taxa@tax_table@.Data)){
  if(ps.taxa@tax_table@.Data[i,2] %in% c("Proteobacteria", "Myxococcota", "Acidobacteriota", "Bdellovibrionota", "RCP2-54", "Desulfobacterota", "Methylomirabilota", "NB1-j", 
                                           "SAR324 clade(Marine group B)", "Dependentiae",
                                           "Nitrospirota", "Entotheonellaeota", "MBNT15")){
    ps.taxa@tax_table@.Data[i,7] <- "PANNAM"
  }
  if(ps.taxa@tax_table@.Data[i,2] %in% c("Actinobacteriota", "Firmicutes", "Armatimonadota", 
                                           "WPS-2")){
    ps.taxa@tax_table@.Data[i,7] <- "Terrabacteria"
  }
  if(ps.taxa@tax_table@.Data[i,2] %in% c("Planctomycetota", "Verrucomicrobiota")){
    ps.taxa@tax_table@.Data[i,7] <- "PVC"
  }
  if(ps.taxa@tax_table@.Data[i,2] %in% c("Bacteroidota", "Latescibacterota", "FCPU426",
                                           "Gemmatimonadota", "Fibrobacterota")){
    ps.taxa@tax_table@.Data[i,7] <- "FCB"
  }
  if(ps.taxa@tax_table@.Data[i,2] %in% c("Patescibacteria", "GAL15", "Bacteria Kingdom", "FCPU426")){
    ps.taxa@tax_table@.Data[i,7] <- "Bacteria Candidate Group"
  }
  if(ps.taxa@tax_table@.Data[i,2] %in% c("Chloroflexi")){
    ps.taxa@tax_table@.Data[i,7] <- "CCD"
  }
  if(ps.taxa@tax_table@.Data[i,2] %in% c("Cyanobacteria")){
    ps.taxa@tax_table@.Data[i,7] <- "CMS"
  }
  if(ps.taxa@tax_table@.Data[i,2] %in% c("Elusimicrobiota", "Spirochaetota")){
    ps.taxa@tax_table@.Data[i,7] <- "Gracilicutes"
  }
  if(ps.taxa@tax_table@.Data[i,2] %in% c("Deinococcota")){
    ps.taxa@tax_table@.Data[i,7] <- "DST"
  }
}

tr.phylum = data.frame(id=tr$tip.label, val = rep(1, 550), group = tree.data$SuperPhylum)
tr.response = data.frame(id=tr$tip.label, val = rep(1, 550), group = tree.data$response)

# Phylogenetic Tree Nodes down to Genus level
df <- as.data.frame(tr.response[,3])
rownames(df) <- tr$tip.label

#Generate Phylogenetic Tree Figure (Figure 4) with nodes colored by Clade, 
# nodes are to genus level. Add in information on indicator species
p <- ggtree(ps.taxa, branch.length='none', layout="circular", open.angle=10) +
  geom_tippoint(aes(color=SuperPhylum), shape = 15, show.legend = TRUE, size = 1)

Fig4 <- gheatmap(p, df, offset=.8, width=.2,
               colnames_angle=95, colnames_offset_y = .25)  +
  scale_fill_manual(values=c("#a97851", "#7851A9", "white")) +
  theme(axis.text=element_blank(),
        axis.title=element_blank(),
        axis.text.x =element_blank(),
        axis.ticks = element_blank(),
        legend.text=element_text(size=8),
        legend.position = "none")

legend <- get_legend(p + guides(color = guide_legend(ncol = 6)) + 
                       theme(legend.position = "bottom", legend.title = element_text(size = 8), legend.text = element_text(size = 8)))

legend.2 <- get_legend(Fig4 + guides(color = guide_legend(nrow = 1)) + 
                         theme(legend.position = "bottom", legend.title = element_text(size = 8), legend.text = element_text(size = 8)))
Legends.Tmp <- plot_grid(legend, legend.2, ncol = 1, align = "h", axis = "lr")
Figure.Tmp <- plot_grid(Fig4, legend.2, ncol= 1, align = "h", rel_heights  = c(1,.1), axis = "lr")

ggsave(filename= paste0(getwd(),"/Figures/Figure4.pdf"), 
       plot = Figure.Tmp, width = 180, height = 240, units=c("mm"), dpi=600) 


# Change in Abundance Tables
otu.table <- as.matrix(otu_table(ps.taxa))
abundance.table <- otu.table@.Data
sum.abund <- colSums(abundance.table)
tr.abundance = data.frame(id=tr$tip.label, val = sum.abund, group = tree.data$Phylum)

# Chi-Squared Contingency Table on ASV repsonses by Clade/Phylum
# Create Table of Indicator Species by Clade
response <- merge(tr.phylum, tr.response, by = "id")
response <- response[,c("group.x", "group.y")]
colnames(response) <- c("Clade", "Indicator")
summary.response <- with(response,table(Clade,Indicator)) #Generate Table

print(summary.response)

#Generate a chi-squared contingency table on response by Clade
x.table <- data.frame(
  Clade = c("Bacteria Candidate Group", "CCD", "CMS", "DST", "FCB", "Gracilicutes", "PANNAM", "PVC", "Terrabacteria"),
  `Burn (-)` = c(0, 6, 1, 0, 7, 1, 31, 20, 19),
  `Burn (0)` = c(25, 26, 2, 0, 54, 3, 205, 34, 74),
  `Burn (+)` = c(1, 0, 1, 1, 5, 0, 23, 2, 9)
)

x.table <- x.table[,-3] #Remove no response (neutral) column
colnames(x.table)[2] <- "Burn (-)" #Add indicator names
colnames(x.table)[3] <- "Burn (+)"
rownames(x.table) <- x.table$Clade #Add rownames (clades)
chi.table <- x.table #generate table
chi.table <- chi.table[,-1] #remove redundant column name

chi.test <- chisq.test(chi.table) #run chi-square test
print(chi.test) #generate chiu-squared results

# Generate Table of changes in genera pre- and post-burn organized by Clade and Phyla
ps.phylum <- ps %>%
  tax_glom(taxrank = "Phylum")

otu.table <- as.matrix(otu_table(ps.phylum)) #Pull out abundance table
abundance.table <- otu.table@.Data

merged.table <- cbind(abundance.table, ps.env[,c("fire", "year")]) #Combine abundances with environmental data
merged.table <- as.data.frame(merged.table)

# Summarize by Fire Status and Pre-/Post- Burn
for(i in 1:32){
  merged.table[,i] <- as.numeric(merged.table[,i]) #Set as Numeric
}

merged.table$fire <- as.factor(merged.table$fire) #Set as Factor
merged.table$year <- as.factor(merged.table$year) #Set as Factor

merged.data <- aggregate(.~fire+year, merged.table, FUN = sum) #Create data frame to append too

# Add proper column names (environment and phyla)
colnames(merged.data) <- c("fire", "year", "Proteobacteria","Acidobacteriota","Actinobacteriota",
                           "Verrucomicrobiota","RCP2-54","Planctomycetota", "WPS-2",
                           "Bacteroidota","Chloroflexi","Myxococcota","Cyanobacteria",
                           "Methylomirabilota","Gemmatimonadota","Dependentiae",
                           "Patescibacteria","Firmicutes","Armatimonadota","Nitrospirota",
                           "Entotheonellaeota","MBNT15","Elusimicrobiota","Desulfobacterota",
                           "NB1-j","Bdellovibrionota","Latescibacterota",
                           "FCPU426","Fibrobacterota","SAR324 clade(Marine group B)",
                           "Deinococcota","Spirochaetota","Bacteria Kingdom", "GAL15")

#Separate Burned and Unburned Plots
burned.data <- merged.data[merged.data$fire == "Burned",]
unburned.data <- merged.data[merged.data$fire == "Unburned",]

#Calculate Change
rownames(burned.data) <- c("Pre-", "Post-")
burned.data <- burned.data[,-1]
burned.data <- burned.data[,-1]
burned.data["Change",] <- ((burned.data["Post-",] - burned.data["Pre-",])/burned.data["Pre-",])*100
burned.data <- round(burned.data)

rownames(unburned.data) <- c("Pre-", "Post-")
unburned.data <- unburned.data[,-1]
unburned.data <- unburned.data[,-1]
unburned.data["Change",] <- ((unburned.data["Post-",] - unburned.data["Pre-",])/unburned.data["Pre-",])*100
unburned.data <- round(unburned.data)

# Export Table of Burned Data
burned.data <- t(burned.data)
total.data <- colSums(burned.data)
total.data[3] <- ((total.data[2] - total.data[1])/total.data[1])*100
xtable(burned.data, type = "latex", file = paste0(getwd(), "/Tables/burnedchange.tex"))

# Export Table of Unburned Data
unburned.data <- t(unburned.data)
total.data <- colSums(unburned.data)
total.data[3] <- ((total.data[2] - total.data[1])/total.data[1])*100
xtable(unburned.data, type = "latex", file = paste0(getwd(), "/Tables/unburnedchange.tex"))

# N.B. Table in Manuscript also has Log-Ratio Change. This code has not been 
# generated in R yet, was originally done manually in LateX. Will be updated to
# reflect the manuscript.



## Supplementary Figures (microViz and Rocca Paper)

# Normality

# Beta Diversity PERMANOVA
pn <- ps.weighted %>%
  tax_agg("Genus") %>%
  dist_calc("bray") %>%
  dist_permanova(
    seed = 1,
    variables = c("pH","Fire", "Fall", "Spring", "Summer"),
    n_processes = 1,
    n_perms = 99)
    
print(pn)

# Heatmap
# Correlation Heatmap
cor_heatmap(ps.phylum)

# Composition Heatmap
comp <- comp_heatmap(ps.phylum,  sample_anno = sampleAnnotation(
  Site = anno_sample_cat("plotID.x", legend_title = "Site"),
  `Burn Status` = anno_sample_cat("fire", legend_title = "Burn Status"),
  Time = anno_sample_cat("Time", legend_title = "Time")
))

pdf(file = paste0(getwd(),"/Figures/FigureSCompHeatmap.pdf"))
comp.heatmap <- ComplexHeatmap::draw(
  object = comp,
  annotation_legend_list = attr(comp, "AnnoLegends"), merge_legends = TRUE)
dev.off()

# Beta-dispersion
# Reuse ps.weighted, with Fire, Summer, Fall, and Spring set as factors
ps.weighted@sam_data$Fire <- as.factor(ps.weighted@sam_data$Fire)
ps.weighted@sam_data$Fall <- as.factor(ps.weighted@sam_data$Fall)
ps.weighted@sam_data$Spring <- as.factor(ps.weighted@sam_data$Spring)
ps.weighted@sam_data$Summer <- as.factor(ps.weighted@sam_data$Summer)

bd <- ps.weighted %>%
  tax_agg("Genus") %>%
  dist_calc("bray") %>%
  dist_bdisp(variables = c("Fire", "Fall", "Spring", "Summer")) %>%
  bdisp_get()

print(bd)
# Iris Plot
# Use Phyla Clade as Grouping
# Filter to opnly post-burn plots
ord <- ps.taxa %>%
  subset_samples(fire == "Burned") %>%
  tax_fix(unknowns = c("Bacteria Candidate Group", "CCD", "CMS", "FCB", "Gracilicutes", "PANNAM", "PVC", "Terrabacteria")) %>%
  tax_agg("Phylum") %>%
  dist_calc("bray") %>%
  ord_calc(method = "PCoA")

op <- ord_plot_iris(
  data = ord,
  n_taxa = 32,
  tax_level = "Phylum",
  anno_colour = "Time",
  anno_colour_style = list(size = 3)) +
  scale_fill_manual(values=(mycolors), name="Phylum")

ggsave(filename= paste0(getwd(),"/Figures/FigureS3.pdf"), 
       plot = op, width = 180, height = 240, units=c("mm"), dpi=600) 

# PhyloD

