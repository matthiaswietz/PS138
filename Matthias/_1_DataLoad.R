
#############################################################
  ###  ARCWATCH - AMPLICON ANALYSIS  ###
#############################################################

# This script: load and format dada2-amplicons 
# Load and format metadata

#setwd("Y:/AWI_MPI/collaborations/Beluga/Rstats")
setwd("C:/Users/mwietz/ownCloud - mwietz@owncloud.mpi-bremen.de/AWI_MPI/FRAM/Arcwatch/paper_Amplicons/Rstats")
#load("iNEXT.Rdata")

library(dplyr)
library(tibble)
library(tidyr)
library(gtools)
library(readr)
library(stringr)
library(ggplot2)
library(forcats)
library(purrr)
library(ampvis2)
#library(phyloseq)
library(vegan)
library(UpSetR)
library(mixOmics)
library(amap)
library(ape)
#library(lubridate)
#library(ANCOMBC)
#library(FEAST)
#library(iNEXT)
library(olsrr)
#library(rstatix)
#library(openair)
#library(ggOceanMaps)
#library(ggspatial)
library(solartime)


#############################################################
### METADATA -- Arcwatch ###
#############################################################

ENV.cao <- read.csv(
  "./metadata/envArcwatch.txt", h=T, 
  sep="\t", stringsAsFactors=F, skipNul=T) %>%
  distinct() %>%
  mutate(
    date=as.Date(date, format="%Y-%m-%d"),
    study="Arcwatch") %>%
   mutate_if(is.numeric, round, 2)

#############################################################
### METADATA -- FRAM ###
#############################################################

info.fram <- read.csv(
  "./metadata/sampleInfo.txt", h=T, 
  sep="\t", stringsAsFactors=F, skipNul=T) %>%
  #filter(mooring!="NK") %>%
  filter(!RAS_id %in% c(
    "03_2017_F4_1","04_2017_EGC_2")) %>%
  dplyr::select(c(
    "RAS_id","lat","locus_tag",
    "date1","date2","date3","date4",
    "mooringFull","mooring")) %>%
  reshape2::melt(id.vars=c(
    "locus_tag","mooringFull",
    "RAS_id","lat","mooring")) %>%
  dplyr::rename(
    sample_id = RAS_id,
    date = value) %>%
  mutate(date=as.Date(date)) 

# Calculate mean date
meanDate <- read.csv(
  "./metadata/sampleInfo.txt", 
  h=T, sep="\t", stringsAsFactors=F, skipNul=T) %>%
  #filter(mooring!="NK") %>%
  filter(!RAS_id %in% c(
    "03_2017_F4_1","04_2017_EGC_2")) %>%
  mutate_at(vars(
    date1, date2, date3, date4),
    as.Date, format="%Y-%m-%d") %>%
  rowwise %>%
  mutate(date = mean.Date(c(
    date1, date2, date3, date4), na.rm=T)) %>%
  #mutate(date = as.character(date)) %>%
  as.data.frame() %>%
  rename(sample_id=RAS_id)

## CTD 2016-2020 ##
CTD1 <- read.table(
  "./metadata/CTD.txt", 
  h=T, sep="\t", stringsAsFactors=F, skipNul=T)

## CTD 2021-2022 ##
# only depth <100m (table contains the two deeper RAS too)
# average by date
CTD2 <- read.table(
  "./metadata/CTD_2122.txt", sep="\t", header=T) %>%
  filter(depth < 100) %>%
  group_by(date, mooringFull) %>%
  summarize(across(where(
    is.numeric), mean, na.rm=T), .groups="drop") 

# Join
CTD.fram <- rbind(CTD1, CTD2)  %>%
  mutate(date=as.Date(date))  %>% 
  mutate(mooring = case_when(
    grepl("^F4", mooringFull) ~ sub("-.*", "", mooringFull),
    grepl("^EGC", mooringFull) ~ sub("-.*", "", mooringFull),
    grepl("^Fevi", mooringFull) ~ sub("-.*", "", mooringFull),
    grepl("^HG-IV-", mooringFull) ~ "HG-IV",
    TRUE ~ NA_character_))

# ICE COVER ##
# add column "Fevi" as duplicate of HG-VI (same mooring at 60m)
Ice.fram <- read.table(
  "./metadata/IceConc.txt",
  h=T, sep = "\t", check.names=F) %>%
  mutate(
    mm = sprintf("%02d", mm), # Add leading zero 
    dd = sprintf("%02d", dd), # Add leading zero 
    date = paste(yyyy, mm, dd, sep = "-")) %>%
  dplyr::select(-c(yyyy, mm, dd)) %>%
  mutate(Fevi=`HG-IV`) %>%
  mutate(date=as.Date(date)) %>%
  reshape2::melt(id.vars=c("date")) %>%
  dplyr::rename(
    mooring = variable, 
    iceConc = value) %>% 
  left_join(dplyr::select(
    info.fram, sample_id, mooring, date, variable)) %>%
  distinct(date, mooring, .keep_all=T)

# Mean iceConc per sampling
Ice.avg <- Ice.fram %>%
  group_by(sample_id) %>%
  summarise(iceConc=mean(iceConc)) %>%
  ungroup

# PAST ICE 
IcePast <- Ice.fram[!duplicated(
  Ice.fram$sample_id, incomparables=NA),] %>%
  fill(sample_id, .direction='up') %>%
  group_by(sample_id) %>%
  dplyr::summarize(icePast=mean(iceConc)) %>%
  drop_na() 

####################################################

# Join all
# Round numerics
# Convert NaN to NA
# Subset to moorings EGC and F4 (plus NegCtr)
# rename for consistency with Arcwatch
ENV.fram <- info.fram %>%
  left_join(CTD.fram) %>%
  left_join(Ice.fram) %>%
  left_join(IcePast) %>%
  distinct() %>%
  group_by(sample_id, locus_tag) %>%
  summarize_if(is.numeric, mean, na.rm=T) %>%
  ungroup() %>% 
  mutate_if(is.numeric, ~ifelse(is.nan(.), NA, .)) %>%
  left_join(meanDate) %>%
  mutate(date=as.Date(date, format="%Y-%m-%d")) %>%
  filter(mooring %in% c("EGC","F4","NegCtr") & locus_tag=="16S") %>%
  #dplyr::select(-c(
  #  "date1","date2","date3","date4","sig","depth",
  #  "CO2","pH","pcr_primer","target_subfragment",
   # "DNA_ng_uL","year","miseq_id","read_fwd",
   # "read_rev","read_fwd_run2","read_rev_run2")) %>%
  mutate(
    gear="RAS",
    study="FRAM") %>%
  mutate_if(is.numeric, round, 2) %>%
  dplyr::rename(
    InsituTemperature = temp,
    PracticalSalinity = sal,
    numeric_depth = depth) %>%
  mutate(Watermass1 = case_when(
    AW_frac > 0.7~"AW",
    PW_frac > 0.7~"PW",
    TRUE~"mixed"))

#############################################################
### COMBINE + FORMAT METADATA ###

# merge; move AWI_id as 1st column (for Ampvis)
ENV <- bind_rows(
  ENV.cao, ENV.fram) %>%
  dplyr::select(sample_id, everything())


#############################################################
  ### ASVs + TAXONOMY ###
#############################################################

ASV <- read.table(
  "../dada/PS138_16S_seqtab_wFRAM.txt",
  h = T,
  sep = "\t",
  check.names = F) %>% 
  filter(rowSums(.>= 3) > 3)

TAX <- read.table(
  "../dada/PS138_16S_tax_wFRAM.txt",
  h = T, 
  sep = "\t", 
  stringsAsFactors = F, 
  row.names = 1)

# match with subsetted ASVs
TAX <- TAX[row.names(ASV),]

# Rename BAC-NAs with last known taxrank + "uc"
k <- ncol(TAX)-1
for (i in 2:k) {
  if (sum(is.na(TAX[, i])) >1) {
    temp <- TAX[is.na(TAX[, i]), ]
    for (j in 1:nrow(temp)) {
      if (sum(is.na(
        temp[j, i:(k+1)])) == length(temp[j, i:(k+1)])) {
        temp[j, i] <- paste(temp[j, (i-1)], " uc", sep = "")
        temp[j, (i+1):(k+1)] <- temp[j, i]
      }
    }
    TAX[is.na(TAX[, i]), ] <- temp}
  if (sum(is.na(TAX[, i]))==1) {
    temp <- TAX[is.na(TAX[, i]), ]
    if (sum(is.na(temp[i:(k+1)])) == length(temp[i:(k+1)])) {
      temp[i] <- paste(temp[(i-1)], " uc", sep="")
      temp[(i+1):(k+1)] <- temp[i]
    }
    TAX[is.na(TAX[, i]),] <- temp
  }
}
TAX[is.na(TAX[, (k+1)]), (k+1)] <- paste(
  TAX[is.na(TAX[, (k+1)]), k], " uc", sep="")

# shorten/modify names
TAX <- TAX %>%
  mutate(across(everything(),~gsub("Clade ","SAR11 Clade ", .))) %>%
  mutate(across(everything(),~gsub(" clade","", .))) %>%
  mutate(across(everything(),~gsub("Candidatus","Cand", .))) %>%
  mutate(across(everything(),~gsub("Roseobacter NAC11-7 lineage","NAC11-7", .))) %>%
  mutate(across(everything(),~gsub(" marine group","", .))) %>%
  mutate(across(everything(),~gsub(" terrestrial group","", .))) %>%
  mutate(across(everything(),~gsub("_CC9902","", .))) %>%
  mutate(across(everything(),~gsub("(Marine group B)","", ., fixed=T))) %>%
  mutate(across(everything(),~gsub("(SAR406)","SAR406", ., fixed=T)))

#######################################

### CONTAMINANT CHECK -- Arcwatch

caoSamples <- ENV %>% filter(
  grepl("PS138", PANGAEA_event) & !grepl("NegCtr", sample_id)) %>% 
  pull(clip_id)

# Average negative controls
caoNeg <- ENV %>% filter(
  grepl("PS138-neg", sample_title)) %>% 
  pull(clip_id)

caoNeg <- round(rowMeans(ASV[, caoNeg]), 0)

# Subtract negative counts 
asv1 = ASV[caoSamples] - caoNeg

##################################

### CONTAMINANT CHECK -- RAS

# Define samples PCR-amplified with 25/30/35 cycles
ras.cyc25 <- ENV.fram %>% filter(
  cycles=="25" & !grepl("NegCtr", sample_id)) %>% 
  pull(clip_id)
ras.cyc30 <- ENV.fram %>% filter(
  cycles=="30" & !grepl("NegCtr", sample_id)) %>% 
  pull(clip_id)
ras.cyc35 <- ENV.fram %>% filter(
  cycles=="35" & !grepl("NegCtr", sample_id)) %>% 
  pull(clip_id)

# Corresponding negative controls
ras.neg25 <- ENV.fram %>% filter(
  cycles=="25" & grepl("NegCtr", sample_id)) %>% 
  pull(clip_id)
ras.neg30 <- ENV.fram %>% filter(
  cycles=="30" & grepl("NegCtr", sample_id)) %>% 
  pull(clip_id)
ras.neg35 <- ENV.fram %>% filter(
  cycles=="35" & grepl("NegCtr", sample_id)) %>% 
  pull(clip_id)

# Calculate means per controls; round
ras.neg25 <- round(rowMeans(ASV[, ras.neg25]), 0)
ras.neg30 <- round(rowMeans(ASV[, ras.neg30]), 0)
ras.neg35 <- round(rowMeans(ASV[, ras.neg35]), 0)

# Subtract negative counts 
asv2 = ASV[ras.cyc25] - ras.neg25
asv3 = ASV[ras.cyc30] - ras.neg30
asv4 = ASV[ras.cyc35] - ras.neg35

#############################

# Rejoin everything
ASV <- cbind(asv1, asv2, asv3, asv4)

# Set neg to zero
ASV[ASV < 0] <- 0

##################################

# Remove Mitochondria / chloroplasts
# Remove further potential contaminants 
TAX <- TAX %>% filter(
    !Family %in% c(
      "Burkholderiaceae", "Weeksellaceae", 
      "Neisseriaceae","Dermacoccaceae",
      "Corynebacteriaceae", "Bacillaceae",
      "Xanthobacteraceae", "Enterococcaceae",
      "Xanthomonadaceae","Carnobacteriaceae",
      "Streptococcaceae","Propionibacteriaceae",
      "Beijerinckiaceae","Mitochondria") &
    !Order %in% c(
      "Chloroplast","Peptostreptococcales-Tissierellales") & 
    !Class %in% c("Bacilli","Clostridia"))

# match ASVs after contaminant removal
ASV <- ASV[row.names(TAX),]

# remove NegCtr columns
#ASV2 <- ASV[, !grepl('negativ|control|NK', names(ASV), ignore.case=T)]

################################

# Factorize relevant parameters
ENV$depth_category <- factor(ENV$depth_category, levels=rev(c(
  "MP","Ice-TS","Ice-BS","UIW","2","10", "chl-max","25", "50", "100","200","500",
  "1000","1500","2000","3000","bottom+20","bottom+5","bottom","vent-plume")))

# Sort  
ASV <- ASV[,mixedsort(names(ASV))]
ENV <- ENV[mixedorder(ENV$clip_id),] %>%
  filter(!grepl("NegCtr", sample_id))

# Match and rename 
ASV <- ASV[,c((match(
  ENV$clip_id, colnames(ASV))))]
colnames(ASV) = ENV$sample_id


#############################################################
 ### TRANSFORM ###
#############################################################

# Relative abundances; round 2 digits
ASV.rel <- as.data.frame(apply(
  ASV, 2, function(x) round(x / sum(x) * 100, 3)))
ASV.hel = as.data.frame(
  apply(ASV, 2, function(x) sqrt(x / sum(x))))

# Ampvis load
ampvis <- data.frame(
  OTU = rownames(ASV),
  ASV, TAX, check.names=F)
ampvis <- amp_load(ampvis, ENV)

# Distance matrix -- all taxa
dist <- t(ASV.hel) %>% 
 as.data.frame %>% 
 vegdist(., method="bray")

#############################################################

# remove temp-data
rm(temp, i, j, k)

