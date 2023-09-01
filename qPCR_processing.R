## Processing qPCR data from Roche LightCycler 480

# load packages
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(gridExtra)

# load result files from qPCR experiments
pcr1 <- read.csv('C:/Users/aggiej/Documents/Fellowship/Tissue culture/MSC_CM_in_vitro/qPCR/PCR1.csv')
pcr2 <- read.csv('C:/Users/aggiej/Documents/Fellowship/Tissue culture/MSC_CM_in_vitro/qPCR/PCR2.csv')
pcr3 <- read.csv('C:/Users/aggiej/Documents/Fellowship/Tissue culture/MSC_CM_in_vitro/qPCR/PCR3.csv')
pcr4 <- read.csv('C:/Users/aggiej/Documents/Fellowship/Tissue culture/MSC_CM_in_vitro/qPCR/PCR4.csv')
pcr5 <- read.csv('C:/Users/aggiej/Documents/Fellowship/Tissue culture/MSC_CM_in_vitro/qPCR/PCR5.csv')

# combine results into one file
pcr <- rbind(pcr1, pcr2, pcr3, pcr4, pcr5)

# leave only relevant columns
pcr <- pcr %>% select('Experiment.Name', 'Position', 'CrossingPoint')

# remove rows with NA Cp values
pcr <- drop_na(pcr)

# extract information from the first column into new columns
pcr <- pcr %>% separate(Experiment.Name, c('Study', 'Gene', 'Plate'), '_')

# drop the redundant column (all data comes from one study)
pcr <- pcr %>% select(!'Study')

# load plate layout data (sample identifiers)
layout <- read.csv('C:/Users/aggiej/Documents/Fellowship/Tissue culture/MSC_CM_in_vitro/qPCR/Plate_layouts.csv')

# load donor data
don <- read.csv('C:/Users/aggiej/Documents/Fellowship/Tissue culture/MSC_CM_in_vitro/Metadata.csv')

# add donor information to the layout dataset
layout <- layout %>% full_join(don, by = 'Donor')

# create a unique well identifier
# upper case identifier refers to position of PCR plate (qPCR technical replicate) 
# lower case on cell culture plate (TC technical replicate)
pcr <- pcr %>% mutate(PlateWell = paste(Plate, Position, sep = ''))
layout <- layout %>% mutate(PlateWell = paste(Plate, Well, sep = ''))

# join qPCR results with layout data
data <- pcr %>% full_join(layout, by = c('PlateWell', 'Plate', 'Position' = 'Well'))

# calculate mean Cp and standard deviation of each triplicate
data_sum <- data %>%
  group_by(Gene, Sample) %>%
  summarise(mean_Cp = mean(CrossingPoint), sd_Cp = sd(CrossingPoint))

# remove row representing empty wells (no Sample information)
data_sum <- data_sum %>% filter(Sample != '')

# subset data from negative controls and assess Cp
data_neg <- data_sum %>% filter(grepl('^NoT', Sample))

# remove negative controls from database
data_sum <- data_sum[!grepl('^NoT', data_sum$Sample),]

# subset samples with sd Cp > 0.2 (conventional quality threshold)
sum_sd <- data_sum %>% filter(sd_Cp > 0.2)

# link samples to the raw data
data_sd <- sum_sd %>% left_join(data, by = c('Gene', 'Sample'), multiple = 'all')
length(unique(data_sd$Sample))

# plot individual Cp data points from these samples
# separate plot for each gene (similar Cp range)
unique(data_sd$Gene)

col1 <- data_sd %>% 
  filter(Gene == 'COL1') %>% 
  ggplot(aes(x = Sample, y = CrossingPoint)) + 
  geom_point() +
  ggtitle('COL1')

col3 <- data_sd %>% 
  filter(Gene == 'COL3') %>% 
  ggplot(aes(x = Sample, y = CrossingPoint)) + 
  geom_point() +
  ggtitle('COL3')

gap <- data_sd %>% 
  filter(Gene == 'GAP') %>% 
  ggplot(aes(x = Sample, y = CrossingPoint)) + 
  geom_point() +
  ggtitle('GAPDH')

mmp1 <- data_sd %>% 
  filter(Gene == 'MMP1') %>% 
  ggplot(aes(x = Sample, y = CrossingPoint)) + 
  geom_point() +
  ggtitle('MMP1')

mmp3 <- data_sd %>% 
  filter(Gene == 'MMP3') %>% 
  ggplot(aes(x = Sample, y = CrossingPoint)) + 
  geom_point() +
  ggtitle('MMP3')

grid.arrange(col1, col3, gap, mmp1, mmp3, nrow = 2)

# create vector of outlying Cp values based on graph inspection
gap
col1
col3
mmp1
mmp3

outl_gap <- c('P1E8')
outl_col1 <- c('P6C7', 'P6C1', 'P6C2', 'P6C3', 'P1F10', 'P1F10')
outl_col3 <- c('P5B12', 'P1E5', 'P1E3')
Outl_mmp1
Outl_mmp3 <- c('P1B4', 'P2E10')