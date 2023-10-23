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

# load metadata (biological replicate information)
don <- read.csv('C:/Users/aggiej/Documents/Fellowship/Tissue culture/MSC_CM_in_vitro/Metadata.csv')

# add metadata to the layout dataset
layout <- layout %>% full_join(don, by = 'Donor')

# create a unique well identifier
# upper case identifier refers to position of PCR plate (qPCR technical replicate) 
# lower case refers to experimental technical replicate (e.g. tissue culture wells)
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

# link samples with sd Cp > 0.2 to the raw data
data_sd <- sum_sd %>% left_join(data, by = c('Gene', 'Sample'), multiple = 'all')
length(unique(data_sd$Sample))

# plot individual Cp data points from these samples
# separate plot for each gene
unique(data_sd$Gene)

p1 <- data_sd %>% group_by(Gene) %>%
  do(plots = ggplot(data = .) + (aes(x = Sample, y = CrossingPoint)) + geom_point() +
       ggtitle(.$Gene)) 

marrangeGrob(grobs = p1$plots, nrow = 3, ncol = 2)

# create vector of outlying Cp values based on graph and data inspection
outl_gap <- c('P1E8')
outl_col1 <- c('P6C7', 'P6C1', 'P6C2', 'P6C3', 'P1F10', 'P1F10')
outl_col3 <- c('P5B12', 'P1E5', 'P1E3')
outl_mmp1 <- c('P5H6', 'P5D10', 'P5D9', 'P5D3', 'P3G10', 'P3D2', 'P1F7')
outl_mmp3 <- c('P1B4', 'P2E10')

# remove outlying values from the database

gap_corr <- data %>%
  filter(Gene == 'GAP' & !(PlateWell %in% outl_gap))

col1_corr <- data %>%
  filter(Gene == 'COL1' & !(PlateWell %in% outl_col1))

col3_corr <- data %>%
  filter(Gene == 'COL3' & !(PlateWell %in% outl_col3))

mmp1_corr <- data %>%
  filter(Gene == 'MMP1' & !(PlateWell %in% outl_mmp1))

mmp3_corr <- data %>%
  filter(Gene == 'MMP3' & !(PlateWell %in% outl_mmp3))

data_corr <- rbind(gap_corr, col1_corr, col3_corr, mmp1_corr, mmp3_corr) 

# calculate mean Cp and standard deviation of each triplicate 
data_corr_sum <- data_corr %>%
  group_by(Gene, Sample) %>%
  summarise(mean_Cp = mean(CrossingPoint), sd_Cp = sd(CrossingPoint))

# filter and link with metadata
data_corr_sum <- data_corr_sum %>% filter(Sample != '')

data_corr_sum <- data_corr_sum[!grepl('^NoT', data_corr_sum$Sample),]

layout2 <- layout %>% select(!c('Plate','Well','PlateWell')) %>% unique()

data_clean <- data_corr_sum %>% left_join(layout2, by = c('Sample'))


# transform categorical variables from integers to factors
data_clean$T1 <- as.factor(data_clean$T1)
data_clean$T2 <- as.factor(data_clean$T2)

# plot distribution of gene raw Cp to assess the effect of experimental conditions
pdf('Raw_Cp.pdf')

p2 <- data_clean %>% group_by(Gene) %>%
  do(plots = ggplot(data = .) + (aes(x = T1, y = mean_Cp, colour = T2)) + geom_point() +
       ggtitle(.$Gene)) 

marrangeGrob(grobs = p2$plots, nrow = 1, ncol = 1)

dev.off()

# change format from long to wide
data_clean <- data_clean %>% select(!sd_Cp) %>%
  pivot_wider(names_from = Gene, values_from = mean_Cp)


# calculate delta Cp value for each gene
data_clean <- data_clean %>% 
  mutate(COL1_RQ = 2^-(COL1 - GAP)) %>%
  mutate(COL3_RQ = 2^-(COL3 - GAP)) %>%
  mutate(MMP1_RQ = 2^-(MMP1 - GAP)) %>%
  mutate(MMP3_RQ = 2^-(MMP3 - GAP))

# remove columns containing raw Cp value
data_clean <- data_clean %>% select(!c('GAP', 'COL1', 'COL3', 'MMP1', 'MMP3'))

# save processed data for statistical analysis
write.csv(data_clean, 'qPCR_MSC_CM_processed.csv')


