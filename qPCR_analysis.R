# Analysis of the effect of MSC CM on gene expression in equine tenocytes

# load packages
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(gridExtra)

# load processed qPCR dataset
data_clean <- read.csv('C:/Users/aggiej/Documents/Fellowship/Tissue culture/MSC_CM_in_vitro/qPCR/qPCR_MSC_CM_processed.csv')

# transform dummy variables from integers to factors
data_clean$T1 <- as.factor(data_clean$T1)
data_clean$T2 <- as.factor(data_clean$T2)

## col1
col1_rq <- data_clean %>%
  ggplot(aes(y = COL1_RQ, x = T1, colour = T2)) + geom_point() + ggtitle('COL1')
col1_rq

hist(data_clean$COL1_RQ)

# without extreme values
col1_rq <- data_clean %>%
  filter(COL1_RQ < 30) %>%
  ggplot(aes(y = COL1_RQ, x = T1, colour = T2)) + geom_point() + ggtitle('COL1')

col1_rq

## col3
col3_rq <- data_clean %>%
  ggplot(aes(y = COL3_RQ, x = T1, colour = T2)) + geom_point() + ggtitle('COL3')
col3_rq

hist(data_clean$COL1_RQ)

# without extreme values
col3_rq <- data_clean %>%
  filter(COL3_RQ <30) %>%
  ggplot(aes(y = COL3_RQ, x = T1, colour = T2)) + geom_point() + ggtitle('COL3')
col3_rq

## mmp1
mmp1_rq <- data_clean %>%
  ggplot(aes(y = MMP1_RQ, x = T1, colour = T2)) + geom_point() + ggtitle('MMP1')
mmp1_rq

hist(data_clean$MMP1_RQ)

## mmp3
mmp3_rq <- data_clean %>%
  ggplot(aes(y = MMP3_RQ, x = T1, colour = T2)) + geom_point() + ggtitle('MMP3')
mmp3_rq

hist(data_clean$MMP3_RQ)




#
#
#


# MODELLING EFFECT OF TREATMENTS ON GENE EXPRESSION
install.packages('emmeans')
library(emmeans)

## col1
model1 <- glm(COL1_RQ ~ T1*T2,
              data = data_clean,
              family = Gamma(link = 'identity'))
summary(model1)

# pairwise comparison conditioning on presence/absence of inflammation
em1 <- emmeans(model1, pairwise ~ T2 | T1)
em1



## col3
model2 <- glm(COL3_RQ ~ T1*T2,
              data = data_clean,
              family = Gamma)
summary(model2)

# pairwise comparison
em2 <- emmeans(model2, pairwise ~ T2 | T1)
em2



## mmp1
model3 <- glm(MMP1_RQ ~ T1*T2,
              data = data_clean,
              family = Gamma(link = 'identity'))
summary(model3)



## mmp3
model4 <- glm(MMP3_RQ ~ T1*T2,
              data = data_clean,
              family = Gamma)
summary(model4)


#
#
#

# EFFECT IN DIFFERENT AGE GROUPS
old <- data_clean %>% filter(RecAge == 'O')
young <- data_clean %>% filter(RecAge == 'Y')

# col1
col1_o <- old %>%
  ggplot(aes(y = COL1_RQ, x = T1, colour = T2)) + geom_jitter() + ggtitle('COL1 Old')
col1_o

col1_y <- young %>%
  ggplot(aes(y = COL1_RQ, x = T1, colour = T2)) + geom_jitter() + ggtitle('COL1 Young')
col1_y

model1_o <- glm(COL1_RQ ~ T1*T2,
                data = old,
                family = Gamma(link = 'identity'))
summary(model1_o)

model1_y <- glm(COL1_RQ ~ T1*T2,
                data = young,
                family = Gamma(link = 'identity'))
summary(model1_y)

# col3
col1_o <- old %>%
  ggplot(aes(y = COL3_RQ, x = T1, colour = T2)) + geom_jitter() + ggtitle('COL3 Old')
col1_o

col1_y <- young %>%
  ggplot(aes(y = COL3_RQ, x = T1, colour = T2)) + geom_jitter() + ggtitle('COL3 Young')
col1_y

model2_o <- glm(COL3_RQ ~ T1*T2,
                data = old,
                family = Gamma(link = 'identity'))
summary(model2_o)

model2_y <- glm(COL3_RQ ~ T1*T2,
                data = young,
                family = Gamma(link = 'identity'))
summary(model2_y)

# mmp1
mmp1_o <- old %>%
  ggplot(aes(y = MMP1_RQ, x = T1, colour = T2)) + geom_jitter() + ggtitle('MMP1 Old')
mmp1_o

mmp1_y <- young %>%
  ggplot(aes(y = MMP1_RQ, x = T1, colour = T2)) + geom_jitter() + ggtitle('MMP1 Young')
mmp1_y

model3_o <- glm(MMP1_RQ ~ T1*T2,
                data = old,
                family = Gamma(link = 'identity'))
summary(model3_o)

model3_y <- glm(MMP1_RQ ~ T1*T2,
                data = young,
                family = Gamma(link = 'identity'))
summary(model3_y)

# mmp3
mmp3_o <- old %>%
  ggplot(aes(y = MMP3_RQ, x = T1, colour = T2)) + geom_jitter() + ggtitle('MMP3 Old')
mmp3_o

mmp3_y <- young %>%
  ggplot(aes(y = MMP3_RQ, x = T1, colour = T2)) + geom_jitter() + ggtitle('MMP3 Young')
mmp3_y

model4_o <- glm(MMP3_RQ ~ T1*T2,
                data = old,
                family = Gamma(link = 'identity'))
summary(model4_o)

em4 <- emmeans(model4_o, pairwise ~ T2 | T1)
em4


model4_y <- glm(MMP3_RQ ~ T1*T2,
                data = young,
                family = Gamma(link = 'identity'))
summary(model4_y)

#
#
#

# Effect of CTHRC1 level in treated group only
data_cm <- data_clean %>% filter(T2 == 1)

# col1
col1_c <- data_cm %>% ggplot(aes(y = COL1_RQ, x = CTHRC1, colour = T1)) +
  geom_point() + ggtitle('COL1') + geom_smooth()

col1_co <- old %>% ggplot(aes(y = COL1_RQ, x = CTHRC1, colour = T1)) +
  geom_point() + ggtitle('Old') + geom_smooth()

col1_cy <- young %>% ggplot(aes(y = COL1_RQ, x = CTHRC1, colour = T1)) +
  geom_point() + ggtitle('Young') + geom_smooth()

grid.arrange(col1_c, col1_co, col1_cy, ncol = 3)


# col3
col3_c <- data_cm %>% ggplot(aes(y = COL3_RQ, x = CTHRC1, colour = T1)) +
  geom_point() + ggtitle('COL3') + geom_smooth()

col3_co <- old %>% ggplot(aes(y = COL3_RQ, x = CTHRC1, colour = T1)) +
  geom_point() + ggtitle('Old') + geom_smooth()

col3_cy <- young %>% ggplot(aes(y = COL3_RQ, x = CTHRC1, colour = T1)) +
  geom_point() + ggtitle('Young') + geom_smooth()

grid.arrange(col3_c, col3_co, col3_cy, ncol = 3)

# mmp1
mmp1_c <- data_cm %>% ggplot(aes(y = MMP1_RQ, x = CTHRC1, colour = T1)) +
  geom_point() + ggtitle('MMP1') + geom_smooth()

mmp1_co <- old %>% ggplot(aes(y = MMP1_RQ, x = CTHRC1, colour = T1)) +
  geom_point() + ggtitle('Old') + geom_smooth()

mmp1_cy <- young %>% ggplot(aes(y = MMP1_RQ, x = CTHRC1, colour = T1)) +
  geom_point() + ggtitle('Young') + geom_smooth()

grid.arrange(mmp1_c, mmp1_co, mmp1_cy, ncol = 3)

# mmp3
mmp3_c <- data_cm %>% ggplot(aes(y = MMP3_RQ, x = CTHRC1, colour = T1)) +
  geom_point() + ggtitle('MMP3') + geom_smooth()

mmp3_co <- old %>% ggplot(aes(y = MMP3_RQ, x = CTHRC1, colour = T1)) +
  geom_point() + ggtitle('Old') + geom_smooth()

mmp3_cy <- young %>% ggplot(aes(y = MMP3_RQ, x = CTHRC1, colour = T1)) +
  geom_point() + ggtitle('Young') + geom_smooth()

grid.arrange(mmp3_c, mmp3_co, mmp3_cy, ncol = 3)
