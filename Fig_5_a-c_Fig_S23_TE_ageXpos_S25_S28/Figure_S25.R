#### SUPPLEMENTARY FIGURE S25 ####
### This script contains the code for the comparison of insertion ages of CRMs in- and outside the centromere ###
library(SNPRelate)
library(ggplot2)
library(dplyr)
library(stringr)
library(gdsfmt)
library(magrittr)
library(tidyverse)
library(hrbrthemes)
library(ggridges)

#### read in data for TE families ####
GROUP_DIR <- ""
setwd(GROUP_DIR)

name <- "Avia Abigail Cereba Quinta"
species <- expression(paste(italic("Sercale cereale"), " Lo7v3"))
family <- "RLG_ScerLo7"
family_name <- expression(bold("Avia Abigail Cereba Quinta"))

ScerLo7_Cereba <- read.table("age_distribution_RLG_ScerLo7_Cereba")
ScerLo7_Cereba$family <- str_extract(ScerLo7_Cereba$V1, "\\S+(?<=_Cereba)")

ScerLo7_Quinta <- read.table("age_distribution_RLG_ScerLo7_Quinta")
ScerLo7_Quinta$family <- str_extract(ScerLo7_Quinta$V1, "\\S+(?<=_Quinta)")

ScerLo7_Abia <- read.table("age_distribution_RLG_ScerLo7_Abia")
ScerLo7_Abia$family <- str_extract(ScerLo7_Abia$V1, "\\S+(?<=_Abia)")

ScerLo7_Abigail <- read.table("age_distribution_RLG_ScerLo7_Abigail")
ScerLo7_Abigail$family <- str_extract(ScerLo7_Abigail$V1, "\\S+(?<=_Abigail)")


all_Lo7 <- rbind(ScerLo7_Abia, ScerLo7_Abigail, ScerLo7_Cereba, ScerLo7_Quinta)

all_Lo7$TEnum <- str_extract(all_Lo7$V1, "\\d\\S-\\d+")
all_Lo7$chr <- str_extract(all_Lo7$TEnum, "\\d\\S")
all_Lo7$TEstart <- as.numeric(str_extract(all_Lo7$TEnum, "(?<=-)\\d+"))*1000 #time 1000 to get the start position from the ID 
all_Lo7$genome <- as.factor(str_extract(all_Lo7$chr, "(?<=\\d)\\S"))

df_sorted <- all_Lo7 %>% arrange(genome)
df_sorted$chr <- factor(df_sorted$chr, levels = rev(unique(df_sorted$chr)))

hist(df_sorted$V6)

##### define location categories (centromere, <10Mb, >10Mb) #####
df_sorted_centromere_1R <- df_sorted[df_sorted$TEstart>296200000 & df_sorted$TEstart<306200000 & df_sorted$chr == "1R",]
df_sorted_centromere_2R <- df_sorted[df_sorted$TEstart>440500000 & df_sorted$TEstart<452500000 & df_sorted$chr == "2R",]
df_sorted_centromere_3R <- df_sorted[df_sorted$TEstart>525100000 & df_sorted$TEstart<537100000 & df_sorted$chr == "3R",]
df_sorted_centromere_4R <- df_sorted[df_sorted$TEstart>383500000 & df_sorted$TEstart<394900000 & df_sorted$chr == "4R",]
df_sorted_centromere_5R <- df_sorted[df_sorted$TEstart>286600000 & df_sorted$TEstart<298000000 & df_sorted$chr == "5R",]
df_sorted_centromere_6R <- df_sorted[df_sorted$TEstart>308400000 & df_sorted$TEstart<319000000 & df_sorted$chr == "6R",]
df_sorted_centromere_7R <- df_sorted[df_sorted$TEstart>436500000 & df_sorted$TEstart<448500000 & df_sorted$chr == "7R",]

all_centromere_chr <- rbind(df_sorted_centromere_1R,df_sorted_centromere_2R,df_sorted_centromere_3R,df_sorted_centromere_4R,
                            df_sorted_centromere_5R,df_sorted_centromere_6R,df_sorted_centromere_7R)
all_centromere_chr$category <- "centromere"

### Chromosome 1R ###
df_sorted_centromere_minus10Mb_1R <- df_sorted[df_sorted$TEstart>286200000 & df_sorted$TEstart<296200000 & df_sorted$chr == "1R",]
df_sorted_centromere_plus10Mb_1R <-  df_sorted[df_sorted$TEstart>306200000 & df_sorted$TEstart<316200000 & df_sorted$chr == "1R",]
df_sorted_centromere_10Mb_1R <- rbind(df_sorted_centromere_minus10Mb_1R, df_sorted_centromere_plus10Mb_1R)
df_sorted_centromere_10Mb_1R$category <- "<10Mb"

df_sorted_centromere_minus20Mb_1R <- df_sorted[df_sorted$TEstart>0 & df_sorted$TEstart<286200000 & df_sorted$chr == "1R",]
df_sorted_centromere_plus20Mb_1R <-  df_sorted[df_sorted$TEstart>316200000 & df_sorted$TEstart<max(df_sorted$TEstart) & df_sorted$chr == "1R",]
df_sorted_centromere_20Mb_1R <- rbind(df_sorted_centromere_minus20Mb_1R, df_sorted_centromere_plus20Mb_1R) 
df_sorted_centromere_20Mb_1R$category <- ">10Mb"

### Chromosome 2R ###
df_sorted_centromere_minus10Mb_2R <- df_sorted[df_sorted$TEstart>430500000 & df_sorted$TEstart<440500000 & df_sorted$chr == "2R",]
df_sorted_centromere_plus10Mb_2R <-  df_sorted[df_sorted$TEstart>452500000 & df_sorted$TEstart<462500000 & df_sorted$chr == "2R",]
df_sorted_centromere_10Mb_2R <- rbind(df_sorted_centromere_minus10Mb_2R, df_sorted_centromere_plus10Mb_2R)
df_sorted_centromere_10Mb_2R$category <- "<10Mb"

df_sorted_centromere_minus20Mb_2R <- df_sorted[df_sorted$TEstart>0 & df_sorted$TEstart<430500000 & df_sorted$chr == "2R",]
df_sorted_centromere_plus20Mb_2R <-  df_sorted[df_sorted$TEstart>462500000 & df_sorted$TEstart<max(df_sorted$TEstart) & df_sorted$chr == "2R",]
df_sorted_centromere_20Mb_2R <- rbind(df_sorted_centromere_minus20Mb_2R, df_sorted_centromere_plus20Mb_2R) 
df_sorted_centromere_20Mb_2R$category <- ">10Mb"

### Chromosome 3R ###
df_sorted_centromere_minus10Mb_3R <- df_sorted[df_sorted$TEstart>515100000 & df_sorted$TEstart<525100000 & df_sorted$chr == "3R",]
df_sorted_centromere_plus10Mb_3R <-  df_sorted[df_sorted$TEstart>537100000 & df_sorted$TEstart<547100000 & df_sorted$chr == "3R",]
df_sorted_centromere_10Mb_3R <- rbind(df_sorted_centromere_minus10Mb_3R, df_sorted_centromere_plus10Mb_3R)
df_sorted_centromere_10Mb_3R$category <- "<10Mb"

df_sorted_centromere_minus20Mb_3R <- df_sorted[df_sorted$TEstart>0 & df_sorted$TEstart<515100000 & df_sorted$chr == "3R",]
df_sorted_centromere_plus20Mb_3R <-  df_sorted[df_sorted$TEstart>547100000 & df_sorted$TEstart<max(df_sorted$TEstart) & df_sorted$chr == "3R",]
df_sorted_centromere_20Mb_3R <- rbind(df_sorted_centromere_minus20Mb_3R, df_sorted_centromere_plus20Mb_3R) 
df_sorted_centromere_20Mb_3R$category <- ">10Mb"

### Chromosome 4R ###
df_sorted_centromere_minus10Mb_4R <- df_sorted[df_sorted$TEstart>373500000 & df_sorted$TEstart<383500000 & df_sorted$chr == "4R",]
df_sorted_centromere_plus10Mb_4R <-  df_sorted[df_sorted$TEstart>394900000 & df_sorted$TEstart<404900000 & df_sorted$chr == "4R",]
df_sorted_centromere_10Mb_4R <- rbind(df_sorted_centromere_minus10Mb_4R, df_sorted_centromere_plus10Mb_4R)
df_sorted_centromere_10Mb_4R$category <- "<10Mb"

df_sorted_centromere_minus20Mb_4R <- df_sorted[df_sorted$TEstart>0 & df_sorted$TEstart<373500000 & df_sorted$chr == "4R",]
df_sorted_centromere_plus20Mb_4R <-  df_sorted[df_sorted$TEstart>404900000 & df_sorted$TEstart<max(df_sorted$TEstart) & df_sorted$chr == "4R",]
df_sorted_centromere_20Mb_4R <- rbind(df_sorted_centromere_minus20Mb_4R, df_sorted_centromere_plus20Mb_4R) 
df_sorted_centromere_20Mb_4R$category <- ">10Mb"

### Chromosome 5R ###
df_sorted_centromere_minus10Mb_5R <- df_sorted[df_sorted$TEstart>276600000 & df_sorted$TEstart<286600000 & df_sorted$chr == "5R",]
df_sorted_centromere_plus10Mb_5R <-  df_sorted[df_sorted$TEstart>298000000 & df_sorted$TEstart<308000000 & df_sorted$chr == "5R",]
df_sorted_centromere_10Mb_5R <- rbind(df_sorted_centromere_minus10Mb_5R, df_sorted_centromere_plus10Mb_5R)
df_sorted_centromere_10Mb_5R$category <- "<10Mb"

df_sorted_centromere_minus20Mb_5R <- df_sorted[df_sorted$TEstart>0 & df_sorted$TEstart<276600000 & df_sorted$chr == "5R",]
df_sorted_centromere_plus20Mb_5R <-  df_sorted[df_sorted$TEstart>308000000 & df_sorted$TEstart<max(df_sorted$TEstart) & df_sorted$chr == "5R",]
df_sorted_centromere_20Mb_5R <- rbind(df_sorted_centromere_minus20Mb_5R, df_sorted_centromere_plus20Mb_5R) 
df_sorted_centromere_20Mb_5R$category <- ">10Mb"

### Chromosome 6R ###
df_sorted_centromere_minus10Mb_6R <- df_sorted[df_sorted$TEstart>298400000 & df_sorted$TEstart<308400000 & df_sorted$chr == "6R",]
df_sorted_centromere_plus10Mb_6R <-  df_sorted[df_sorted$TEstart>319000000 & df_sorted$TEstart<329000000 & df_sorted$chr == "6R",]
df_sorted_centromere_10Mb_6R <- rbind(df_sorted_centromere_minus10Mb_6R, df_sorted_centromere_plus10Mb_6R)
df_sorted_centromere_10Mb_6R$category <- "<10Mb"

df_sorted_centromere_minus20Mb_6R <- df_sorted[df_sorted$TEstart>0 & df_sorted$TEstart<298400000 & df_sorted$chr == "6R",]
df_sorted_centromere_plus20Mb_6R <-  df_sorted[df_sorted$TEstart>329000000 & df_sorted$TEstart<max(df_sorted$TEstart) & df_sorted$chr == "6R",]
df_sorted_centromere_20Mb_6R <- rbind(df_sorted_centromere_minus20Mb_6R, df_sorted_centromere_plus20Mb_6R) 
df_sorted_centromere_20Mb_6R$category <- ">10Mb"

### Chromosome 7R ###
df_sorted_centromere_minus10Mb_7R <- df_sorted[df_sorted$TEstart>426500000 & df_sorted$TEstart<436500000 & df_sorted$chr == "7R",]
df_sorted_centromere_plus10Mb_7R <-  df_sorted[df_sorted$TEstart>448500000 & df_sorted$TEstart<458500000 & df_sorted$chr == "7R",]
df_sorted_centromere_10Mb_7R <- rbind(df_sorted_centromere_minus10Mb_7R, df_sorted_centromere_plus10Mb_7R)
df_sorted_centromere_10Mb_7R$category <- "<10Mb"

df_sorted_centromere_minus20Mb_7R <- df_sorted[df_sorted$TEstart>0 & df_sorted$TEstart<426500000 & df_sorted$chr == "7R",]
df_sorted_centromere_plus20Mb_7R <-  df_sorted[df_sorted$TEstart>458500000 & df_sorted$TEstart<max(df_sorted$TEstart) & df_sorted$chr == "7R",]
df_sorted_centromere_20Mb_7R <- rbind(df_sorted_centromere_minus20Mb_7R, df_sorted_centromere_plus20Mb_7R) 
df_sorted_centromere_20Mb_7R$category <- ">10Mb"

#combine the location categories 
df_sorted_centromere_10Mb <- rbind(df_sorted_centromere_10Mb_1R, df_sorted_centromere_10Mb_2R, df_sorted_centromere_10Mb_3R, df_sorted_centromere_10Mb_4R,
                                   df_sorted_centromere_10Mb_5R, df_sorted_centromere_10Mb_6R, df_sorted_centromere_10Mb_7R)
df_sorted_centromere_20Mb <- rbind(df_sorted_centromere_20Mb_1R, df_sorted_centromere_20Mb_2R, df_sorted_centromere_20Mb_3R, df_sorted_centromere_20Mb_4R,
                                   df_sorted_centromere_20Mb_5R, df_sorted_centromere_20Mb_6R, df_sorted_centromere_20Mb_7R)

all <- rbind(all_centromere_chr, df_sorted_centromere_10Mb, df_sorted_centromere_20Mb)
str(all)
any(duplicated(all$TEnum))#check if unique TEs

##### check if localisation category worked #####
all$chr <- factor(all$chr , levels = c("1R", "2R", "3R", "4R", "5R", "6R", "7R"))
all$category <- factor(all$category , levels = c("centromere", "<10Mb", ">10Mb"))

ggplot(data = all)+
  geom_point(aes(y=V6, x=TEstart/1000000, color= category))+
  ylab("Age [myr]") +
  xlab("Position") +
  labs(color = "category") +
  facet_grid(chr~.) +
  xlim(c(0, 1050)) +
  scale_color_manual(values = c("#02818a", "#67a9cf", "#bdc9e1")) +
  theme_bw()

##### define age categories (0-0.1, 0.1-0.2, ...) #####
all_cat <- all %>%
  mutate(age_cat = case_when(V6 >= 0 & V6 <=0.1 ~ "0-0.1", 
                           V6 > 0.1 & V6 <= 0.25  ~ "0.1-0.25",
                           V6 > 0.25 & V6 <= 0.5  ~ "0.25-0.5",
                           V6 > 0.5 & V6 <= 0.75  ~ "0.5-0.75",
                           V6 > 0.75 & V6 <= 1  ~ "0.75-1",
                           TRUE ~ ">1")) 

###### Plot full-length copies count over age category colored by location category ######
all_cat$age_cat <- factor(all_cat$age_cat, levels = c("0-0.1","0.1-0.25", "0.25-0.5", "0.5-0.75", "0.75-1", ">1"))
all_cat$category <- factor(all_cat$category, levels = c("centromere", "<10Mb", ">10Mb"))

ggplot(data = all_cat) +
    geom_bar(aes(x=age_cat, fill= category), position = "dodge") +
    scale_fill_manual(values = c("#02818a", "#67a9cf", "#bdc9e1")) +
  theme_ipsum(base_size = 20.5, axis_title_size = 25.5, grid_col = "black", axis_col = "black", strip_text_size = 25, ) +
  ylab("Full-length TE copies") +
  xlab("Age [myr]")

###### Plot full-length copies count over age category colored by location category and faceted by TE-family ###### 
#count the occurrence of matches for the three categories
plot_data <- all_cat %>%
  count(age_cat, category, family)

ggplot(plot_data, aes(x = age_cat, y = n, fill = category)) +
  geom_col(position = position_dodge2(preserve = "single", padding = 0.09), aes(group = category)) +
  theme_ipsum(base_size = 13, axis_title_size = 18, grid_col = "black", axis_col = "black", strip_text_size = 18, strip_text_face = "italic" ) +
  xlab("Insertion age [myr]") +
  ylab("Full-length TE copies") +
  facet_grid(~family) +
  scale_fill_manual(values = c("#02818a", "#67a9cf", "#bdc9e1")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

#count the occurrence of matches for the four categories
plot_data <- all_cat %>%
  count(age_cat, category, family, chr)

ggplot(plot_data, aes(x = age_cat, y = n, fill = category)) +
  geom_col(position = position_dodge2(preserve = "single", padding = 0.09), aes(group = category)) +
  theme_ipsum(base_size = 20.5, axis_title_size = 25.5, grid_col = "black", axis_col = "black", strip_text_size = 25, strip_text_face = "italic") +
  xlab("Age [myr]") +
  ylab("Full-length TE copies") +
  facet_grid(chr~family) +
  scale_fill_manual(values = c("#02818a", "#67a9cf", "#bdc9e1"))


##### calculate the percentage #####
(plot_data[plot_data$family == "RLG_ScerLo7_Abia" & plot_data$age_cat == "0-0.1" & plot_data$category == "centromere",]$n)/(plot_data[plot_data$family == "RLG_ScerLo7_Abia" & plot_data$age_cat == "0-0.1" & plot_data$category == "centromere",]$n + plot_data[plot_data$family == "RLG_ScerLo7_Abia" & plot_data$age_cat == "0-0.1" & plot_data$category == "<10Mb",]$n +
  plot_data[plot_data$family == "RLG_ScerLo7_Abia" & plot_data$age_cat == "0-0.1" & plot_data$category == ">10Mb",]$n)*100

