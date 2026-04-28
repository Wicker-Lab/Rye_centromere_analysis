#### SUPPLEMENTARY FIGURE S28 ####
#### This script contains the code for the comparison Weining H1 and H2 to Lo7 ####
library(tidyverse)
library(ggridges)
library(ggplot2)
library(dplyr)
library(stringr)
library(magrittr)
library(ggraph)
library(patchwork)

custom_theme2 <- theme_classic(base_size = 10) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(color = "black"),
    axis.line = element_line(color = "black"),
    strip.text = element_text(size = 12),
    panel.grid.major = element_line(color = "lightgray"),  # add back grid
    panel.grid.minor = element_line(color = "lightgray"),   # optional
    panel.background = element_rect(fill = NA, color = NA),
    plot.background = element_rect(fill = NA, color = NA)
  )

OUTPUT_DIR <- ""
DATA_DIR <- ""
setwd(DATA_DIR)

#### SUPPLEMENTARY 28 PANEL A-E ####
#### chromosome 7: Weining H1 x Lo7 ####
##### Preparation Weining #####
setwd(DATA_DIR)
Weining_H1_Cereba <- read.table("age_distribution_LTR_Cereba_H1", skip = 1)
Weining_H1_Cereba$family <- str_extract(Weining_H1_Cereba$V1, "\\S+(?<=_Cereba)")

Weining_H1_Quinta <- read.table("age_distribution_LTR_Quinta_H1", skip = 1)
Weining_H1_Quinta$family <- str_extract(Weining_H1_Quinta$V1, "\\S+(?<=_Quinta)")

Weining_H1_Abia <- read.table("age_distribution_LTR_Abia_H1", skip = 1)
Weining_H1_Abia$family <- str_extract(Weining_H1_Abia$V1, "\\S+(?<=_Abia)")

Weining_H1_Abigail <- read.table("age_distribution_LTR_Abigail_H1", skip = 1)
Weining_H1_Abigail$family <- str_extract(Weining_H1_Abigail$V1, "\\S+(?<=_Abigail)")

all_Weining_H1 <- rbind(Weining_H1_Abia, Weining_H1_Abigail, Weining_H1_Cereba, Weining_H1_Quinta)

all_Weining_H1$TEnum <- str_extract(all_Weining_H1$V1, "\\d\\S-\\d+")
all_Weining_H1$chr <- str_extract(all_Weining_H1$TEnum, "\\d\\S")
all_Weining_H1$TEstart <- as.numeric(str_extract(all_Weining_H1$TEnum, "(?<=-)\\d+"))*1000 #time 1000 to get the start position from the ID 
all_Weining_H1$genome <- as.factor(str_extract(all_Weining_H1$chr, "(?<=\\d)\\S"))

df_sorted_H1 <- all_Weining_H1 %>% arrange(genome)
df_sorted_H1$chr <- factor(df_sorted_H1$chr, levels = rev(unique(df_sorted_H1$chr)))

chr7R_Weining_H1 <- df_sorted_H1[df_sorted_H1$chr == "7R",]
str(chr7R_Weining_H1)
chr7R_Weining_H1$V6 <- as.numeric(chr7R_Weining_H1$V6)

##### gaps #####
gaps <- read.table("Scer_Weining_v3_H1_chr7R_gaps", sep = ";")
gaps$start <- round(gaps$V1/1000000, 2)
gaps$end <- round(gaps$V2/1000000, 2)
gaps <- gaps[gaps$start>420 & gaps$end<490,]

p_right <- ggplot(data = chr7R_Weining_H1, aes(y = TEstart/1000000, x = V6, color = family)) +
  geom_point(shape=20, size=0.8)  +
  xlab("Age [myr]") +
  ylab("Position [Mb] H1") +
  labs(color = "Family") +
  scale_color_manual(values=c( "#e66101", "#fdb863","#5e3c99", "#b2abd2")) + 
  scale_x_continuous(limits = c(0,3), breaks = scales::pretty_breaks(n = 10)) +
  custom_theme2 +
  scale_y_reverse(limits = c(485, 435), breaks = scales::pretty_breaks(n = 10)) +
  geom_rect(data = gaps, aes(xmin = 2.8, xmax = 3, ymin = start, ymax = end), inherit.aes = FALSE, fill = NA, color = "black", size = 0.3)+
  geom_hline(yintercept = 457.00, color = "black", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 462.40, color = "black", linetype = "dashed", size = 0.5)
p_right

##### Preparation Lo7 #####
Lo7_Cereba <- read.table("age_distribution_RLG_ScerLo7_Cereba")
Lo7_Cereba$family <- str_extract(Lo7_Cereba$V1, "\\S+(?<=_Cereba)")

Lo7_Quinta <- read.table("age_distribution_RLG_ScerLo7_Quinta")
Lo7_Quinta$family <- str_extract(Lo7_Quinta$V1, "\\S+(?<=_Quinta)")

Lo7_Abia <- read.table("age_distribution_RLG_ScerLo7_Abia")
Lo7_Abia$family <- str_extract(Lo7_Abia$V1, "\\S+(?<=_Abia)")

Lo7_Abigail <- read.table("age_distribution_RLG_ScerLo7_Abigail")
Lo7_Abigail$family <- str_extract(Lo7_Abigail$V1, "\\S+(?<=_Abigail)")

all_Lo7 <- rbind(Lo7_Abia, Lo7_Abigail, Lo7_Cereba, Lo7_Quinta)

all_Lo7$TEnum <- str_extract(all_Lo7$V1, "\\d\\S-\\d+")
all_Lo7$chr <- str_extract(all_Lo7$TEnum, "\\d\\S")
all_Lo7$TEstart <- as.numeric(str_extract(all_Lo7$TEnum, "(?<=-)\\d+"))*1000 #time 1000 to get the start position from the ID 
all_Lo7$genome <- as.factor(str_extract(all_Lo7$chr, "(?<=\\d)\\S"))

df_sorted_H2 <- all_Lo7 %>% arrange(genome)
df_sorted_H2$chr <- factor(df_sorted_H2$chr, levels = rev(unique(df_sorted_H2$chr)))

chr7R_Lo7 <- df_sorted_H2[df_sorted_H2$chr == "7R",]
str(chr7R_Lo7)
chr7R_Lo7$V6 <- as.numeric(chr7R_Lo7$V6)

##### gaps #####
gaps <- read.table("Scer_Lo7_V3_chr7R_gaps", sep = ";")
gaps$start <- round(gaps$V1/1000000, 2)
gaps$end <- round(gaps$V2/1000000, 2)
gaps <- gaps[gaps$start>425 & gaps$end<490,]

p_top <- ggplot(data = chr7R_Lo7, aes(x = TEstart/1000000, y = V6, color = family)) +
  geom_point(shape=20, size=0.8)  +
  ylab("Age [myr]") +
  xlab("Position [Mb] Lo7") +
  labs(color = "Family") +
  scale_color_manual(values=c("#e66101", "#fdb863", "#5e3c99", "#b2abd2")) + 
  scale_x_continuous(limits = c(430, 480), breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(limits = c(0, 3), breaks = scales::pretty_breaks(n = 6))+
  custom_theme2 +
  lapply(1:nrow(gaps), function(i) {
    annotation_custom(
      grob = grid::rectGrob(gp = grid::gpar(col = "black", fill = NA, lwd = 0.3)),
      xmin = gaps$start[i], xmax = gaps$end[i], ymin = 3.1, ymax = 3.3)}) +
  coord_cartesian(clip = "off") +#do this to ensure the gaps are plotted although they are outside of the coordinate system
  geom_vline(xintercept = 451.07, color = "black", linetype = "dashed", size = 0.2) +
  geom_vline(xintercept = 455.79, color = "black", linetype = "dashed", size = 0.2)
p_top

##### coliniarity plot Weining_H1 x Lo7 #####
comparison <- read.table("log_chr_comp_Scer_Lo7_V3_chr7R__x__Scer_Weining_v3_H1_chr7R.tab", skip=1)
comparison$Lo7 <- (comparison$V2)/1e+06
comparison$Weining_H1 <- (comparison$V3)/1e+06
comparison$Identity <- as.numeric(str_extract(comparison$V5, "\\d+.\\d+|\\d+"))


p_center <- ggplot(data = comparison, aes(x = Lo7, y = Weining_H1, color = Identity)) +
  geom_point(size=0.1) +
  custom_theme2 +
  scale_x_continuous(limits = c(430, 480), position = "bottom", breaks = scales::pretty_breaks(n = 10)) +
  scale_y_reverse(limits = c(485, 435), position = "left", breaks = scales::pretty_breaks(n = 10)) +
  scale_colour_gradient(low = "gray80", high = "black", space = "Lab", na.value = "grey50", guide = "colourbar", aesthetics = "colour") +
  geom_vline(xintercept = 451.07, color = "black", linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = 455.79, color = "black", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 457.01, color = "black", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 462.40, color = "black", linetype = "dashed", size = 0.5)
p_center

(p_top | p_right) / #add p_right as plot_spacer to ensure that axis line up
  (p_center | p_right) +
  plot_layout(widths = c(4, 4), heights = c(4, 4), guides = 'collect')

##### ChIP seq CENH3 Lo7 #####
bedgraph_Lo7 <- read.table("Lo7_rep_input_100kb_mapping_parameters.log2ratio.bedgraph")
bedgraph_Lo7_chr7 <- bedgraph_Lo7[bedgraph_Lo7$V1 == "chr7R",]

p_CENH3_Lo7 <- ggplot(data = bedgraph_Lo7_chr7, aes(x=V3/1000000, y=V4)) +
  geom_line(color = "black",  linewidth = 0.3) +
  custom_theme2 +
  ylab("log2 CENH3/Input") +
  xlab("") +
  coord_cartesian(xlim = c(430, 480),  ylim = c(-2, 6)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))

p_CENH3_Lo7

(p_CENH3_Lo7 | p_CENH3_Lo7) / 
  (p_top | p_right) / 
  (p_center | p_right) + 
  plot_layout(widths = c(4, 4), heights = c(2, 4, 4), guides = 'collect')

##### ChIP seq CENH3 Weining H1 #####
bedgraph_Weining_H1 <- read.table("WeiningV3.H1_log2ratio.bedgraph")
bedgraph_Weining_H1_chr7 <- bedgraph_Weining_H1[bedgraph_Weining_H1$V1 == "chr7R",]

p_CENH3_Weining_H1 <- ggplot(data = bedgraph_Weining_H1_chr7, aes(x = V3/1000000, y = V4)) +
  geom_line(color = "black", linewidth = 0.3) +
  custom_theme2 +
  ylab("log2 CENH3/Input") +
  xlab("") +
  coord_flip() +
  scale_x_reverse(limits = c(485, 435),
                  breaks = scales::pretty_breaks(n = 10))
p_CENH3_Weining_H1

(p_CENH3_Lo7 | p_CENH3_Lo7 | p_CENH3_Weining_H1) / #add in p_CENH3_Lo7 and p_CENH3_Weining_H1 as spacers to retain exact scales
  (p_top | p_right | p_CENH3_Weining_H1) / #add in p_right and p_CENH3_Weining_H1 as spacers to retain exact scales
  (p_center | p_right | p_CENH3_Weining_H1) + 
  plot_layout(widths = c(4, 4, 4), heights = c(4, 4, 4), guides = 'collect')

#### chromosome 7: Weining H2 x Lo7 ####
##### Preparation Weining #####
setwd(DATA_DIR)
Weining_H2_Cereba <- read.table("age_distribution_LTR_Cereba_H2")
Weining_H2_Cereba$family <- str_extract(Weining_H2_Cereba$V1, "\\S+(?<=_Cereba)")

Weining_H2_Quinta <- read.table("age_distribution_LTR_Quinta_H2")
Weining_H2_Quinta$family <- str_extract(Weining_H2_Quinta$V1, "\\S+(?<=_Quinta)")

Weining_H2_Abia <- read.table("age_distribution_LTR_Abia_H2")
Weining_H2_Abia$family <- str_extract(Weining_H2_Abia$V1, "\\S+(?<=_Abia)")

Weining_H2_Abigail <- read.table("age_distribution_LTR_Abigail_H2")
Weining_H2_Abigail$family <- str_extract(Weining_H2_Abigail$V1, "\\S+(?<=_Abigail)")

all_Weining_H2 <- rbind(Weining_H2_Abia, Weining_H2_Abigail, Weining_H2_Cereba, Weining_H2_Quinta)

all_Weining_H2$TEnum <- str_extract(all_Weining_H2$V1, "\\d\\S-\\d+")
all_Weining_H2$chr <- str_extract(all_Weining_H2$TEnum, "\\d\\S")
all_Weining_H2$TEstart <- as.numeric(str_extract(all_Weining_H2$TEnum, "(?<=-)\\d+"))*1000 #time 1000 to get the start position from the ID 
all_Weining_H2$genome <- as.factor(str_extract(all_Weining_H2$chr, "(?<=\\d)\\S"))

df_sorted_H2 <- all_Weining_H2 %>% arrange(genome)
df_sorted_H2$chr <- factor(df_sorted_H2$chr, levels = rev(unique(df_sorted_H2$chr)))

chr7R_Weining_H2 <- df_sorted_H2[df_sorted_H2$chr == "7R",]
str(chr7R_Weining_H2)
chr7R_Weining_H2$V6 <- as.numeric(chr7R_Weining_H2$V6)

##### gaps #####
gaps <- read.table("Scer_Weining_v3_H2_chr7R_gaps", sep = ";")
gaps$start <- round(gaps$V1/1000000, 2)
gaps$end <- round(gaps$V2/1000000, 2)
gaps <- gaps[gaps$start>420 & gaps$end<490,]

p_right <- ggplot(data = chr7R_Weining_H2, aes(y = TEstart/1000000, x = V6, color = family)) +
  geom_point(shape=20, size=0.8)  +
  xlab("Age [myr]") +
  ylab("Position [Mb] H1") +
  labs(color = "Family") +
  scale_color_manual(values=c( "#e66101", "#fdb863","#5e3c99", "#b2abd2")) + 
  scale_x_continuous(limits = c(0,3), breaks = scales::pretty_breaks(n = 10)) +
  custom_theme2 +
  scale_y_reverse(limits = c(500, 400), breaks = scales::pretty_breaks(n = 10)) +
  geom_rect(data = gaps, aes(xmin = 2.8, xmax = 3, ymin = start, ymax = end), inherit.aes = FALSE, fill = NA, color = "black", size = 0.3)
p_right

##### Preparation Lo7 #####
Lo7_Cereba <- read.table("age_distribution_RLG_ScerLo7_Cereba")
Lo7_Cereba$family <- str_extract(Lo7_Cereba$V1, "\\S+(?<=_Cereba)")

Lo7_Quinta <- read.table("age_distribution_RLG_ScerLo7_Quinta")
Lo7_Quinta$family <- str_extract(Lo7_Quinta$V1, "\\S+(?<=_Quinta)")

Lo7_Abia <- read.table("age_distribution_RLG_ScerLo7_Abia")
Lo7_Abia$family <- str_extract(Lo7_Abia$V1, "\\S+(?<=_Abia)")

Lo7_Abigail <- read.table("age_distribution_RLG_ScerLo7_Abigail")
Lo7_Abigail$family <- str_extract(Lo7_Abigail$V1, "\\S+(?<=_Abigail)")

all_Lo7 <- rbind(Lo7_Abia, Lo7_Abigail, Lo7_Cereba, Lo7_Quinta)

all_Lo7$TEnum <- str_extract(all_Lo7$V1, "\\d\\S-\\d+")
all_Lo7$chr <- str_extract(all_Lo7$TEnum, "\\d\\S")
all_Lo7$TEstart <- as.numeric(str_extract(all_Lo7$TEnum, "(?<=-)\\d+"))*1000 #time 1000 to get the start position from the ID 
all_Lo7$genome <- as.factor(str_extract(all_Lo7$chr, "(?<=\\d)\\S"))

df_sorted_H2 <- all_Lo7 %>% arrange(genome)
df_sorted_H2$chr <- factor(df_sorted_H2$chr, levels = rev(unique(df_sorted_H2$chr)))

chr7R_Lo7 <- df_sorted_H2[df_sorted_H2$chr == "7R",]
str(chr7R_Lo7)
chr7R_Lo7$V6 <- as.numeric(chr7R_Lo7$V6)

##### gaps #####
gaps <- read.table("Scer_Lo7_V3_chr7R_gaps", sep = ";")
gaps$start <- round(gaps$V1/1000000, 2)
gaps$end <- round(gaps$V2/1000000, 2)
gaps <- gaps[gaps$start>425 & gaps$end<490,]

p_top <- ggplot(data = chr7R_Lo7, aes(x = TEstart/1000000, y = V6, color = family)) +
  geom_point(shape=20, size=0.8)  +
  ylab("Age [myr]") +
  xlab("Position [Mb] Lo7") +
  labs(color = "Family") +
  scale_color_manual(values=c("#e66101", "#fdb863", "#5e3c99", "#b2abd2")) + 
  scale_x_continuous(limits = c(400, 500), breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(limits = c(0, 3), breaks = scales::pretty_breaks(n = 6))+
  custom_theme2 +
  lapply(1:nrow(gaps), function(i) {
    annotation_custom(
      grob = grid::rectGrob(gp = grid::gpar(col = "black", fill = NA, lwd = 0.3)),
      xmin = gaps$start[i], xmax = gaps$end[i], ymin = 3.1, ymax = 3.3)}) +
  coord_cartesian(clip = "off") #do this to ensure the gaps are plotted although they are outside of the coordinate system
p_top

#### SUPPLEMENTARY S28 PANEL F ####
##### coliniarity plot Weining_H2 x Lo7 #####
comparison <- read.table("log_chr_comp_Scer_Lo7_V3_chr7R__x__Scer_Weining_v3_H2_chr7R.tab", skip=1)
comparison$Lo7 <- (comparison$V2)/1e+06
comparison$Weining_H2 <- (comparison$V3)/1e+06
comparison$Identity <- as.numeric(str_extract(comparison$V5, "\\d+.\\d+|\\d+"))

custom_theme2_largerTEXT <- theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 22),
    axis.text = element_text(color = "black"),
    axis.line = element_line(color = "black"),
    strip.text = element_text(size = 22),
    panel.grid.major = element_line(color = "lightgray"),  # add back grid
    panel.grid.minor = element_line(color = "lightgray"),   # optional
    panel.background = element_rect(fill = NA, color = NA),
    plot.background = element_rect(fill = NA, color = NA)
  )


ggplot(data = comparison, aes(x = Lo7, y = Weining_H2, color = Identity)) +
  geom_point(size=0.1) +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(size = 22),
    axis.text = element_text(color = "black"),
    axis.line = element_line(color = "black"),
    strip.text = element_text(size = 22),
    panel.grid.major = element_line(color = "lightgray"),  # add back grid
    panel.grid.minor = element_line(color = "lightgray"),   # optional
    panel.background = element_rect(fill = NA, color = NA),
    plot.background = element_rect(fill = NA, color = NA)
  )+
  scale_x_continuous(limits = c(430, 480), position = "bottom", breaks = scales::pretty_breaks(n = 10)) +
  scale_y_reverse(limits = c(485, 435), position = "left", breaks = scales::pretty_breaks(n = 10)) +
  scale_colour_gradient(low = "gray80", high = "black", space = "Lab", na.value = "grey50", guide = "colourbar", aesthetics = "colour")
