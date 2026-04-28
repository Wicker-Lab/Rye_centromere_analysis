#### SUPPLEMENTARY FIGURE S23 ####
### This script contains the code for the comparison of EDTA annotations ###
library(ggplot2)
library(stringr)
library(patchwork)
library(dplyr)

custom_theme <- theme_classic(base_size = 10) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(color = "black"),
    axis.line = element_line(color = "black"),
    strip.text = element_text(size = 12),
    panel.grid.major = element_line(color = "gray"),  # add back grid
    panel.grid.minor = element_line(color = "white"),   # optional
    panel.background = element_rect(fill = NA, color = NA),
    plot.background = element_rect(fill = NA, color = NA)
  )

OUTPUT_DIR <- ""

#### data ####
EDTA_Lo7 <- read.table("Scer_Lo7_V3_centromere_EDTA_summary_for_R", header = T)
EDTA_Lo7$chr_nr <- str_extract(EDTA_Lo7$chr, "(?<=chr)\\S+(?=R)")
EDTA_Lo7$chromosome <- str_extract(EDTA_Lo7$chr, "(?<=chr)\\S+")
EDTA_Lo7$accession <- str_extract(EDTA_Lo7$genome, "(?<=Scer_)(.*)")


EDTA_TA299 <- read.table("Tmon_TA299_centromere_EDTA_summary_for_R", header = T)
EDTA_TA299$chr_nr <- str_extract(EDTA_TA299$chr, "(?<=chr)\\S+(?=A)")
EDTA_TA299$chromosome <- str_extract(EDTA_TA299$chr, "(?<=chr)\\S+")
EDTA_TA299$accession <- str_extract(EDTA_TA299$genome, "(?<=Tmon_)(.*)")

  
#### plot by chromosome ####
EDTA_both <- rbind(EDTA_Lo7,EDTA_TA299)
ggplot(EDTA_both, aes(fill=family, y=percent, x=chr_nr)) + 
  geom_bar(position="stack", stat="identity") +
  custom_theme +
  ylab("EDTA results [%]") +
  xlab("Chromosome") +
  labs(fill = "Family") +
  scale_fill_manual(values=c("#e66101","#fdb863","#5e3c99", "#b2abd2")) +
  facet_wrap(~accession)

p1 <- ggplot(EDTA_Lo7, aes(fill=family, y=percent, x=chromosome)) + 
  geom_bar(position="stack", stat="identity") +
  custom_theme +
  ylab("EDTA results [%]") +
  xlab(expression(paste(italic("S. cereale"), " Lo7v3"))) +
  labs(fill = "Family") +
  scale_fill_manual(values=c("#e66101","#fdb863","#5e3c99", "#b2abd2")) +
  scale_y_continuous(limits = c(0, 99), breaks = scales::pretty_breaks(n = 20))

p2 <- ggplot(EDTA_TA299, aes(fill=family, y=percent, x=chromosome)) + 
  geom_bar(position="stack", stat="identity") +
  custom_theme +
  ylab("") +
  xlab(expression(paste(italic("T. monococcum"), " TA299"))) +
  labs(fill = "Family") +
  scale_fill_manual(values=c("#e66101","#fdb863","#5e3c99", "#b2abd2")) +
  scale_y_continuous(limits = c(0, 99), breaks = scales::pretty_breaks(n = 20), position = "right") +
  theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank()   
  )

(p1 | p2) + plot_layout(widths = c(2, 2), heights = c(4, 4), guides = 'collect')

#### with CS-IAAS genome ####
##### data #####
EDTA_CS <- read.table("EDTA_summary_CS-IAAS", header = T)
EDTA_CS$chr_nr <- str_extract(EDTA_CS$chr, "(?<=chr)\\S+(?=\\S)")
EDTA_CS$chromosome <- str_extract(EDTA_CS$chr, "(?<=chr)\\S+")
EDTA_CS$genome <- as.factor(str_extract(EDTA_CS$chromosome, "(?<=\\d)\\S+"))

df_sorted <- EDTA_CS %>% arrange(genome)
df_sorted$chromosome <- factor(df_sorted$chromosome, levels = unique(df_sorted$chromosome))

ggplot(df_sorted, aes(fill=family, y=percent, x=chromosome)) + 
  geom_bar(position="stack", stat="identity") +
  custom_theme +
  ylab("EDTA results [%]") +
  xlab("Chromosome") +
  labs(fill = "Family") +
  scale_fill_manual(values=c("#e66101","#fdb863","#5e3c99", "#b2abd2")) +
  scale_y_continuous(limits = c(0,99), breaks =  scales::pretty_breaks(n = 20))

p3 <- ggplot(EDTA_CS[EDTA_CS$genome == "A",], aes(fill=family, y=percent, x=chromosome)) + 
  geom_bar(position="stack", stat="identity") +
  custom_theme +
  ylab("EDTA results [%]") +
  xlab("") +
  labs(fill = "Family") +
  scale_fill_manual(values=c("#e66101","#fdb863","#5e3c99", "#b2abd2")) +
  scale_y_continuous(limits = c(0, 99), breaks = scales::pretty_breaks(n = 20))

p4 <- ggplot(EDTA_CS[EDTA_CS$genome == "B",], aes(fill=family, y=percent, x=chromosome)) + 
  geom_bar(position="stack", stat="identity") +
  custom_theme +
  ylab("") +
  xlab(expression(paste(italic("T. aestivum"), " CS IAAS"))) +
  labs(fill = "Family") +
  scale_fill_manual(values=c("#e66101","#fdb863","#5e3c99", "#b2abd2")) +
  scale_y_continuous(limits = c(0, 99), breaks = scales::pretty_breaks(n = 20), position = "right") +
  theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank()
  )

p5 <- ggplot(EDTA_CS[EDTA_CS$genome == "D",], aes(fill=family, y=percent, x=chromosome)) + 
  geom_bar(position="stack", stat="identity") +
  custom_theme +
  ylab("") +
  xlab("") +
  labs(fill = "Family") +
  scale_fill_manual(values=c("#e66101","#fdb863","#5e3c99", "#b2abd2")) +
  scale_y_continuous(limits = c(0, 99), breaks = scales::pretty_breaks(n = 20), position = "right") +
  theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank()   
  )

(p3|p4|p5) + plot_layout(widths = c(2,2,2), heights = c(4,4,4), guides = 'collect')& theme(legend.position = "top")

(p1 | p2 | plot_spacer()) /
  (p3 | p4 | p5 ) + plot_layout(widths = c(2,2,2), heights = c(6,6), guides = 'collect')

