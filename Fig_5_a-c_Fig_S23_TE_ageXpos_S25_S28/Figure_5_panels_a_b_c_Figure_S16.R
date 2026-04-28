#### FIGURE 5 A-C and SUPPLEMENTARY FIGURE S16 ####
### the variables in the different sections of this script are NOT unique
### => make sure to run one chromosome at a time and save output!

#### FIGURE 5 PANEL A ####
library(magrittr)
library(tidyverse)
library(ggridges)
library(ggplot2)
library(dplyr)
library(stringr)
library(patchwork)

##### Preparation #####

OUTPUT_DIR <- ""
CHR_SIZE_DIR <- ""

# adjust per group
GROUP_DIR <- ""
setwd(GROUP_DIR)

name <- "Abia Abigail Cereba Quinta"
species <- expression(paste(italic("Sercale cereale"), " Lo7v3"))
#species <- expression(paste(italic("Sercale cereale"), " Weining"))
family <- "RLG_ScerLo7"
#family <- "RLG_Weining"
family_name <- expression(bold("Abia Abigail Cereba Quinta"))

ScerLo7_Cereba <- read.table("age_distribution_RLG_ScerLo7_Cereba", skip=1)
ScerLo7_Cereba$family <- str_extract(ScerLo7_Cereba$V1, "\\S+(?<=_Cereba)")

ScerLo7_Quinta <- read.table("age_distribution_RLG_ScerLo7_Quinta", skip=1)
ScerLo7_Quinta$family <- str_extract(ScerLo7_Quinta$V1, "\\S+(?<=_Quinta)")

ScerLo7_Abia <- read.table("age_distribution_RLG_ScerLo7_Abia", skip=1)
ScerLo7_Abia$family <- str_extract(ScerLo7_Abia$V1, "\\S+(?<=_Abia)")

ScerLo7_Abigail <- read.table("age_distribution_RLG_ScerLo7_Abigail", skip=1)
ScerLo7_Abigail$family <- str_extract(ScerLo7_Abigail$V1, "\\S+(?<=_Abigail)")

all_Lo7 <- rbind(ScerLo7_Abia, ScerLo7_Abigail, ScerLo7_Cereba, ScerLo7_Quinta)

all_Lo7$TEnum <- str_extract(all_Lo7$V1, "\\d\\S-\\d+")
all_Lo7$chr <- str_extract(all_Lo7$TEnum, "\\d\\S")
all_Lo7$TEstart <- as.numeric(str_extract(all_Lo7$TEnum, "(?<=-)\\d+"))*1000 #time 1000 to get the start position from the ID 
all_Lo7$genome <- as.factor(str_extract(all_Lo7$chr, "(?<=\\d)\\S"))

df_sorted <- all_Lo7 %>% arrange(genome)
df_sorted$chr <- factor(df_sorted$chr, levels = rev(unique(df_sorted$chr)))

##### import sizes of the chromosomes #####
setwd(CHR_SIZE_DIR)
sizeList <- read.table("Scer_Lo7_V3_genome_size_list", header = TRUE)
sizeList <- sizeList %>% arrange(Name)
sizeList$chr <- str_extract(sizeList$Name, "(?<=chr)\\d\\S")
sizeList$chromosome <- str_extract(sizeList$chr, "\\d(?<=\\S)")
sizeList$genome <- as.factor(str_extract(sizeList$chr, "(?<=\\d)\\S"))
sizeList <- sizeList %>% arrange(genome)

sizeList <- mutate(sizeList, chr_nr=1:length(sizeList$Name))

sizeList <- sizeList %>% arrange(desc(chr_nr))
sizeList <- mutate(sizeList, chr_nr_sorted=1:length(sizeList$Name))

all_Lo7$TEnum <- str_extract(all_Lo7$V1, "\\d\\S-\\d+")
all_Lo7$chr <- str_extract(all_Lo7$TEnum, "\\d\\S")
all_Lo7$TEstart <- as.numeric(str_extract(all_Lo7$TEnum, "(?<=-)\\d+"))*1000 #time 1000 to get the start position from the ID 
all_Lo7$genome <- as.factor(str_extract(all_Lo7$chr, "(?<=\\d)\\S"))

df_sorted <- all_Lo7 %>% arrange(genome)
df_sorted$chr <- factor(df_sorted$chr, levels = rev(unique(df_sorted$chr)))
df_sorted$chromosome_full <- str_extract(df_sorted$V1, "chr\\d+[A-Za-z]")
df_sorted <- merge(df_sorted, sizeList, by = "chr")
df_sorted <- df_sorted[df_sorted$V6<=3.0,]


##### import centromere position #####
Lo7_centromere_pos <- read.table("centromere_positions_Lo7", header = T)
Lo7_centromere_pos$start_full <- Lo7_centromere_pos$start*1e+06
Lo7_centromere_pos$end_full <- Lo7_centromere_pos$end*1e+06
Lo7_centromere_pos$chr <- str_extract(Lo7_centromere_pos$Lo7, "(?<=chr)\\d\\S")
df_sorted_centromere_pos <-merge(df_sorted, Lo7_centromere_pos, by = "chr", all.x = T)

ggplot(data = df_sorted_centromere_pos, aes(x = TEstart/1000000, y = chr, color = family)) +
  geom_point(shape=108, size=4.5, alpha= 0.4) +
  labs(color = "Families") +
  xlab("Position [Mb]") +
  ylab("") +
  scale_color_manual(values=c("#e66101","#fdb863","#5e3c99", "#b2abd2")) +
  scale_x_continuous(limits = c(0, 1200), breaks = scales::pretty_breaks(n = 10)) +
  geom_rect(aes(xmin=0, xmax=Length/1000000, ymin=chr_nr_sorted-0.32, ymax=chr_nr_sorted+0.38), color ="black", fill = NA, linewidth=0.1, inherit.aes = F) +
  geom_rect(aes(xmin=start, xmax=end, ymin=chr_nr_sorted+0.44, ymax=chr_nr_sorted+0.48), color = "black", fill = "black", linewidth=0.1, inherit.aes = F) +
  #  geom_point(aes(x = center, y = chr_nr_sorted+0.44), shape = 6, size=1, color = "black") +
  ggtitle(species, family_name) +
  theme_classic(base_size = 12) +
  facet_grid(family~.)

ggplot(data = df_sorted_centromere_pos, aes(x = TEstart, y = chr, color = family)) +
  geom_point(shape=108, size=4.5, alpha= 0.4) +
  labs(color = "Families") +
  xlab("Position [Mb]") +
  ylab("") +
  scale_color_manual(values=c("#e66101","#fdb863","#5e3c99", "#b2abd2")) +
  scale_x_continuous(limits = c(0, 12e08), breaks = scales::pretty_breaks(n = 10)) +
  geom_rect(aes(xmin=0, xmax=Length, ymin=chr_nr_sorted-0.32, ymax=chr_nr_sorted+0.38), color ="black", fill = NA, linewidth=0.1, inherit.aes = F) +
  geom_rect(aes(xmin=start_full, xmax=end_full, ymin=chr_nr_sorted+0.44, ymax=chr_nr_sorted+0.48), color = "black", fill = "black", linewidth=0.1, inherit.aes = F) +
  #  geom_point(aes(x = center, y = chr_nr_sorted+0.44), shape = 6, size=1, color = "black") +
  ggtitle(species, family_name) +
  theme_minimal() +
  facet_grid(family~.)

#### FIGURE 5 PANEL B ####

##### DEFINE THE CHROMOSOME ----------------------------------------------------------------------------------------------------------
chromosome_of_choice <- as.character("chr4R")
start_cent <- 350
end_cent <- 430
#-------------------------------------------------------------------------------------------------------------------------------------

##### import ChIP-seq data #####
bedgraph <- read.table("Lo7_rep_input_100kb_mapping_parameters.log2ratio.bedgraph")

# subset chromosomes
chr_choice <- df_sorted[df_sorted$chromosome_full == chromosome_of_choice,]
chr_choice <- chr_choice[chr_choice$TEstart>((start_cent*1000000)-4100000) & chr_choice$TEstart<((end_cent*1000000)+4100000),]#filter TEs in displayed region
chr_choice$V6 <- as.numeric(chr_choice$V6)
# subset bedgraph
bedgraph_choice<- bedgraph[bedgraph$V1 == chromosome_of_choice,]

# read in gap data
gaps <- read.table(paste0(CHR_SIZE_DIR, "Scer_Lo7_V3_", chromosome_of_choice,"_gaps"), sep = ";")
gaps$start <- round(gaps$V1/1000000, 0)
gaps$end <- round(gaps$V2/1000000, 0)
gaps <- gaps[gaps$start>start_cent-5 & gaps$end<end_cent+5,]


p1 <- ggplot(data = chr_choice, aes(x = TEstart/1000000, y = V6, color = family)) +
  geom_point(shape=20, size=0.6)  +
  ylab("Age [myr]") +
  xlab("") +
  labs(color = "Family", tag = chromosome_of_choice)+
  scale_color_manual(values=c("#e66101","#fdb863","#5e3c99", "#b2abd2")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(limits = c(0, 3)) +
  theme_classic(base_size = 12) +
  lapply(1:nrow(gaps), function(i) {
    annotation_custom(
      grob = grid::rectGrob(gp = grid::gpar(col = "black", fill = NA, lwd = 1.2)),
      xmin = gaps$start[i], xmax = gaps$end[i], ymin = 3.2, ymax = 3.4)}) +
  coord_cartesian(xlim = c(start_cent, end_cent), clip = "off") +
  geom_segment(x = start_cent, xend = end_cent, y = 3.3, yend = 3.3, color = "black", linewidth = 0.4)

p1

p2 <- ggplot(data = bedgraph_choice, aes(x=V3/1000000, y=V4)) +
  geom_line(color = "black",  linewidth = 0.4) +
  theme_classic() +
  ylab("log2 CENH3/Input") +
  xlab("Position [Mb]") +
  coord_cartesian(xlim = c(start_cent, end_cent),  ylim = c(-2, 6)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))

p2

layout <- "
A
A
A
A
A
B
B
"

p1+p2+plot_layout(design = layout)


#### FIGURE 5 PANEL C ####

##### DEFINE THE CHROMOSOME ----------------------------------------------------------------------------------------------------------
chromosome_of_choice <- as.character("chr7R")
start_cent <- 400
end_cent <- 480
#-------------------------------------------------------------------------------------------------------------------------------------

##### import ChIP-seq data #####
bedgraph <- read.table("Lo7_rep_input_100kb_mapping_parameters.log2ratio.bedgraph")

# subset chromosomes
chr_choice <- df_sorted[df_sorted$chromosome_full == chromosome_of_choice,]
chr_choice <- chr_choice[chr_choice$TEstart>((start_cent*1000000)-4100000) & chr_choice$TEstart<((end_cent*1000000)+4100000),]#filter TEs in displayed region

# subset bedgraph
bedgraph_choice<- bedgraph[bedgraph$V1 == chromosome_of_choice,]

# read in gap data
gaps <- read.table(paste0(CHR_SIZE_DIR, "Scer_Lo7_V3_", chromosome_of_choice,"_gaps"), sep = ";")
gaps$start <- round(gaps$V1/1000000, 0)
gaps$end <- round(gaps$V2/1000000, 0)
gaps <- gaps[gaps$start>start_cent-5 & gaps$end<end_cent+5,]


p1 <- ggplot(data = chr_choice, aes(x = TEstart/1000000, y = V6, color = family)) +
  geom_point(shape=20, size=0.6)  +
  ylab("Age [myr]") +
  xlab("Position [Mb]") +
  labs(color = "Family")+
  scale_color_manual(values=c("#e66101","#fdb863","#5e3c99", "#b2abd2")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_classic(base_size = 12) +
  lapply(1:nrow(gaps), function(i) {
    annotation_custom(
      grob = grid::rectGrob(gp = grid::gpar(col = "black", fill = NA, lwd = 1.2)),
      xmin = gaps$start[i], xmax = gaps$end[i], ymin = 3.2, ymax = 3.3)}) +
  coord_cartesian(xlim = c(start_cent, end_cent), clip = "off") +
  geom_segment(x = start_cent, xend = end_cent, y = 3.25, yend = 3.25, color = "black", linewidth = 0.4) +
  geom_segment(x = 458, y = -0.15, xend = 458, yend = 3, linetype = "solid", color = "gray75", size = 0.06) + #A1
  annotate("text", x = 458, y = 3.05, label = "A1", color = "gray65", size = 7/.pt) +
  geom_segment(x = 466.5, y = -0.15, xend = 466.5, yend = 3, linetype = "solid", color = "gray75", size = 0.06) + #A2
  annotate("text", x = 466.5, y = 3.05, label = "A2", color = "gray65", size = 7/.pt) +
  geom_segment(x = 416.5, y = -0.15, xend = 416.5, yend = 2.9, linetype = "solid", color = "gray75", size = 0.06) + #B1
  annotate("text", x = 416.5, y = 2.95, label = "B1", color = "gray65", size = 7/.pt) +
  geom_segment(x = 433.5, y = -0.15, xend = 433.5, yend = 2.9, linetype = "solid", color = "gray75", size = 0.06) + #B2
  annotate("text", x = 433.5, y = 2.95, label = "B2", color = "gray65", size = 7/.pt) +
  geom_segment(x = 451.07, y = -0.15, xend = 451.07, yend = 2.75, linetype = "solid", color = "gray75", size = 0.06) + #C1
  annotate("text", x = 451.07, y = 2.8, label = "C1", color = "gray65", size = 7/.pt) +
  geom_segment(x = 455.79, y = -0.15, xend = 455.79, yend = 2.75, linetype = "solid", color = "gray75", size = 0.06) + #C2
  annotate("text", x = 455.79, y = 2.8, label = "C2", color = "gray65", size = 7/.pt) +
  theme(legend.position = "none")
p1

p2 <- ggplot(data = bedgraph_choice, aes(x=V3/1000000, y=V4)) +
  geom_line(color = "black",  linewidth = 0.4) +
  theme_classic() +
  ylab("log2 CENH3/Input") +
  xlab("") +
  coord_cartesian(xlim = c(start_cent, end_cent),  ylim = c(-2, 6)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))

p2

layout <- "
B
A
A
A
A
"

p1+p2+plot_layout(design = layout)

#### SUPPLEMENTARY FIGURE 16 ####

##### DEFINE THE CHROMOSOME ----------------------------------------------------------------------------------------------------------
#choose display windows based on centromere positions
chromosome_of_choice <- as.character("chr1R")
start_cent <-  260
end_cent <-  340

chromosome_of_choice <- as.character("chr2R")
start_cent <- 410
end_cent <- 490

chromosome_of_choice <- as.character("chr3R")
start_cent <- 480
end_cent <- 560

chromosome_of_choice <- as.character("chr4R")
start_cent <- 350
end_cent <- 430

chromosome_of_choice <- as.character("chr5R")
start_cent <- 250
end_cent <- 330

chromosome_of_choice <- as.character("chr6R")
start_cent <- 270
end_cent <- 350

chromosome_of_choice <- as.character("chr7R")
start_cent <- 400
end_cent <- 480

#-------------------------------------------------------------------------------------------------------------------------------------

##### import ChIP-seq data #####
bedgraph <- read.table("Lo7_rep_input_100kb_mapping_parameters.log2ratio.bedgraph")

# subset chromosomes
chr_choice <- df_sorted[df_sorted$chromosome_full == chromosome_of_choice,]
chr_choice <- chr_choice[chr_choice$TEstart>((start_cent*1000000)-4100000) & chr_choice$TEstart<((end_cent*1000000)+4100000),]#filter TEs in displayed region

# subset bedgraph
bedgraph_choice <- bedgraph[bedgraph$V1 == chromosome_of_choice,]

# read in gap data
gaps <- read.table(paste0(CHR_SIZE_DIR, "Scer_Lo7_V3_", chromosome_of_choice,"_gaps"), sep = ";")
gaps$start <- round(gaps$V1/1000000, 0)
gaps$end <- round(gaps$V2/1000000, 0)
gaps <- gaps[gaps$start>start_cent-5 & gaps$end<end_cent+5,]

p1 <- ggplot(data = chr_choice, aes(x = TEstart/1000000, y = V6, color = family)) +
  geom_point(shape=20, size=0.6)  +
  ylab("Age [myr]") +
  xlab("") +
  labs(color = "Family", tag = chromosome_of_choice)+
  scale_color_manual(values=c("#e66101","#fdb863","#5e3c99", "#b2abd2")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_classic(base_size = 12) +
  lapply(1:nrow(gaps), function(i) {annotation_custom(grob = grid::rectGrob(gp = grid::gpar(col = "black", fill = NA, lwd = 1.2)), xmin = gaps$start[i], xmax = gaps$end[i], ymin = 3.15, ymax = 3.35)}) + #use this line for everything exept chr6R, because there are no gaps in the centromere to annotate 
  coord_cartesian(xlim = c(start_cent, end_cent), clip = "off") +
  geom_segment(x = start_cent, xend = end_cent, y = 3.25, yend = 3.25, color = "black", linewidth = 0.4)
p1

p2 <- ggplot(data = bedgraph_choice, aes(x=V3/1000000, y=V4)) +
  geom_line(color = "black",  linewidth = 0.4) +
  theme_classic() +
  ylab("log2 CENH3/Input") +
  xlab("Position [Mb]") +#only for 4R and 7R
  coord_cartesian(xlim = c(start_cent, end_cent),  ylim = c(-2, 6)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))

p2

layout <- "
A
A
A
B
B
"

png(paste0(OUTPUT_DIR, "/", "Sercale_cereale_ChIP_Abia_Abigail_Cereba_Quinta_ageXpos_colorblindfrindly_", chromosome_of_choice, "_v2.png"), width = 7, height = 4, units = "in", res = 600)
plot <- p1+p2+plot_layout(design = layout)
plot
assign(chromosome_of_choice, plot, envir = .GlobalEnv)
dev.off()

(chr1R | chr5R)/
  (chr2R | chr6R)/
  (chr3R | chr7R)/
  (chr4R | guide_area()) + plot_layout(widths = c(1,1,1,1), heights = c(2.5,2.5,2.5,2.5), guides = 'collect') &
  theme(legend.position = "bottom")
