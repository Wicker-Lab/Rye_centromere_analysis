

library("ggplot2")
library("ggfortify") #maybe not needed ?
library("ggrepel")
library("dplyr")
library("tidyverse")
library(sjPlot)
library("patchwork")
library("heatmap3")
library("heatmaply")

setwd("~/data/dir_rye_genome/dir_Scer_genome_2024/dir_R_scripts")

getwd()





# TE feature count ranking (TE expressio levels), used for Fig, S31a and S31b -----------

top <-20

# for wheat, 12 reps, used for Fig. S31b -----------------------
infile <- "feature_counts_compiled_Taes_R_RLX"
df <-read.table(infile, sep='\t', header=T)
head(df)

spec <- "Taes"

# add avg2 colums
df$avg2 <- rowMeans(df[ , 3:14], na.rm = TRUE)
df$median2 <- apply(df[ , 3:14], 1, median, na.rm = TRUE)

# exclude short elements
dfl <- df[df$len>1000,]

# order by avg expression
df_ord <- dfl[order(dfl$median2,decreasing=TRUE),]

top20 <- df_ord[0:top,0:19] # for wheat
df_long <- pivot_longer(top20,names_to="rep",values_to="fc",cols=3:14) # for wheat
head(df_long)

# add color flag for centromeric/other TEs 
df_long$TE_group <- "O"

df_long <- df_long %>%
  mutate(TE_group = ifelse(grepl("Cereba|Abia|Quinta|Abigail", ID), "Centromeric", "Other"))


head(df_long)

df_long$ID <- factor(df_long$ID, levels=rev(unique(df_long$ID)))
colors <- c("Centromeric" = "red", "Other" = "gray40")


p <- ggplot(df_long, aes(y = ID, x = fc, color = TE_group)) +
  geom_boxplot() +  # keep boxes white inside
  scale_color_manual(values = colors) +
  geom_jitter(width=1,height=0.3,alpha=1,size=0.3,color='blue')+
  xlab("Expression level [RPKB]")+
  ylab("TE family")+
  theme_classic()+
  theme(axis.text.y = element_text(hjust = 0, face = "italic"))
p

# write plot to output file
outfile <- paste(spec,"_TE_expression_v2_RLX_FR_top",top,sep='')
path <- getwd()
out_png <- paste(outfile,".png",sep='')
wide <- 7
high <- wide/35*top
ggsave(p, filename=out_png, path = path, width=wide,height=high,dpi=600)





# for rye (Scer), 12 reps, used for Figure S31a --------------------
infile <- "feature_counts_compiled_Scer_FR_RLX"
df <-read.table(infile, sep='\t', header=T)
head(df)

spec <- "Scer"

# add avg2 colums
df$avg2 <- rowMeans(df[ , 3:14], na.rm = TRUE)
df$median2 <- apply(df[ , 3:14], 1, median, na.rm = TRUE)

# exclude short elements
dfl <- df[df$len>1000,]

# order by avg expression
df_ord <- dfl[order(dfl$median2,decreasing=TRUE),]

top20 <- df_ord[0:top,0:19] # for wheat
df_long <- pivot_longer(top20,names_to="rep",values_to="fc",cols=3:14) # for Scer 6 repl



# add color flag for centromeric/other TEs 
df_long$TE_group <- "O"

df_long <- df_long %>%
  mutate(TE_group = ifelse(grepl("Cereba|Abia|Quinta|Abigail", ID), "Centromeric", "Other"))


head(df_long)

df_long$ID <- factor(df_long$ID, levels=rev(unique(df_long$ID)))
colors <- c("Centromeric" = "red", "Other" = "gray40")

p <- ggplot(df_long, aes(y = ID, x = fc, color = TE_group)) +
  geom_boxplot() +  # keep boxes white inside
  scale_color_manual(values = colors) +
  geom_jitter(width=1,height=0.3,alpha=1,size=0.3,color='blue')+
  xlab("Expression level [RPKB]")+
  ylab("TE family")+
  theme_classic()+
  theme(axis.text.y = element_text(hjust = 0, face = "italic"))
p

# write plot to output file
outfile <- paste(spec,"_TE_expression_v2_RLX_FR_top",top,sep='')
path <- getwd()
out_png <- paste(outfile,".png",sep='')
wide <- 7
high <- wide/35*top
ggsave(p, filename=out_png, path = path, width=wide,height=high,dpi=600)








# TE insertion age distribution analysis, used for Fig. S31c------------------
infile <- "all_TE_centromere_pos_vs_age_vs_pos"
df <-read.table(infile, sep='\t', header=T)
head(df)

# make counts for all families in subgenomes
n_labels <- df[df$species != 'HvulH',] %>%
  group_by(fam3) %>%
  summarise(n = n(), .groups = "drop")
n_labels <- as.data.frame(n_labels)
head(n_labels)

# calc medians of insertions age  for all families in subgenomes
medians <- df[df$species != 'HvulH',] %>%
  group_by(fam3) %>%
  summarise(median = median(age, na.rm=TRUE))

medians <- as.data.frame(medians)
medians

# calc minimum of age
age_min <- df[df$species != 'HvulH',] %>%
  group_by(fam3) %>%
  summarise(median = min(age, na.rm=TRUE))

age_min <- as.data.frame(age_min)
age_min

xmax <- 9

p <- ggplot(df[df$species != 'HvulH',], aes(y = fam3, x = age)) +
#p <- ggplot(df, aes(y = fam3, x = age)) +
  geom_jitter(width=0.1,height=0.3,alpha=0.6,size=0.2,aes(color=deltaC/1000))+
  scale_color_gradient2(low = "yellow", mid="red",high = "black",midpoint =25) +
  geom_boxplot(alpha=0.1) +  # keep boxes white inside
  xlim(0,xmax)+
  geom_text(data = n_labels, aes(x = (xmax*0.95),y = fam3,label = paste0("n= ",n)),hjust=0,vjust=0, size = 3  ) +
  geom_text(data = medians, aes(x = (xmax*0.8),y = fam3,label = paste0("med= ",median)),hjust =0, vjust=0, size = 3  ) +
  geom_text(data = age_min, aes(x = (xmax*0.65),y = fam3,label = paste0("min= ",median)),hjust =0, vjust=0, size = 3  ) +
  xlab("Insertion age [Myr]")+
  ylab("TEheatmap3")+
  scale_y_discrete(limits = rev)+
  theme_classic()

p


# write plot to output file
outfile <- "TE_insertion_ages_centromeric"
path <- getwd()
out_png <- paste(outfile,".png",sep='')
wide <- 11
high <- 6
ggsave(p, filename=out_png, path = path, width=wide,height=high,dpi=600)


