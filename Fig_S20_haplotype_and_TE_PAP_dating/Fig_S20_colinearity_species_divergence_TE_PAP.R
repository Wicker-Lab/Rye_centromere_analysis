

library("ggplot2")
library("dplyr")
library("tidyverse")
library("patchwork")
library("heatmap3")
library("heatmaply")

setwd("~/data/dir_rye_genome/dir_Scer_genome_2024/dir_R_scripts")

getwd()



# plot for haplotype divergence times overall (used for figure S25c)  ---------------- 

colors <- c("#0088FF","#FF0000","#FF5500","#FFBB00","#000000","#00FF00")

infile <- "chr_comp_date_log_all_win400"
df <-read.table(infile, sep='\t', header=T)
head(df)

# calc medians and means to show in plots
medians <- aggregate(age ~ set_ID, data = df, FUN = median)
medians <- as.data.frame(medians)
medians

means <- aggregate(age ~ set_ID, data = df, FUN = mean)
means <- as.data.frame(means)
means

med <- merge(medians,means, by="set_ID")
med$medi <- round(med$age.x,4)
med$avg <- round(med$age.y,4)

med

# make counts for groups 
ns <- aggregate(age ~ set_ID, data = df, FUN = length)
ns <- as.data.frame(ns)
ns

p <- ggplot(df, aes(y = set_ID, x = age)) +
  geom_boxplot() +
  geom_jitter(width=0,height=0.3,alpha=0.4,size=1,aes(color=set_ID))+
  scale_color_manual(values = colors) +
  #geom_violin(alpha=0.5)+
  xlab("Haplotype divergence [MYA]")+
  ylab("Data set")+
  xlim(0,0.2)+
  # add median values
  geom_text(
    data = med,
    aes(x = medi, y = set_ID,label=paste(medi,", ",avg),hjust=0.1,vjust=-5.7),
    inherit.aes = FALSE,
    color = "black",
    size = 3
  ) +
  
  # add counts per group
  geom_text(
    data = ns,
    aes(x = 0.17, y = set_ID,label=paste("n=",age),hjust=0.1,vjust=0),
    inherit.aes = FALSE,
    color = "black",
    size = 3
  ) +
  
  #xlim(0,0.015)+
  theme_classic()
p




# write plot to output file, half page (8 cm) 
outfile <- paste(infile,"_box_H1_x_H2_divergence_overall",sep='')
path <- getwd()
out_png <- paste(outfile,".png",sep='')
wide <- 7
high <- 2.7
ggsave(p, filename=out_png, path = path, width=wide,height=high,dpi=600)










# plot for TE insertion times for individual regions, used for Fig S 27d---------------- 

colors <- c("#0088FF","#0000FF","#AA0000","#FF0000","#FF5500","#FFBB00")

infile <- "chr_comp_TE_pos_age_log_all"
df <-read.table(infile, sep='\t', header=T)
head(df)

df$pap <- paste(df$set_ID,"_",df$set_ID2,sep='')
head(df)

# calc medians and means to show in plots
medians <- aggregate(age ~ pap, data = df, FUN = median)
medians <- as.data.frame(medians)
medians

means <- aggregate(age ~ pap, data = df, FUN = mean)
means <- as.data.frame(means)
means

med <- merge(medians,means, by="pap")
med$medi <- round(med$age.x,4)
med$avg <- round(med$age.y,4)

med

# make counts for groups 
ns <- aggregate(age ~ pap, data = df, FUN = length)
ns <- as.data.frame(ns)
ns

p <- ggplot(df, aes(y = pap, x = age)) +
  geom_boxplot() +
  geom_jitter(width=0,height=0.3,alpha=0.7,size=0.3,aes(color=pap))+
  scale_color_manual(values = colors) +
  #geom_violin(alpha=0.5)+
  xlab("TE insertion age [MYA]")+
  ylab("Data set")+
  xlim(0,0.2)+
  
  # add median values
  geom_text(
    data = med,
    aes(x = age.x, y = pap,label=paste(medi,", ",avg),hjust=0.1,vjust=-3.7),
    inherit.aes = FALSE,
    color = "black",
    size = 3
  ) +
  
  # add counts per group
  geom_text(
    data = ns,
    aes(x = 0.17, y = pap,label=paste("n=",age),hjust=1,vjust=0),
    inherit.aes = FALSE,
    color = "black",
    size = 3
  ) +
  
  #xlim(0,0.015)+
  theme_classic()
p


# write plot to output file, half page (8 cm) 
outfile <- paste(infile,"_box_H1_x_H2_divergence_regions",sep='')
path <- getwd()
out_png <- paste(outfile,".png",sep='')
wide <- 8.1
high <- 3.6
ggsave(p, filename=out_png, path = path, width=wide,height=high,dpi=600)




# make a few stats on young insertions (not used for figure)
infile <- "chr_comp_TE_pos_age_log_all_v2"
df <-read.table(infile, sep='\t', header=T)
head(df)

young_H2 <- df[df$age==0 & df$level >3,]
old_H2 <- df[df$age>0 & df$level >3,]

young_H1 <- df[df$age==0 & df$level <=3,]
old_H1 <- df[df$age>0 & df$level <=3,]

cent_H2_Abia <- df[df$pos1>435000000 & df$pos2 <=448000000 & df$ID=="RLG_Abia",]
cent_H2_Abigail <- df[df$pos1>435000000 & df$pos2 <=448000000 & df$ID=="RLG_Abigail",]


df$pap <- paste(df$set_ID,"_",df$set_ID2,sep='')
head(df)






# plot annotation TE PAP, useed for Fig 27e---------------------------------

infile <- "annotation_TE_PAP_Lo7_x_H2_x_H1"
df <-read.table(infile, sep='\t', header=T)
head(df)

# make counts for groups 
ns <- aggregate(age ~ level, data = df, FUN = length)
ns <- as.data.frame(ns)
ns

colors <- c("#FF0000","#0000FF")

p <- ggplot(df, aes(x = pos1/1000000, y = level)) +
  geom_point(shape='|',size=4,aes(color=fam)) +
  scale_color_manual(values = colors) +
  #geom_violin(alpha=0.5)+
  xlab("Position on Lo7 chr7R [Mb]")+
  ylab("Data set")+
  ylim(1,6)+
  # add counts per group
  geom_text(
    data = ns,
    aes(x = 485, y = level,label=paste("n=",age),hjust=2,vjust=0),
    inherit.aes = FALSE,
    color = "black",
    size = 3
  ) +
  
  
  theme_classic()
p




# write plot to output file 
outfile <- paste(infile,"_annotation_plot",sep='')
path <- getwd()
out_png <- paste(outfile,".png",sep='')
wide <- 11
high <- 5
ggsave(p, filename=out_png, path = path, width=wide,height=high,dpi=600)



# linear plot for age estimates along chromosomes (not used for figure) ----------------
colors <- c("#0088FF","#FF0000","#FF5500","#FFBB00")

infile <- "chr_comp_date_log_all_win400"
df <-read.table(infile, sep='\t', header=T)
head(df)

p <- ggplot(df[df$flag=='centL',], aes(x = pos1,y=age, color=age)) +
  geom_point(size=0.5)+
  theme_classic()
p

