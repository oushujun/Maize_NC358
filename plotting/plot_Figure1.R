library(tidyr)
library(ggplot2)
library(wesanderson)
library(scales)
library(reshape2)
library(plyr) 
library(grid)
library(gridExtra)
library(lattice)
library(tibble)
#install.packages("devtools")
#devtools::install_github("thomasp85/patchwork")
library(patchwork)

#plot order
ID = c("21k_20x" ,"21k_30x" ,"21k_40x" ,"21k_50x" ,"21k_60x" ,"21k_75x", "11k_50x" ,"16k_50x")

###########
## BUSCO ##
###########
datadir="/Users/oushujun/Google Drive/study/Maize research/NC358/basic stats"
setwd(datadir)

BUSCO = read.table("NC358_BUSCO.txt",header=T)
str(BUSCO)

#make plot
BUSCO_p = BUSCO %>% 
  gather(variable, value, Complete) %>% 
  ggplot(aes(x = Assembly, y = value, fill = "#5F9EA0")) + 
  geom_bar(stat = "identity",colour="black") +
  scale_fill_manual(values = c("#5F9EA0"),labels = c("BUSCO (%)"),guide = guide_legend(reverse=TRUE)) + theme_classic() +
  labs(x =" ", y = "Complete BUSCO (%)")+
  theme(axis.title.y = element_text(size=18, face="bold"),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(face="bold", size=14, angle=35,vjust=1,hjust=0.95),
        legend.text = element_text(size=18),
        legend.title=element_blank()) +
  theme(strip.text.y = element_text(size = 17)) +
  theme(legend.position="none") +
  coord_cartesian(ylim=c(0 ,100)) +
  scale_x_discrete(limits=ID) + 
  theme(strip.text.x = element_text(size = 17))
BUSCO_p


#################################
## Contig and Bionano conflict ##
#################################

ylim.prim <- c(0, 700)   # in this example, Bionano.conflict
ylim.sec <- c(0, 12000)    # in this example, Contig number

#y-axes conversion
b <- diff(ylim.prim)/diff(ylim.sec)
a <- b*(ylim.prim[1] - ylim.sec[1]) + 20

cols <- c("Bionano"="#5F9EA0", "Contig"="#d8b365")
contig_p = ggplot(BUSCO, aes(x=Assembly)) + 
  geom_col(aes(y = Bionano.conflict, fill = "Bionano"), colour="black") +
  geom_line(aes(y = a + Contig.number*b, group=1, color = "Contig"), size = 1.5) + 
  scale_y_continuous("Bionano conflict", sec.axis = sec_axis(~ (. - a)/b, name = "Contig number")) +
  scale_x_discrete(limits=ID) + 
  labs(x =" ", size=14, face="bold")+
  coord_cartesian(ylim=c(0 ,700)) +
  scale_colour_manual(values=cols) + scale_fill_manual(values=cols) +
  theme_classic() +
  theme(axis.title.y = element_text(size=20, face="bold"),
        axis.text.y = element_text(size=18, angle=90, hjust=0.5), 
        axis.title.y.right = element_text(size=20, face="bold", angle=90),
        axis.text.y.right = element_text(size=18, angle=90, hjust=0.5), 
        axis.text.x = element_text(face="bold", size=20, angle=35,vjust=1,hjust=0.95),
        legend.text = element_text(size=18),
        legend.title=element_blank()) +
  theme(strip.text.y = element_text(size = 17)) +
  theme(legend.position = c(0.15, 1),
        legend.box = "horizontal",
        legend.justification = c("left", "top"))
  
contig_p


#########
## LAI ##
#########
datadir="/Users/oushujun/Google Drive/study/Maize research/NC358/TEs"
setwd(datadir)

LAI = read.table('NC358_all_pseudomolecules_v1.fasta.LAI', header=T, sep = "\t")
str(LAI)

LAI_p = LAI %>%
  ggplot(aes(x = Genome, y = LAI)) + 
  geom_boxplot(colour="black") +
  theme_classic() +
  labs(x =" ", y = "Regional LAI (3Mb)",size=14, face="bold")+
  theme(axis.title.y = element_text(size=20, face="bold"),
        axis.text.y = element_text(size=18), 
        axis.text.x = element_text(face="bold", size=20, angle=35,vjust=1,hjust=0.95),
        legend.text = element_text(size=18),
        legend.title=element_blank()) +
  theme(strip.text.y = element_text(size = 17)) +
  theme(legend.position="none") +
  scale_x_discrete(limits=ID) + 
  theme(strip.text.x = element_text(size = 17))
LAI_p


#############
## RNA-seq ##
#############
datadir="/Users/oushujun/Google Drive/study/Maize research/NC358/RNA"
setwd(datadir)

RNA = read.table('nc358_RNASeq.txt', header=T, sep = "\t")
str(RNA)

RNA_p = ggplot(RNA, aes(x=Assembly, y=uniquely_mapped_percent)) + 
  geom_boxplot(colour="black")+
  geom_jitter(position=position_jitter(0.2)) + 
  theme_classic() +
  labs(x =" ", y = "Uniquely mapped (%)",size=14, face="bold")+
  theme(axis.title.y = element_text(size=20, face="bold"),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(face="bold", size=20, angle=35,vjust=1,hjust=0.95),
        legend.text = element_text(size=18),legend.title=element_blank()) +
  theme(strip.text.y = element_text(size = 17)) +
  theme(legend.position="none") +
  scale_x_discrete(limits=ID) + 
  coord_cartesian(ylim=c(50 ,100)) +
  theme(strip.text.x = element_text(size = 17))
RNA_p


#########
## CPU ##
#########
#raw data
CPU <- tibble(
  raw_dp = c(20.0771809, 29.99473011, 40.04825189, 50.1167232, 60.19903304, 75.28813814),
  raw_cpu = c(1563, 4162, 6363, 10693, 12386, 32950),
  correct_dp = c(11.0514682, 21.1812726, 29.0648353, 36.5075149, 39.1364557, 44.4026851),
  correct_cpu = c(1860, 4036, 5959, 7914, 8849, 11520)
)
str(CPU)

#curve fitting using https://mycurvefit.com/
raw <- data.frame(dp = c(20:76))
raw$cpu <- 20603100000 + (3136.685 - 20603100000)/(1 + (raw$dp/1932.377)^4.148144)
correct <- data.frame(dp = c(10:45))
correct$cpu <- 6438752000 + (1284.689 - 6438752000)/(1 + (correct$dp/56334.74)^1.872455)

cols <- c("Falcon correction"="#649EFC", "Canu assembly"="#F57670", "Curve fitting"="grey")
CPU_p = ggplot() + 
  geom_line(raw, mapping=aes(x=dp, y=cpu, colour="Curve fitting"), linetype = 2, size=1) +
  geom_point(CPU, mapping=aes(x=raw_dp, y=raw_cpu, colour="Falcon correction"),size=5) +
  geom_line(correct, mapping=aes(x=dp, y=cpu, colour="Curve fitting"), linetype = 2, size=1) +
  geom_point(CPU, mapping=aes(x=correct_dp, y=correct_cpu, colour="Canu assembly"),size=5) +
  coord_cartesian(xlim=c(0,80), ylim=c(0,40000)) +
  scale_y_continuous(breaks=c(0, 20000, 40000)) +
  scale_colour_manual(values=cols, breaks=names(cols)) +
  guides(fill = guide_legend(reverse=TRUE)) +
  labs(x ="Read depth (x)\n\n", y = "CPU core hour")+
  theme_classic() +
  theme(axis.title.x = element_text(size=20, face="bold",vjust=-1),
        axis.title.y = element_text(size=20, face="bold"),
        axis.text.y = element_text(size=18, angle=90, hjust=0.5), 
        axis.text.x = element_text(size=18, vjust=1,hjust=0.5),
    #    axis.text.x = element_text(face="bold", size=22, vjust=1,hjust=0.5),
        legend.text = element_text(size=16),legend.title=element_blank()) +
  theme(strip.text.x = element_text(size = 17)) +
  theme(strip.text.y = element_text(size = 17)) +
  theme(legend.position = c(0.05, 1),
        legend.background = element_rect(fill=NA),
        legend.justification = c("left", "top"))


CPU_p


################
## Bionano SV ##
################
datadir="/Users/oushujun/Google Drive/study/Maize research/NC358/optical map/"
setwd(datadir)
sv_files <- list.files(path = "Pseudomolecules_SV_3.4",pattern="*_3.4.sv")
total<-data.frame(matrix(ncol=0,nrow=0))  
for (files in sv_files)
{
  svs <- read.table(paste0("./Pseudomolecules_SV_3.4/", files),stringsAsFactors=FALSE)
  svs_id <- gsub(pattern = "_", replacement = "_",substr(files, 7, 13))
  #print(svs_id)
  svs_id_df <- cbind(svs, rep(c(svs_id), times = dim(svs)[1])) 
  #print(dim(svs_id_df))
  total <- rbind(total,svs_id_df)
}

colnames(total) <- c("Chromosome", "Start","End","SV_Type","Size","Line")
total$SV_Type<- replace(total$SV_Type, total$SV_Type == "deletion", "test")
total$SV_Type<- replace(total$SV_Type, total$SV_Type == "insertion", "deletion")
total$SV_Type<- replace(total$SV_Type, total$SV_Type == "test", "insertion")
total<-total[- grep("inv", total$SV_Type),]
str(total)

SV_p = ggplot(total, aes(x=Line, y=Size/1000,color=SV_Type,fill=SV_Type)) +
  geom_jitter(cex=1.3, position =position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) + 
  scale_y_continuous(trans='log10') +
  scale_x_discrete(limits=ID) + 
  theme_classic() +  
  theme(axis.title.y = element_text(size=20, face="bold"),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(face="bold", size=20, angle=35,vjust=1,hjust=0.95),
        legend.text = element_text(size=18),
        legend.title=element_blank(),
        legend.spacing.x = unit(0.3, 'cm')) +
  theme(legend.position="top") +
  labs(x =" ", y = "Bionano SV Size (kb)",size=14, face="bold") +
  geom_boxplot(color="black",lwd=0.3, outlier.colour = NA,alpha = 0.0,position = position_dodge(0.8)) +scale_color_manual(values=c("#5F9EA0","#CD853F","#003366")) +
  theme(legend.position = c(0.25, 1),
        legend.direction = "horizontal",
        legend.justification = c("left", "top"))
SV_p



######################
## make final plots ##
######################

pdf("Figure1_2.pdf", width=15,height=16,pointsize=12, paper='special')
lay = rbind(c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
             c(NA,NA,NA,NA,NA,1,1,1,2,2,2,2),
             c(NA,NA,NA,NA,NA,3,3,3,3,3,3,3))
lay2 = rbind(c(1,1,1,2,2,3,3,3,3))
lay3 = rbind(c(1),
             c(1),
             c(1),
             c(2))
g1 = arrangeGrob(grobs = list(BUSCO_p,contig_p,LAI_p), layout_matrix = lay)
g2 = arrangeGrob(grobs = list(RNA_p,CPU_p,SV_p), layout_matrix = lay2)
grid.arrange(g1, g2, layout_matrix=lay3)
dev.off()

