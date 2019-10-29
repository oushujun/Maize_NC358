library(tidyr)
library(ggplot2)
library(wesanderson)
library(scales)
library(reshape2)
library(plyr) 
library(grid)
library(gridExtra)
library(lattice)
#install.packages("devtools")
#devtools::install_github("thomasp85/patchwork")
library(patchwork)

#plot order
ID = c("21k_20x" ,"21k_30x" ,"21k_40x" ,"21k_50x" ,"21k_60x" ,"21k_75x", "11k_50x" ,"16k_50x")

##############
## LTR>26kb ##
##############
datadir="/Users/oushujun/Google Drive/study/Maize research/NC358/TEs/Len stat/"
setwd(datadir)

Fra_all = read.table("NC358_all_pseudomolecules_v1.fasta.mod.out.len",header=T)

#plot with jitter boxplot for 50kb > len > 20kb using ggplot
Fra_all_filter = subset(Fra_all, Len >= 26000)
gg_box_jitter = ggplot(data=Fra_all_filter, aes(x=Genome, y=Len)) + 
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2))
gg_box_jitter

head(Fra_all_filter)

Fra_all_p = Fra_all_filter %>% 
  gather(variable, value, Len) %>% 
  ggplot(aes(x = Genome, y = value)) + 
  geom_boxplot(colour="black") +
  geom_jitter(position=position_jitter(0.2)) + theme_classic() +
  labs(x =" ", y = "LTR sequence length",size=14, face="bold")+
  theme(axis.title.y = element_text(size=22, face="bold"),axis.text.y = element_text(size=18), axis.text.x = element_text(face="bold", size=22, angle=35,vjust=1,hjust=0.95),legend.text = element_text(size=18),legend.title=element_blank()) +
  theme(strip.text.y = element_text(size = 17)) +
  theme(legend.position="top") +
  scale_x_discrete(limits=ID) + 
  theme(strip.text.x = element_text(size = 17))
Fra_all_p


datadir="/Users/oushujun/Google Drive/study/Maize research/NC358/optical map/"
setwd(datadir)


############
## All TE ##
############
#repeats = read.csv('NC358_sum.csv')
repeats = read.table('NC358_repeat_sum.txt', header=T)
repeats$knob180 <- as.numeric(repeats$knob180) / 1000000
repeats$NOR <- as.numeric(repeats$NOR) / 1000000
repeats$CentC <- as.numeric(repeats$CentC) / 1000000
repeats$TR.1 <- as.numeric(repeats$TR.1) / 1000000
repeats$subtelomere <- as.numeric(repeats$subtelomere) / 1000000

colnames(repeats)[colnames(repeats)=="knob180"] <- "eKnob180"
colnames(repeats)[colnames(repeats)=="TR.1"] <- "dTR.1"
colnames(repeats)[colnames(repeats)=="NOR"] <- "cNOR"
colnames(repeats)[colnames(repeats)=="CentC"] <- "bCentC"
colnames(repeats)[colnames(repeats)=="subtelomere"] <- "asubtelomere"

head(repeats)
repeats_p = repeats %>% 
  gather(variable, value, eKnob180,dTR.1,cNOR,bCentC,asubtelomere) %>% 
  ggplot(aes(x = Assembly, y = value, fill = variable)) + 
  geom_bar(stat = "identity",colour="black") +
  scale_fill_manual(values = c("white","grey","#d8b365","#5F9EA0","#697184"),labels = c("Subtelomere","CentC", "NOR", "TR-1", "Knob180"),guide = guide_legend(reverse=TRUE)) + theme_classic() +
  labs(x =" ", y = "Repeat contents (Mb)",size=14, face="bold")+
  theme(axis.title.y = element_text(size=22, face="bold"),axis.text.y = element_text(size=18), axis.text.x = element_text(face="bold", size=22, angle=35,vjust=1,hjust=0.95),legend.text = element_text(size=18),legend.title=element_blank()) +
  scale_y_continuous(labels=comma) +
  theme(strip.text.y = element_text(size = 17)) +
  theme(legend.position="top") +
  scale_x_discrete(limits=ID) + 
  theme(strip.text.x = element_text(size = 17))
repeats_p


###########
### LTR ###
###########
datadir="/Users/oushujun/Google Drive/study/Maize research/NC358/TEs"
setwd(datadir)
LTR = read.table('NC358.LTR.txt', header=T)

#combine chrs to genome
data = data.frame(matrix(ncol = 4))
names(data) = c("Assembly", "Assembly.size", "Ngap", "Bionano.size")
for(i in ID){
  data = na.omit(rbind(data, c(paste(i), colSums(LTR[LTR$Assembly==i, 5:7]))))
}
data

#redefine data type
LTR = data
LTR$Assembly.size = as.numeric(LTR$Assembly.size)
LTR$LTR = as.numeric(LTR$LTR)
LTR$Ngap = as.numeric(LTR$Ngap)
LTR$Bionano.size = as.numeric(LTR$Bionano.size)

#convert to percentage
LTR$LTR <- 100*as.numeric(LTR$Assembly.size) / LTR$Bionano.size
LTR$Ngap <- 100*as.numeric(LTR$Ngap) / LTR$Bionano.size
head(LTR)

#make plot
LTR_p = LTR %>% 
  gather(variable, value, Ngap, LTR) %>% 
  ggplot(aes(x = Assembly, y = value, fill = variable, group=c(LTR$Ngap, LTR$LTR))) + 
  geom_bar(stat = "identity",colour="black") +
  scale_fill_manual(values = c("#5F9EA0","white"),labels = c("LTR","Ngap"),guide = guide_legend(reverse=F)) + theme_classic() +
  labs(x =" ", y = "Bionano Percentage (%)",size=14, face="bold")+
  theme(axis.title.y = element_text(size=22, face="bold"),axis.text.y = element_text(size=18), axis.text.x = element_text(face="bold", size=22, angle=35,vjust=1,hjust=0.95),legend.text = element_text(size=18),legend.title=element_blank()) +
  scale_y_continuous(labels=comma, breaks=c(0,25,50,75,100)) + 
  coord_cartesian(ylim=c(0 ,110)) +
  theme(strip.text.y = element_text(size = 17)) +
  theme(legend.position="top") +
  scale_x_discrete(limits=ID) + 
  theme(strip.text.x = element_text(size = 17))
LTR_p


###########
## knobs ##
###########
datadir="/Users/oushujun/Google Drive/study/Maize research/NC358/optical map/"
setwd(datadir)
knobs = read.table('NC358_knobs.txt', header=T)

knobs$Other <- 100*(knobs$Size - knobs$Knob180 - knobs$TR.1 - knobs$NGap) / knobs$Bionano.size
knobs$Knob180 <- 100*as.numeric(knobs$Knob180) / knobs$Bionano.size
knobs$TR.1 <- 100*as.numeric(knobs$TR.1) / knobs$Bionano.size
knobs$NGap <- 100*as.numeric(knobs$NGap) / knobs$Bionano.size

knobs$Chromosome<- replace(knobs$Chromosome, knobs$Chromosome == "4", "Chr4 Knob: 2.72Mb")
knobs$Chromosome<- replace(knobs$Chromosome, knobs$Chromosome == "7", "Chr7 Knob: 2.38Mb")
colnames(knobs)[colnames(knobs)=="Knob180"] <- "eKnob180"
colnames(knobs)[colnames(knobs)=="TR.1"] <- "dTR.1"
colnames(knobs)[colnames(knobs)=="Other"] <- "cOther"
colnames(knobs)[colnames(knobs)=="NGap"] <- "bNGap"

head(knobs)
knobs_p = knobs %>% 
  gather(variable, value, eKnob180, dTR.1, cOther, bNGap) %>% 
  ggplot(aes(x = Assembly, y = value, fill = variable)) + 
  geom_bar(stat = "identity",colour="black") +
  scale_fill_manual(values = c("white","grey","#d8b365","#5F9EA0"),labels = c("Ngap", "Other", "TR-1", "Knob180"),guide = guide_legend(reverse=TRUE)) + theme_classic() +
  labs(x =" ", y = "Bionano Percentage (%)",size=14, face="bold")+
  theme(axis.title.y = element_text(size=22, face="bold"),axis.text.y = element_text(size=18), axis.text.x = element_text(face="bold", size=22, angle=35,vjust=1,hjust=0.95),legend.text = element_text(size=16),legend.title=element_blank()) +
  scale_y_continuous(labels=comma) + facet_grid(Chromosome ~ .) +
  theme(strip.text.y = element_text(size = 17)) +
  theme(legend.position="top") +
  scale_x_discrete(limits=ID) + 
  facet_wrap(~Chromosome,ncol=1) +
  theme(strip.text.x = element_text(size = 17))
knobs_p


###########
## gene ##
###########
datadir="/Users/oushujun/Google Drive/study/Maize research/NC358/optical map/"
setwd(datadir)
gene = read.table('NC358_tandemGene.txt', header=T)
str(gene)

gene$Bionano.size = as.numeric(gene$Bionano.size)
gene$Assembled <- 100*as.numeric(gene$Assembled - gene$Gap) / gene$Bionano.size
gene$NGap <- 100*as.numeric(gene$Gap) / gene$Bionano.size

colnames(gene)[colnames(gene)=="Assembled"] <- "Assembled"
colnames(gene)[colnames(gene)=="NGap"] <- "NGap"

gene$Locus = as.character(gene$Locus)
gene$Locus = replace(gene$Locus, gene$Locus == "zein", "zein: 62kb")
gene$Locus = replace(gene$Locus, gene$Locus == "Rp1-D", "Rp1-D: 536kb")

head(gene)

gene_p = gene %>% 
  gather(variable, value, Assembled, NGap) %>% 
  ggplot(aes(x = Assembly, y = value, fill = variable, group=c(gene$Gap, gene$Assembled))) + 
  geom_bar(stat = "identity",colour="black") +
  scale_fill_manual(values = c("#5F9EA0","white"),labels = c("Tandem gene array","Ngap"), guide = guide_legend(reverse=F)) + theme_classic() +
  labs(x =" ", y = "Bionano Percentage (%)",size=14, face="bold")+
  theme(axis.title.y = element_text(size=22, face="bold"), axis.text.y = element_text(size=18), axis.text.x = element_text(face="bold", size=22, angle=35,vjust=1,hjust=0.95),legend.text = element_text(size=18),legend.title=element_blank()) +
  scale_y_continuous(labels=comma) + facet_grid(Locus ~ .) + 
  coord_cartesian(ylim=c(0 ,100)) +
  theme(strip.text.y = element_text(size = 17)) +
  theme(legend.position="top") +
  scale_x_discrete(limits=ID) + 
  facet_wrap(~Locus,ncol=1) +
  theme(strip.text.x = element_text(size = 17))
gene_p


##############
## telomere ##
##############
telomere = read.table('NC358_telomere7mer.txt', header=T)
str(telomere)

telomere$All = replace(telomere$All, telomere$All == "All", "Telomere 7-mer")
telomere$All = as.numeric(telomere$All)

head(telomere)
str(telomere)

#make plot
telomere_p = telomere %>% 
  gather(variable, value, All) %>% 
  ggplot(aes(x = Assembly, y = value, fill = "#5F9EA0")) + 
  geom_bar(stat = "identity",colour="black") +
  scale_fill_manual(values = c("#5F9EA0"),labels = c("Telomere 7-mer"),guide = guide_legend(reverse=TRUE)) + theme_classic() +
  labs(x =" ", y = "Telomere 7-mer count",size=14, face="bold")+
  theme(axis.title.y = element_text(size=22, face="bold"),axis.text.y = element_text(size=18), axis.text.x = element_text(face="bold", size=22, angle=35,vjust=1,hjust=0.95),legend.text = element_text(size=18),legend.title=element_blank()) +
  theme(strip.text.y = element_text(size = 17)) +
  theme(legend.position="top") +
  scale_x_discrete(limits=ID) + 
  theme(strip.text.x = element_text(size = 17))
telomere_p


###########
## centC ##
###########
centC = read.table('NC358_centC.txt', header=T)

#combine chrs to genome
data = data.frame(matrix(ncol = 5))
names(data) = c("Assembly", "CentC", "Ngap", "Assembly.size", "Bionano.size")
for(i in ID){
  data = na.omit(rbind(data, c(paste(i), colSums(centC[centC$Assembly==i, 5:8]))))
}
data

#redefine data type
centC = data
centC$Assembly.size = as.numeric(centC$Assembly.size)
centC$CentC = as.numeric(centC$CentC)
centC$Ngap = as.numeric(centC$Ngap)
centC$Bionano.size = as.numeric(centC$Bionano.size)

#convert to percentage
centC$Other <- 100*(as.numeric(centC$Assembly.size) - as.numeric(centC$CentC) - as.numeric(centC$Ngap)) / centC$Bionano.size
centC$CentC <- 100*as.numeric(centC$CentC) / centC$Bionano.size
centC$Ngap <- 100*as.numeric(centC$Ngap) / centC$Bionano.size
head(centC)

#make plot
centC_p = centC %>% 
  gather(variable, value, CentC, Other, Ngap) %>% 
  ggplot(aes(x = Assembly, y = value, fill = variable)) + 
  geom_bar(stat = "identity",colour="black") +
  scale_fill_manual(values = c("white","grey","#5F9EA0"),labels = c("Ngap","Other","CentC"),guide = guide_legend(reverse=TRUE)) + theme_classic() +
  labs(x =" ", y = "Bionano Percentage (%)",size=14, face="bold")+
  theme(axis.title.y = element_text(size=22, face="bold"),axis.text.y = element_text(size=18), axis.text.x = element_text(face="bold", size=22, angle=35,vjust=1,hjust=0.95),legend.text = element_text(size=18),legend.title=element_blank()) +
  scale_y_continuous(labels=comma) + 
  coord_cartesian(ylim=c(0 ,100)) +
  theme(strip.text.y = element_text(size = 17)) +
  theme(legend.position="top") +
  scale_x_discrete(limits=ID) + 
  theme(strip.text.x = element_text(size = 17))
centC_p


#####################
## ChIP centromere ##
#####################
cent_chip = na.omit(read.table('NC358_ChIP_centromere.txt', header=T))
head(cent_chip)

#combine chrs to genome
data = data.frame(matrix(ncol = 4))
names(data) = c("Assembly", "Assembly.size", "CentC", "Ngap")
for(i in ID){
  data = na.omit(rbind(data, c(paste(i), colSums(cent_chip[cent_chip$Assembly==i, 6:8]))))
}
data

#redefine data type
cent_chip = data
cent_chip$Assembly.size = as.numeric(cent_chip$Assembly.size)/1000000
cent_chip$CentC = as.numeric(cent_chip$CentC)/1000000
cent_chip$Ngap = as.numeric(cent_chip$Ngap)/1000000
cent_chip$Other <- cent_chip$Assembly.size - cent_chip$CentC - cent_chip$Ngap

head(cent_chip)

#make plot
cent_chip_p = cent_chip %>% 
  gather(variable, value, CentC, Other, Ngap) %>% 
  ggplot(aes(x = Assembly, y = value, fill = variable)) + 
  geom_bar(stat = "identity",colour="black") +
  scale_fill_manual(values = c("white","grey","#5F9EA0"),labels = c("Ngap","Other","CentC"),guide = guide_legend(reverse=TRUE)) + theme_classic() +
  labs(x =" ", y = "Centromere size (Mb)",size=14, face="bold")+
  theme(axis.title.y = element_text(size=22, face="bold"),axis.text.y = element_text(size=18), axis.text.x = element_text(face="bold", size=22, angle=35,vjust=1,hjust=0.95),legend.text = element_text(size=18),legend.title=element_blank()) +
  scale_y_continuous(labels=comma) + 
  theme(strip.text.y = element_text(size = 17)) +
  theme(legend.position="top") +
  scale_x_discrete(limits=ID) + 
  theme(strip.text.x = element_text(size = 17))
cent_chip_p


#########
## NOR ##
#########
NOR = read.table('NC358_NOR.txt', header=T)

#convert to percentage
NOR$Bionano.size = as.numeric(NOR$Bionano.size)
NOR$NOR <- 100*as.numeric(NOR$Assembled) / NOR$Bionano.size
NOR$Ngap <- 100*as.numeric(NOR$Ngap) / NOR$Bionano.size
head(NOR)

#make plot
NOR_p = NOR %>% 
  gather(variable, value, NOR, Ngap) %>% 
  ggplot(aes(x = Assembly, y = value, fill = variable)) + 
  geom_bar(stat = "identity",colour="black") +
  scale_fill_manual(values = c("white","#5F9EA0"),labels = c("Ngap","NOR"),guide = guide_legend(reverse=TRUE)) + theme_classic() +
  labs(x =" ", y = "Bionano Percentage (%)",size=14, face="bold")+
  theme(axis.title.y = element_text(size=22, face="bold"),axis.text.y = element_text(size=18), axis.text.x = element_text(face="bold", size=22, angle=35,vjust=1,hjust=0.95),legend.text = element_text(size=18),legend.title=element_blank()) +
  scale_y_continuous(labels=comma) + 
  coord_cartesian(ylim=c(0 ,100)) +
  theme(strip.text.y = element_text(size = 17)) +
  theme(legend.position="top") +
  scale_x_discrete(limits=ID) + 
  theme(strip.text.x = element_text(size = 17))
NOR_p


########################
### make final plots ###
########################

pdf("Figure2_1.pdf", width=15,height=5,pointsize=12, paper='special')
grid.arrange(repeats_p,
             Fra_all_p, ncol=2)
dev.off()

pdf("Figure2_2.pdf", width=15,height=11,pointsize=12, paper='special')
grid.arrange(telomere_p,
             cent_chip_p,
             centC_p,
             NOR_p,
             gene_p,
             knobs_p, ncol=3)
dev.off()

pdf("Figure2.pdf", width=15,height=15,pointsize=12, paper='special')
lay <- rbind(c(1,1,1,2,2,2),
             c(3,3,4,4,5,5),
             c(6,6,7,7,8,8))
grid.arrange(grobs = list(repeats_p, Fra_all_p, telomere_p, cent_chip_p, centC_p, NOR_p, gene_p, knobs_p),
             layout_matrix = lay)

plot = (repeats_p | Fra_all_p) / (telomere_p | cent_chip_p | centC_p) / (NOR_p | gene_p | knobs_p) +
  plot_annotation(tag_levels = "A", theme(element_text(size=22, face="bold")))
plot

dev.off()

pdf("Figure2_3.pdf", width=15,height=15,pointsize=12, paper='special')
lay <- rbind(c(1,1,1,2,2,2),
             c(3,3,4,4,5,5),
             c(6,6,7,7,8,8))
grid.arrange(grobs = list(repeats_p, Fra_all_p, LTR_p, telomere_p, centC_p, NOR_p, gene_p, knobs_p),
             layout_matrix = lay)
dev.off()



