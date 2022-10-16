#Load the packages
library (ggplot2)
library (ggpubr)
library (dplyr)
library (tidyr)
library(stringr)
library("vegan")
library(readr)
library(reshape2)
library(RColorBrewer)
library("tibble")
library(lmerTest)
library(lme4)

#Basic Process
ARGenes <- read.delim("ARGs_table.txt", header=TRUE) %>% 
separate (., Resistance.Type.Gene, into=c( "Resistance_Type", "Gene"), sep = ";")%>%
  select( -DSgene) %>% 
  unite (., "united",  Resistance_Type:Gene:Sample, sep = ";") %>%  
aggregate(., Abundance~united, sum) %>% 
separate (., "united", into = c ("Resistance_Type", "Gene", "Sample"), sep = ";") 
ARGTotal_counts = ARGenes%>%select( -Gene, -Resistance_Type) %>% 
 aggregate(., Abundance~Sample, sum)
colnames(ARGTotal_counts) <- c("Sample", "Total_ARGs")
#Include metadata from the different studies
metadata <- read_delim("metadata.csv", ",", escape_double = FALSE, trim_ws = TRUE)
Taxonomy =  read.delim("global.taxonomy.table.txt", header=TRUE, sep = ",")
Bacterial_Taxonomy =  filter(Taxonomy, str_detect(Taxonomy, "Bacteria"))
Bacterial_Taxonomy =Bacterial_Taxonomy %>% gather(key ="Sample", value = "Abundance", -Taxonomy)
#Bacterial Counts
Bcounts = Bacterial_Taxonomy %>% select(Sample, Abundance)  %>% 
  aggregate(., Abundance~Sample, sum)
colnames(Bcounts) <- c("Sample", "Bacterial_Counts")
Fungal_Taxonomy =  filter(Taxonomy, str_detect(Taxonomy, "Fungi")) %>% 
 gather(., key ="Sample", value = "Abundance", -Taxonomy) %>% 
separate (., "Taxonomy", into = c ("Domain", "Kingdom", "Phylum", 
                                   "Class", "Order", "Family", "Genus", "Species"), sep = ";")
Fcounts = Fungal_Taxonomy %>% select(Sample, Abundance) %>% 
 aggregate( ., Abundance~Sample, sum)
colnames(Fcounts) <- c("Sample", "Fungal_Counts")
total_counts = full_join(Fcounts, Bcounts, by="Sample")
total_table = full_join(total_counts, ARGenes, by="Sample") %>% 
 full_join(., metadata, by="Sample") %>% 
filter (., Sequencing %in% c ("Metagenome")) %>% 
 mutate(., RA = (Abundance/Bacterial_Counts)) %>% 
filter (., Abundance > 1)
total_table = unite(total_table, "united", Gene, Location, sep = ",")
total_table = total_table[unsplit(table(total_table$united), total_table$united) > 3, ]
total_table = separate(total_table, "united", into = c("Gene", "Location"), sep = ",")
NMDS_Location = select(total_table, Location, Location_Sample, Gene, RA)
NMDS_Location2 = pivot_wider(NMDS_Location, names_from = Gene, values_from = RA, values_fill = 0.00000001, values_fn = sum)
CategoricalVar. = select(NMDS_Location2, Location, Location_Sample)
gathered = gather(NMDS_Location2, -Location, -Location_Sample, key = Gene, value=RA)
#Figure 1: ARG abundance  of selected ARGs
gathered2 = filter (gathered, Gene %in% c("blaTEM", "aph(3')", "aph(3'')", "sul1", "sul2", "blaOXA"))
ARGboxplot=ggplot(data = gathered2, aes (y=log10(RA), x=Gene, fill=Location ))  + geom_boxplot() + #geom_errorbar(stat = "summary", position = "dodge", width = 1) +
  theme_pubr(legend = c("top")) + scale_fill_brewer(palette = "Set1") +
  theme( axis.text.y = element_text( size=40, face = "bold", colour="black")) +
  theme(axis.text.x = element_text(size = 30, face = "bold", colour="black", angle = 0, hjust = 0)) +
  theme(legend.text = element_text( size=40, face="bold", colour = "black")) + 
  theme( axis.title.y = element_text( size=40, face = "bold", colour="black")) + 
  theme(axis.title.x = element_text(size = 40, face = "bold", colour="black")) +
  theme(legend.title =  element_text( size=40, face="bold", colour = "black"))  + ylab(" log10 ARG hits/Bacterial SSU") + xlab("Gene") + stat_compare_means (method="kruskal", label= "p.signif", size=10)
ARGboxplot

# NMDS ARG profile, Figure 2A
NMDS_Location2 = select(NMDS_Location2, -Location, -Location_Sample)
NMDS_Location2 = log10(NMDS_Location2)
mMDS<-metaMDS(NMDS_Location2, distance="euclidean", k=2, trymax=90)
NMDS_1 <- mMDS$points[,1] 
NMDS_2 <- mMDS$points[,2]
NMDS.plot<-cbind(NMDS_Location2, NMDS_1, NMDS_2, CategoricalVar.)
p<-ggplot(NMDS.plot, aes(NMDS_1, NMDS_2, color=Location))+
  geom_point( size=10) +   theme_pubr() + theme( axis.text.x = element_text( size=40, face = "bold", colour="black"))   +
  theme( axis.text.y = element_text( size=40, face = "bold", colour="black")) + 
  theme(axis.text.x = element_text(size = 40, face = "bold", colour="black")) +
  theme(legend.text = element_text( size=40, face="bold", colour = "black")) + 
  theme(legend.title = element_text(size = 40, face="bold", colour = "black")) +
  theme(axis.title = element_text(size =40, face="bold")) + 
  theme(plot.title = element_text(size=10, face="bold", hjust = 0.5)) +
  theme(axis.title.x = element_text(size=10, face="bold")) +
  theme(axis.title.x =  element_text(size = 40, face = "bold", colour="black"))  + scale_colour_brewer (palette = "Set1")
p

#PERMANOVA tests for FIGURE 2A
dist <- vegdist(NMDS_Location2, method='euclidean')
Comparison.div<-adonis2(dist~CategoricalVar.$Location, data=NMDS_Location2, permutations = 999999, method="euclidean")
Comparison.div

gathered2 = gathered
gathered2$RA = log10(gathered2$RA)
#Mean values
mean_location = dcast(gathered2, Gene~Location, value.var = "RA", mean)
#SD values
sd_location = dcast(gathered2, Gene~Location, value.var = "RA", sd)

#Select grouping based on Bacterial Orders
Taxonomy =  read.delim("global.taxonomy.table.txt", header=TRUE, sep = ",")
Bacterial_Taxonomy =  filter(Taxonomy, str_detect(Taxonomy, "Bacteria"))
Bacterial_Taxonomy =Bacterial_Taxonomy %>% gather(key ="Sample", value = "Abundance", -Taxonomy)
#Bacterial Counts
Bacterial_Orders =   Bacterial_Taxonomy%>% separate(., Taxonomy, into=c("Kingdom", 
                                                                        "Phylum",
                                                                        "Class", 
                                                                        "Order", 
                                                                        "Family",
                                                                        "Genus", 
                                                                        "Species"), sep=";")%>% 
unite ( "Order", Kingdom:Order, sep = ";") %>%
dcast (., Order~Sample, value.var = "Abundance", sum )  %>%
 gather (., -Order, key=Sample, value=Abundance_Bacteria)

total_counts = full_join(Bacterial_Orders, Bcounts, by="Sample")
total_table = full_join(total_counts, metadata, by="Sample") %>%
 filter (., Sequencing %in% c ("Metagenome")) %>%
mutate(., RA = ((Abundance_Bacteria)/Bacterial_Counts )*100) %>%
 unite(., "united", Order, Location, sep =",") %>%
 filter (., Abundance_Bacteria > 1)
total_table = total_table[unsplit(table(total_table$united), total_table$united) > 3, ]
total_table = separate(total_table, "united", into = c("Order", "Location"), sep = ",")

#Figure 2B. NMDS with Bacterial Abundance over different countries
NMDS_Location = select(total_table, Location, Location_Sample, Order, RA)
NMDS_Location2 = pivot_wider(NMDS_Location, names_from = Order, values_from = RA, values_fill = 0, values_fn = sum)
CategoricalVar. = select(NMDS_Location2, Location, Location_Sample)
gathered = gather(NMDS_Location2, -Location, -Location_Sample, key = Order, value=RA)
NMDS_Location2 = select(NMDS_Location2, -Location, -Location_Sample)
NMDS_Location2 = log10(NMDS_Location2+10^-8)+8
NMDS_Location2 = na.omit(NMDS_Location2)
mMDS3<-metaMDS(NMDS_Location2, distance="bray", k=2, trymax=35)
NMDS_1 <- mMDS3$points[,1] 
NMDS_2 <- mMDS3$points[,2]
NMDS.plot<-cbind(NMDS_Location2, NMDS_1, NMDS_2, CategoricalVar.)
NMDS.plot$Location = factor(NMDS.plot$Location, c("China", "Japan", "Saudi Arabia", "USA", "Germany"))
p<-ggplot(NMDS.plot, aes(NMDS_1, NMDS_2, color=Location))+
  geom_point( size=10) +   theme_pubr() + theme( axis.text.x = element_text( size=40, face = "bold", colour="black"))   +
  theme( axis.text.y = element_text( size=40, face = "bold", colour="black")) + 
  theme(axis.text.x = element_text(size = 40, face = "bold", colour="black")) +
  theme(legend.text = element_text( size=40, face="bold", colour = "black")) + 
  theme(legend.title = element_text(size = 40, face="bold", colour = "black")) +
  theme(axis.title = element_text(size =40, face="bold")) + 
  theme(plot.title = element_text(size=10, face="bold", hjust = 0.5)) +
  theme(axis.title.x = element_text(size=10, face="bold")) +
  theme(axis.title.x =  element_text(size = 40, face = "bold", colour="black"))  + scale_colour_brewer (palette = "Set1")
p

dist <- vegdist(NMDS_Location2, method='bray')
#PERMANOVA Test
Comparison.div<-adonis2(dist~CategoricalVar.$Location, data=NMDS_Location2, permutations = 99999, method="bray")
Comparison.div

Bacterial_Taxonomy =  filter(Taxonomy, str_detect(Taxonomy, "Bacteria"))
Bacterial_Taxonomy =Bacterial_Taxonomy%>% gather(key ="Sample", value = "Abundance", -Taxonomy)
Bcounts = Bacterial_Taxonomy %>% select(Sample, Abundance)
Bcounts =  aggregate( data=Bcounts, Abundance~Sample, sum)
colnames(Bcounts) <- c("Sample", "Bacterial_Counts")

# Mantel and Procrustes Tests/Figure 2C
ARGenes <- read.delim("ARGs_table.txt", header=TRUE)
ARGenes =   ARGenes%>% separate (Resistance.Type.Gene, into=c( "Resistance_Type", "Gene"), sep = ";")
ARGenes = ARGenes%>%select( -DSgene)
ARGenes =   ARGenes%>% unite ( "united",  Resistance_Type:Gene:Sample, sep = ";") 
ARGenes =  aggregate( data=ARGenes, Abundance~united, sum)
ARGenes =   ARGenes%>% separate ( "united", into = c ("Resistance_Type", "Gene", "Sample"), sep = ";") 
ARGTotal_counts = ARGenes%>%select( -Gene, -Resistance_Type)
ARGTotal_counts =  aggregate( data=ARGTotal_counts, Abundance~Sample, sum)
colnames(ARGTotal_counts) <- c("Sample", "Total_ARGs")
metadata <- read_delim("metadata.csv", ",", escape_double = FALSE, trim_ws = TRUE)
Taxonomy =  read.delim("global.taxonomy.table.txt", header=TRUE, sep = ",")
Bacterial_Taxonomy =  filter(Taxonomy, str_detect(Taxonomy, "Bacteria"))
Bacterial_Taxonomy =Bacterial_Taxonomy %>% gather(key ="Sample", value = "Abundance", -Taxonomy)
Bcounts = Bacterial_Taxonomy %>% select(Sample, Abundance)
Bcounts =  aggregate( data=Bcounts, Abundance~Sample, sum)
colnames(Bcounts) <- c("Sample", "Bacterial_Counts")
Fungal_Taxonomy =  filter(Taxonomy, str_detect(Taxonomy, "Fungi"))
Fungal_Taxonomy =Fungal_Taxonomy%>% gather(key ="Sample", value = "Abundance", -Taxonomy)
Fungal_Taxonomy =   Fungal_Taxonomy%>% separate ( "Taxonomy", into = c ("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") 
Fcounts = Fungal_Taxonomy %>% select(Sample, Abundance)
Fcounts =  aggregate( data=Fcounts, Abundance~Sample, sum)
metadata <- read_delim("metadata.csv", ",", escape_double = FALSE, trim_ws = TRUE)
colnames(Fcounts) <- c("Sample", "Fungal_Counts")
total_counts = full_join(Fcounts, Bcounts, by="Sample")
total_table = full_join(total_counts, ARGenes, by="Sample")
total_table = full_join(total_table, metadata, by="Sample")
total_table = filter (total_table, Sequencing %in% c ("Metagenome"))
total_table = na.omit(total_table)
total_table = mutate(total_table, RA = (Abundance/Bacterial_Counts))
total_table = filter (total_table, Abundance > 1)
total_table = unite(total_table, "united", Gene, Location, sep = ",")
total_table = total_table[unsplit(table(total_table$united), total_table$united) > 3, ]
total_table = separate(total_table, "united", into = c("Gene", "Location"), sep = ",")

NMDS_Location = select(total_table, Location, Location_Sample, Gene, Sample, RA)
NMDS_Location2 = pivot_wider(NMDS_Location, names_from = Gene, values_from = RA, values_fill = 0.00000001, values_fn = sum)
CategoricalVar. = select(NMDS_Location2, Location, Location_Sample )
gathered = gather(NMDS_Location2, -Location, -Location_Sample, key = Gene, value=RA)
NMDS_Location2ARG = select(NMDS_Location2, -Location, -Location_Sample)
#NMDS_Location2ARG = log10(NMDS_Location2)

Bacterial_Taxonomy =  filter(Taxonomy, str_detect(Taxonomy, "Bacteria"))
Bacterial_Taxonomy =Bacterial_Taxonomy%>% gather(key ="Sample", value = "Abundance", -Taxonomy)
Bacterial_Taxonomy =   Bacterial_Taxonomy%>% separate ( "Taxonomy", into = c ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species" ), sep = ";") 
Bcounts = Bacterial_Taxonomy %>% select(Sample, Abundance)
Bcounts =  aggregate( data=Bcounts, Abundance~Sample, sum)
colnames(Bcounts) <- c("Sample", "Bacterial_Counts")
Bacterial_Orders =   Bacterial_Taxonomy%>% unite ( "Order", Kingdom:Order, sep = ";") 

Bacterial_Orders =    dcast (Bacterial_Orders, Order~Sample, value.var = "Abundance", sum ) 
Bacterial_Orders = gather (Bacterial_Orders, -Order, key=Sample, value=Abundance_Bacteria)

metadata <- read_delim("metadata.csv", ",", escape_double = FALSE, trim_ws = TRUE)
total_counts = full_join(Bacterial_Orders, Bcounts, by="Sample")
total_table = full_join(total_counts, metadata, by="Sample")
total_table = filter (total_table, Sequencing %in% c ("Metagenome"))
total_table = mutate(total_table, RA = ((Abundance_Bacteria)/Bacterial_Counts )*100)
total_table = unite(total_table, "united", Order, Location, sep =",")
total_table = filter (total_table, Abundance_Bacteria > 1)
total_table = total_table[unsplit(table(total_table$united), total_table$united) > 3, ]
total_table = separate(total_table, "united", into = c("Order", "Location"), sep = ",")
NMDS_Location = select(total_table, Location, Sample, Location_Sample, Order, RA)
NMDS_Location2 = pivot_wider(NMDS_Location, names_from = Order, values_from = RA, values_fill = 0, values_fn = sum)
CategoricalVar. = select(NMDS_Location2, Location, Location_Sample)
gathered = gather(NMDS_Location2, -Location, -Location_Sample, key = Order, value=RA)
NMDS_Location2 = select(NMDS_Location2, -Location, -Location_Sample)
#NMDS_Location2 = log10(NMDS_Location2) + 8
NMDS_Location2Bact = na.omit(NMDS_Location2)
joined = full_join(NMDS_Location2ARG, NMDS_Location2Bact,  by="Sample")
joined = na.omit(joined
)
bact_abund = joined[, 26:ncol(joined)]
arg_abund =   joined[, 2:25]
arg_abund = log10(arg_abund + 8)
dist.abund = vegdist(bact_abund, method = "bray")
dist.abund2 = vegdist(arg_abund, method = "euclidean")
abund_geo  = mantel(dist.abund, dist.abund2, method = "spearman", permutations = 999, na.rm = TRUE)
abund_geo
abund_geo2  = protest(dist.abund, dist.abund2, method = "spearman", permutations = 999, na.rm = TRUE)
abund_geo2  
pro.test <- procrustes(dist.abund,dist.abund2)

#PROCRUSTES plot 
ctest <- data.frame(rda1=pro.test$Yrot[,1],
                    rda2=pro.test$Yrot[,2],xrda1=pro.test$X[,1],
                    xrda2=pro.test$X[,2])

ggplot(ctest) +
  geom_point(aes(x=xrda1, y=xrda2), size=8, colour="skyblue") +
  geom_point(aes(x=rda1, y=rda2), size=8, colour="darkred") +
  geom_segment(aes(x=xrda1,y=xrda2,xend=rda1,yend=rda2),arrow=arrow(length=unit(0.2,"cm")))


#ARG correlations with Fungal/Bacterial ratio; Fig. 1 
ARGenes <- read.delim("ARGs_table.txt", header=TRUE)
ARGenes =   ARGenes%>% separate (Resistance.Type.Gene, into=c( "Resistance_Type", "Gene"), sep = ";")
ARGenes = ARGenes%>%select( -DSgene)
ARGenes =   ARGenes%>% unite ( "united",  Resistance_Type:Gene:Sample, sep = ";") 
ARGenes =  aggregate( data=ARGenes, Abundance~united, sum)
ARGenes =   ARGenes%>% separate ( "united", into = c ("Resistance_Type", "Gene", "Sample"), sep = ";") 
Taxonomy =  read.delim("global.taxonomy.table.txt", header=TRUE, sep = ",")
Bacterial_Taxonomy =  filter(Taxonomy, str_detect(Taxonomy, "Bacteria"))
Bacterial_Taxonomy =Bacterial_Taxonomy%>% gather(key ="Sample", value = "Abundance", -Taxonomy)
Bacterial_Taxonomy =   Bacterial_Taxonomy%>% separate ( "Taxonomy", into = c ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") 
Bcounts = Bacterial_Taxonomy %>% select(Sample, Abundance)
Bcounts =  aggregate( data=Bcounts, Abundance~Sample, sum)
colnames(Bcounts) <- c("Sample", "Bacterial_Counts")
Fungal_Taxonomy =  filter(Taxonomy, str_detect(Taxonomy, c("Fungi")))
Fungal_Taxonomy =Fungal_Taxonomy%>% gather(key ="Sample", value = "Abundance", -Taxonomy)
Fungal_Taxonomy =   Fungal_Taxonomy%>% separate ( "Taxonomy", into = c ("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") 
Fcounts = Fungal_Taxonomy %>% select(Sample, Abundance)
Fcounts =  aggregate( data=Fcounts, Abundance~Sample, sum)
metadata <- read_delim("metadata.csv", ",", escape_double = FALSE, trim_ws = TRUE)
colnames(Fcounts) <- c("Sample", "Fungal_Counts")
total_counts = full_join(Fcounts, Bcounts, by="Sample")
total_table = full_join(total_counts, ARGenes, by="Sample")
total_table = full_join(total_table, metadata, by="Sample")
total_table = filter (total_table, Sequencing %in% c ("Metagenome"))
total_table = filter (total_table, Gene %in% c ("blaTEM", "aph(3')", "aph(3'')", "sul1", "sul2", "blaOXA"))
total_table = mutate(total_table, RA = ((Abundance)/Bacterial_Counts ))
total_table = mutate(total_table, RAFungi = ((Fungal_Counts)/Bacterial_Counts ))
total_table = filter (total_table, Fungal_Counts > 1)
total_table = filter (total_table, Abundance > 1)

total_table$RAFungi = log10(total_table$RAFungi )
total_table$RA = log10(total_table$RA )

#Figure 3 Fungal/bacterial SSUs against ARGs
plot_fu = ggscatter( data = total_table, y="RA", x="RAFungi", add = "reg.line", conf.int = FALSE, size = 4, color = "royalblue4", add.params = list (color = "royalblue4", fill = "#c8d4e3" ))  + 
  stat_cor(method = "spearman", size=10)  + 
  xlab("log10 Fungal SSU hits/Total Bacterial SSU hits") + ylab("log10 ARG hits/Total Bacterial SSU hits")  + theme( axis.text.y = element_text( size=25, face = "bold", colour="black")) + theme(axis.text.x = element_text(size = 25, face = "bold", colour="black")) +
  theme(legend.text = element_text( size=25, face="bold", colour = "black")) + theme( axis.title.y = element_text( size=25, face = "bold", colour="black")) + theme(axis.title.x = element_text(size = 25, face = "bold", colour="black")) +
  theme(legend.title =  element_text( size=25, face="bold", colour = "black")) +
  facet_wrap(~Gene, ncol = 3) + theme(strip.text = element_text(size=15,  face="bold", colour = "white")) + theme( strip.background = element_rect( fill="skyblue4"))
plot_fu
#Select original_study for performing a mixed effect model: (ARG_Relative_Abundance~Fungal/Bacterial_SSU + (1|Original_Study))

library(stringr)
library(lme4)
library(lmerTest)
library(lmtest)

total_table$Original_Study= str_sub(total_table$Sample, end=-3)

mixed_effect_blatem = lmer(data=total_table%>%filter(., Gene=="blaTEM"), RA~RAFungi+(1|Original_Study))
anova(mixed_effect_blatem)
qqnorm(resid(mixed_effect_blatem))
qqline(resid(mixed_effect_blatem))
shapiro.test(resid(mixed_effect_blatem))

mixed_effect_blaOXA = lmer(data=total_table%>%filter(., Gene=="blaOXA"), RA~RAFungi+(1|Original_Study))
anova(mixed_effect_blaOXA)
qqnorm(resid(mixed_effect_blaOXA))
qqline(resid(mixed_effect_blaOXA))
shapiro.test(resid(mixed_effect_blaOXA))


#Figure S3 Actinobacteria to Fungal/Bacterial SSU 
ARGenes <- read.delim("ARGs_table.txt", header=TRUE)
ARGenes =   ARGenes%>% separate (Resistance.Type.Gene, into=c( "Resistance_Type", "Gene"), sep = ";")
ARGenes = ARGenes%>%select( -DSgene)
ARGenes =   ARGenes%>% unite ( "united",  Resistance_Type:Gene:Sample, sep = ";") 
ARGenes =  aggregate( data=ARGenes, Abundance~united, sum)
ARGenes =   ARGenes%>% separate ( "united", into = c ("Resistance_Type", "Gene", "Sample"), sep = ";") 
Taxonomy =  read.delim("global.taxonomy.table.txt", header=TRUE, sep = ",")
Bacterial_Taxonomy =  filter(Taxonomy, str_detect(Taxonomy, "Bacteria"))
Bacterial_Taxonomy =Bacterial_Taxonomy%>% gather(key ="Sample", value = "Abundance", -Taxonomy)
Bacterial_Taxonomy =   Bacterial_Taxonomy%>% separate ( "Taxonomy", into = c ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") 
Bcounts = Bacterial_Taxonomy %>% select(Sample, Abundance)
Bcounts =  aggregate( data=Bcounts, Abundance~Sample, sum)
colnames(Bcounts) <- c("Sample", "Bacterial_Counts")
Act_Taxonomy =  filter(Taxonomy, str_detect(Taxonomy, c("Actinobacteria")))
Act_Taxonomy =Act_Taxonomy%>% gather(key ="Sample", value = "Abundance", -Taxonomy)
Act_Taxonomy =   Act_Taxonomy%>% separate ( "Taxonomy", into = c ("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") 
ActCounts = Act_Taxonomy %>% select(Sample, Abundance)
ActCounts =  aggregate( data=ActCounts, Abundance~Sample, sum)
metadata <- read_delim("metadata.csv", ",", escape_double = FALSE, trim_ws = TRUE)
colnames(ActCounts) <- c("Sample", "Act_Counts")
total_counts = full_join(ActCounts, Bcounts, by="Sample")
total_table = full_join(total_counts, ARGenes, by="Sample")
total_table = full_join(total_table, metadata, by="Sample")
total_table = filter (total_table, Sequencing %in% c ("Metagenome"))
total_table = filter (total_table, Gene %in% c ("blaTEM", "aph(3')", "aph(3'')", "sul1", "sul2", "blaOXA"))
total_table = mutate(total_table, RA = ((Abundance)/Bacterial_Counts ))
total_table = mutate(total_table, RAFungi = ((Act_Counts)/Bacterial_Counts ))
total_table = filter (total_table, Act_Counts > 1)
total_table = filter (total_table, Abundance > 1)
total_table$RAFungi = log10(total_table$RAFungi )
total_table$RA = log10(total_table$RA )

#Correlations of Actinobacteria with ARGs
plot_act = ggscatter( data = total_table, y="RA", x="RAFungi", add = "reg.line", conf.int = FALSE, size = 4, color = "royalblue4", add.params = list (color = "royalblue4", fill = "#c8d4e3" ))  + 
  stat_cor(method = "spearman", size=10)  + 
  xlab("log10 Actinobacteria SSU hits/Total Bacterial SSU rRNA hits") + ylab("log10 ARG hits/Total Bacterial SSU hits")  + theme( axis.text.y = element_text( size=25, face = "bold", colour="black")) + theme(axis.text.x = element_text(size = 25, face = "bold", colour="black")) +
  theme(legend.text = element_text( size=25, face="bold", colour = "black")) + theme( axis.title.y = element_text( size=25, face = "bold", colour="black")) + theme(axis.title.x = element_text(size = 25, face = "bold", colour="black")) +
  theme(legend.title =  element_text( size=25, face="bold", colour = "black")) +
  facet_wrap(~Gene, ncol = 3) + theme(strip.text = element_text(size=15,  face="bold", colour = "white")) + theme( strip.background = element_rect( fill="skyblue4"))
plot_act


#Co-occurence network Figure S4
ARGenes <- read.delim("ARGs_table.txt", header=TRUE)
ARGenes =   ARGenes%>% separate (Resistance.Type.Gene, into=c( "Resistance_Type", "Gene"), sep = ";")
ARGenes = ARGenes%>%select( -DSgene)
ARGenes =   ARGenes%>% unite ( "united",  Resistance_Type:Gene:Sample, sep = ";") 
ARGenes =  aggregate( data=ARGenes, Abundance~united, sum)
ARGenes =   ARGenes%>% separate ( "united", into = c ("Resistance_Type", "Gene", "Sample"), sep = ";") 
ARGTotal_counts = ARGenes%>%select( -Gene, -Resistance_Type)
ARGTotal_counts =  aggregate( data=ARGTotal_counts, Abundance~Sample, sum)
colnames(ARGTotal_counts) <- c("Sample", "Total_ARGs")
metadata <- read_delim("metadata.csv", ",", escape_double = FALSE, trim_ws = TRUE)

Taxonomy =  read.delim("global.taxonomy.table.txt", header=TRUE, sep = ",")
Bacterial_Taxonomy =  filter(Taxonomy, str_detect(Taxonomy, "Bacteria"))
Bacterial_Taxonomy =Bacterial_Taxonomy %>% gather(key ="Sample", value = "Abundance", -Taxonomy)
Bcounts = Bacterial_Taxonomy %>% select(Sample, Abundance)
Bcounts =  aggregate( data=Bcounts, Abundance~Sample, sum)
colnames(Bcounts) <- c("Sample", "Bacterial_Counts")
metadata <- read_delim("metadata.csv", ",", escape_double = FALSE, trim_ws = TRUE)
colnames(Fcounts) <- c("Sample", "Fungal_Counts")
total_table = full_join(Bcounts, ARGenes, by="Sample")
total_table = full_join(total_table, metadata, by="Sample")
total_table = filter (total_table, Sequencing %in% c ("Metagenome"))
total_table = mutate(total_table, RA = (Abundance/Bacterial_Counts))
total_table = filter (total_table, Abundance > 1)
total_table = unite(total_table, "united", Gene, Location, sep = ",")
total_table = separate(total_table, "united", into = c("Gene", "Location"), sep = ",")
NMDS_Location = select(total_table, Location, Location_Sample, Gene, Sample, RA)
NMDS_Location2 = pivot_wider(NMDS_Location, names_from = Gene, values_from = RA, values_fill = 0, values_fn = sum)
CategoricalVar. = select(NMDS_Location2, Location, Location_Sample )
gathered = gather(NMDS_Location2, -Location, -Location_Sample, key = Gene, value=RA)
NMDS_Location2ARG = select(NMDS_Location2, -Location, -Location_Sample)
Taxonomy =  read.delim("global.taxonomy.table.txt", header=TRUE, sep = ",")
Bacterial_Taxonomy =  filter(Taxonomy, str_detect(Taxonomy, "Bacteria"))
Bacterial_Taxonomy =Bacterial_Taxonomy%>% gather(key ="Sample", value = "Abundance", -Taxonomy)
Bacterial_Taxonomy =   Bacterial_Taxonomy%>% separate ( "Taxonomy", into = c ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species" ), sep = ";") 
Bcounts = Bacterial_Taxonomy %>% select(Sample, Abundance)
Bcounts =  aggregate( data=Bcounts, Abundance~Sample, sum)
colnames(Bcounts) <- c("Sample", "Bacterial_Counts")
Bacterial_Orders =    dcast (Bacterial_Taxonomy, Phylum+Class+Order~Sample, value.var = "Abundance", sum ) 
Bacterial_Orders = gather (Bacterial_Orders, -Phylum,-Class,-Order, key=Sample, value=Abundance_Bacteria)
metadata <- read_delim("metadata.csv", ",", escape_double = FALSE, trim_ws = TRUE)
total_counts = full_join(Bacterial_Orders, Bcounts, by="Sample")
total_table = full_join(total_counts, metadata, by="Sample")
total_table = filter (total_table, Sequencing %in% c ("Metagenome"))
total_table = mutate(total_table, RA = ((Abundance_Bacteria)/Bacterial_Counts )*100)
total_table = filter (total_table, Abundance_Bacteria > 1)
total_table = total_table[unsplit(table(total_table$Order), total_table$Order) > 10, ]
NMDS_Location = select(total_table, Location, Sample, Location_Sample, Order, RA)
NMDS_Location$Order[NMDS_Location$Order==""] <- "Unclassified Bacteria"
NMDS_Location2 = pivot_wider(NMDS_Location, names_from = Order, values_from = RA, values_fill = 0, values_fn = sum)

NMDS_Location2 = select(NMDS_Location2, -Location, -Location_Sample)
NMDS_Location2Bact = na.omit(NMDS_Location2)

Taxonomy =  read.delim("global.taxonomy.table.txt", header=TRUE, sep = ",")
Bacterial_Taxonomy =  filter(Taxonomy, str_detect(Taxonomy, "Bacteria"))
Bacterial_Taxonomy =Bacterial_Taxonomy%>% gather(key ="Sample", value = "Abundance", -Taxonomy)
Bacterial_Taxonomy =   Bacterial_Taxonomy%>% separate ( "Taxonomy", into = c ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") 
Bcounts = Bacterial_Taxonomy %>% select(Sample, Abundance)
Bcounts =  aggregate( data=Bcounts, Abundance~Sample, sum)
colnames(Bcounts) <- c("Sample", "Bacterial_Counts")
Fungal_Taxonomy =  filter(Taxonomy, str_detect(Taxonomy, c("Fungi")))
Fungal_Taxonomy =Fungal_Taxonomy%>% gather(key ="Sample", value = "Abundance", -Taxonomy)
Fungal_Taxonomy =   Fungal_Taxonomy%>% separate ( "Taxonomy", into = c ("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") 
Fcounts = Fungal_Taxonomy %>% select(Sample, Abundance)
Fcounts =  aggregate( data=Fcounts, Abundance~Sample, sum)
metadata <- read_delim("metadata.csv", ",", escape_double = FALSE, trim_ws = TRUE)
colnames(Fcounts) <- c("Sample", "Fungal_Counts")
total_counts = full_join(Fcounts, Bcounts, by="Sample")
total_table = full_join(total_counts, metadata, by="Sample")
total_table = filter (total_table, Sequencing %in% c ("Metagenome")) %>% mutate(., `RFBR`=Fungal_Counts/Bacterial_Counts) %>% filter(., Fungal_Counts>1)
fungal_SSU=select(total_table, Sample, `RFBR`)

joined = full_join(NMDS_Location2ARG, NMDS_Location2Bact,  by="Sample")%>%
full_join(., fungal_SSU, by="Sample")
joined = na.omit(joined)
joined = column_to_rownames(joined, var="Sample")
joined = select(joined, -`Unclassified Bacteria`)

joined[joined==0] <-NA
list_cor = list()
list_p  = list()
log_ratio <- function(data)
{
  # Compute the logarithmus
  log_data <- log(data)
  # Calculate exponential function of column-wise mean values for finite log transformed data
  gm <- exp(mean(log_data[is.finite(log_data)]))
  # Compute the logarithmus
  log_gm <- log(gm)
  # Take the difference of both log-transformed datasets
  data <- log_data - log_gm
  # Return the new OTU table
  return(data)
}

library(igraph)
for (i in unique(1:185)){
  j=i+1
  for (g in unique(j:186)){
    print(paste(i,g))
    file1=joined[c(i,g)]
    file2=na.omit(file1)
    file2[[1]]=log_ratio(file2[[1]])
    file2[[2]]=log_ratio(file2[[2]])
    
    if (length(row.names(file2))>25){
      print(colnames(file1))
      correl= cor.test(file2[[1]],file2[[2]],  method = "spearman")
      list_cor[paste0(colnames(file1[1]),"_", colnames(file1[2]))]=correl$estimate
      list_p[paste0(colnames(file1[1]),"_", colnames(file1[2]))]=correl$p.value
      
    }}}


lista_cor=as.data.frame(unlist(list_cor))
lista_p=as.data.frame(unlist(list_p))
common_paramterers=cbind(lista_cor, lista_p)     
common_paramterers$p=p.adjust(common_paramterers$`unlist(list_p)`, method = "BH")
common_paramterers$col=rownames(common_paramterers)
common_paramterers = filter(common_paramterers, p<0.05)
common_paramterers$corr=common_paramterers$`unlist(list_cor)`

df = select(common_paramterers, col, corr)
df = separate(df, col = "col", into = c("row", "col"), sep="_")
dflist= unique(df$row, df$col)
name <- unique(dflist)
bb=as.data.frame(name)
code = dcast(Bacterial_Taxonomy, Phylum+Order~., value.var = "Abundance", sum)%>%
select(., Phylum, Order) %>% mutate(name=Order)%>%
  filter(., Order%in%(unique(bb$name)))
length(unique(code$Order))

tbl_nodes <- data.frame(name) %>%
  full_join(., code, by="name")%>%filter(., name!="NA"&name!="")
tbl_nodes$Group = tbl_nodes$Phylum
for (i in 1:20){
  if(paste(tbl_nodes[i,4])=="NA"){
    print("ok")
    tbl_nodes[i,4]= tbl_nodes[i,1]  }
}
  
tbl_nodes=select(tbl_nodes,name, Group)%>% filter(., name%in%unique(name))
length(unique(na$name))

df$Correlation=df$corr
df$Correlation[df$Correlation>0]<-"Positive"
df$Correlation[df$Correlation<0]<-"Negative"

g <- graph_from_data_frame(df, directed = FALSE)
weight.scale <- c(1, g$corr)^3
coords <- layout_nicely(g)
a = plot.igraph(g, edge.color=ifelse(df$corr > 0, "#227a37","#e31a1c"),    vertex.label.font=2,  
            vertex.label.cex=0.65,   vertex.label.color = "darkblue", 
            vertex.color="#dfe1f5", vertex.size=0.8, edge.width=E(g)$weight, layout = coords)

library(ggraph)
#Figure S4, Full Network
ggraph(g, layout = 'auto') + 
  geom_edge_link(aes(colour =Correlation)) + 
  geom_node_point(size=4) + theme_pubr() + geom_node_label(aes(label = name)) + coord_flip() +
  xlab("") + ylab("")  + theme( axis.text.y = element_text( size=25, face = "bold", colour="white")) + 
  theme(axis.text.x = element_text(size = 25, face = "bold", colour="white")) +
  theme(legend.text = element_text( size=25, face="bold", colour = "black")) + 
  theme( axis.title.y = element_text( size=25, face = "bold", colour="black")) + theme(axis.title.x = element_text(size = 25, face = "bold", colour="black")) +
  theme(legend.title =  element_text( size=25, face="bold", colour = "black")) + scale_colour_discrete(palette = "Set1", name="Correlation")
  


ddlist= filter(df, row%in%c("RFBR", "aph(3')", "blaTEM", "sul1", "blaOXA")|col%in%c("RFBR", "aph(3')", "blaTEM", "sul1", "blaOXA"))
labels=unique(ddlist$row, ddlist$col)
g <- graph_from_data_frame(ddlist, directed = F)
weight.scale <- c(1, g$corr)
coords <- layout.auto(g)
#Figure 4, Isolated Network
plotl = plot.igraph(g, edge.color= "skyblue",    vertex.label.font=3,  
                vertex.label.cex=1,   vertex.label.color = "black", 
                vertex.color="#dfe1f5", vertex.size=1, edge.width=E(g)$weight, layout = coords)




# Relative Abunance Figure S1
b=ggplot(data = gathered, aes (x=Gene, y=log10(RA)))  +  geom_jitter(color="royalblue4")   + theme_pubr(legend = c("right")) + 
  theme( axis.text.y = element_text( size=25, face = "bold", colour="black")) +
  theme(axis.text.x = element_text(size = 19, face = "bold", colour="black", angle = 90)) +
  theme(legend.text = element_text( size=25, face="bold", colour = "black")) + 
  theme( axis.title.y = element_text( size=25, face = "bold", colour="black")) + 
  theme(axis.title.x = element_text(size = 25, face = "bold", colour="black")) +
  theme(legend.title =  element_text( size=25, face="bold", colour = "black"))  + ylab(" log10 ARG hits/Bacterial SSU") + xlab("Gene") #+ stat_summary(fun = median, fun.min= median, fun.max = median, geom = "crossbar", width=0.3, color="darkred")  

b 






# Figure S3 Phyla & Fungal/Bactaerial SSU ratio
Taxonomy =  read.delim("global.taxonomy.table.txt", header=TRUE, sep = ",")
Bacterial_Taxonomy =  filter(Taxonomy, str_detect(Taxonomy, "Bacteria"))
Bacterial_Taxonomy =Bacterial_Taxonomy%>% gather(key ="Sample", value = "Abundance", -Taxonomy)
Bacterial_Taxonomy =   Bacterial_Taxonomy%>% separate ( "Taxonomy", into = c ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") 
Bcounts = Bacterial_Taxonomy %>% select(Sample, Abundance)
Bcounts =  aggregate( data=Bcounts, Abundance~Sample, sum)
colnames(Bcounts) <- c("Sample", "Bacterial_Counts")
Fungi_Taxonomy =  filter(Taxonomy, str_detect(Taxonomy, c("Fungi")))
Fungi_Taxonomy =Fungi_Taxonomy%>% gather(key ="Sample", value = "Abundance", -Taxonomy)
Fungi_Taxonomy =   Fungi_Taxonomy%>% separate ( "Taxonomy", into = c ("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") 
FungiCounts = Fungi_Taxonomy %>% select(Sample, Abundance)
FungiCounts =  aggregate( data=FungiCounts, Abundance~Sample, sum)
Bacteria_Taxonomy =  filter(Taxonomy, str_detect(Taxonomy, c("Bacteria")))
Bacteria_Taxonomy =Bacteria_Taxonomy%>% gather(key ="Sample", value = "Abundance", -Taxonomy)
Bacteria_Taxonomy =   Bacteria_Taxonomy%>% separate ( "Taxonomy", into = c ( "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") 
BacteriaCounts = Bacteria_Taxonomy %>% select(Sample, Abundance, Phylum)
BacteriaCounts =  aggregate( data=BacteriaCounts, Abundance~Phylum+Sample, sum)
metadata <- read_delim("metadata.csv", ",", escape_double = FALSE, trim_ws = TRUE)
colnames(FungiCounts) <- c("Sample", "FCounts")
total_counts = full_join(FungiCounts, Bcounts, by="Sample")
total_table = full_join(total_counts, BacteriaCounts, by="Sample")
total_table = full_join(total_table, metadata, by="Sample")
total_table = filter (total_table, Sequencing %in% c ("Metagenome"))
total_table = mutate(total_table, RA = ((Abundance)/Bacterial_Counts ))
total_table = mutate(total_table, RAFungi = ((FCounts)/Bacterial_Counts ))
total_table = filter (total_table, FCounts > 1)
total_table = filter (total_table, Abundance > 1)
total_table$RAFungi = log10(total_table$RAFungi )
total_table$RA = log10(total_table$RA )

#Correlations of Bacterial Phyla with Fungal/Bacterial SSS Ratio Figure S2
plot_fp = ggscatter( data = filter(total_table, Phylum!=""), y="RA", x="RAFungi", add = "reg.line",
                     conf.int = FALSE, size = 0.5, color = "royalblue4", add.params = list (color = "royalblue4", fill = "#c8d4e3" ))  + 
  stat_cor(method = "spearman", size=4, colour="black")  + 
  xlab("log10 Fungal SSU hits/Total Bacterial SSU rRNA hits") + 
  ylab("log10 Relative Abundance (SSU%) ")  + theme( axis.text.y = element_text( size=25, face = "bold", colour="black")) + theme(axis.text.x = element_text(size = 25, face = "bold", colour="black")) +
  theme(legend.text = element_text( size=25, face="bold", colour = "black")) + theme( axis.title.y = element_text( size=25, face = "bold", colour="black")) + theme(axis.title.x = element_text(size = 25, face = "bold", colour="black")) +
  theme(legend.title =  element_text( size=25, face="bold", colour = "black")) +
  facet_wrap(~Phylum, ncol = 5) + theme(strip.text = element_text(size=15,  face="bold", colour = "white")) + theme( strip.background = element_rect( fill="skyblue4"))
plot_fp


