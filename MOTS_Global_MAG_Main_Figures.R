library("magrittr")
library("ggplot2")
library("forcats")
library("dplyr")
library("phyloseq")
library("ggpubr")
library("tidyr")
library("Rmisc")
library("microViz")
library("readxl")
library("tidyverse")
library("conflicted")
library("ggExtra")
library("cowplot")
library("gridExtra")
library("ggforce")
library("vegan")
library('rstatix')
library("ggpmisc")
library("data.table")
library("tibble")
library(grid)

My_Theme = theme(axis.title.x = element_text(face="bold",size=24),
                 axis.text.x = element_text(colour = "black", size=22), 
                 axis.text.y = element_text(colour = "black", size=22),
                 axis.title.y = element_text(face="bold",size=24),
                 plot.title = element_text(size = 24),
                 legend.title =element_text(face="bold",size = 14),
                 legend.text = element_text(size = 14),
                 legend.key.size = unit(1, "cm"),
                 strip.text.x = element_text(size=22, face="bold"),
                 strip.text.y = element_text(size=22, face="bold"),
                 panel.border = element_rect(fill = NA, colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank())

library(microViz)
brewerPlus <- distinct_palette()


library(viridis)
library(ggExtra)

#import data from Table S1 with all bins 
All_bins <- read.csv("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/raw_data/MOTS_bin_base_data.csv",fill = TRUE, header = TRUE, sep = ",")

All_bins$Watermass = factor(All_bins$Watermass, levels = c("NW", "STW", "Front", "SAW", "Deep"))

#making baseplot
Baseplot <- ggplot(All_bins, aes(x=Completeness, y=Contamination, color = Watermass)) + 
  geom_point(aes(size=Bin_Size.Mbp.)) + 
  scale_color_viridis(discrete=TRUE)+
  theme_light()+
  My_Theme+
  xlab("Completeness") + 
  ylab("Contamination")+
  labs(colour = "Water Mass", size = "Bin Size (Mbp)")+
  guides(color = guide_legend(override.aes = list(size = 5) ) )
Baseplot

#add marginal density plots
density_plot <- ggMarginal(Baseplot, type = "histogram", binwidth=1)

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/comp_cont_binsize.pdf",width=7,height=5) # Open a new pdf file

#view
density_plot

dev.off()
#Save as PDF

mean(All_bins$Contamination)
range(All_bins$Contamination)
mean(All_bins$Completeness)
range(All_bins$Completeness)

#facetted
baseplot_facet<-Baseplot +facet_grid(.~Watermass)
baseplot_facet

#load data
library(ggforce)

dodge <- position_dodge(width=0.5)  # move dots .01 to the left and right to avoid overlap


qual_plot<-ggplot(All_bins, aes(x=Completeness, y=Contamination, colour=Watermass, group=Watermass)) +
  geom_point(aes(size=Bin_Size.Mbp.),alpha = 0.5)+
  scale_color_viridis(discrete=TRUE)+
  labs(colour = "Water Mass", size = "Bin Size (Mbp)", tag = "A")+
  theme_light()+
  My_Theme+
  theme(legend.position="none")+ 
  guides(colour = guide_legend(override.aes = list(size=5, alpha = 0.5))) +
  xlim(40, 110) +
  scale_y_continuous(breaks=c(0, 4, 8, 12))
qual_plot 


library(cowplot)
#Save the legend 
conflicts_prefer(cowplot::get_legend)
legend <- get_legend(qual_plot)

#Remove the legend from qual_plot
qual_plot <- qual_plot + theme(legend.position="none")      



library(ggridges)
library(ggpubr)

binsize_ridge<-ggplot(All_bins, aes(x=Bin_Size.Mbp., y=Watermass, fill=Watermass))+
  geom_density_ridges(jittered_points = TRUE,alpha = 0.5) +
  scale_y_discrete(drop=FALSE)+
  scale_color_viridis(discrete=TRUE)+
  theme_light()+
  My_Theme+
  xlab("Bin Size (Mbp)") + ylab("Water Mass")+ 
  expand_limits(y= c(1, length(levels(All_bins$Watermass)) + 1.5))+
  theme(legend.position="none")      

binsize_ridge

bincont_ridge<-ggplot(All_bins, aes(x=Contamination, y=Watermass, fill=Watermass))+
  geom_density_ridges(jittered_points = TRUE,alpha = 0.5) +
  scale_y_discrete(drop=FALSE)+
  scale_color_viridis(discrete=TRUE)+
  theme_light()+
  My_Theme+
  labs(tag = "B") +
  xlab("Contamination") + ylab("Water Mass")+ 
  expand_limits(y= c(1, length(levels(All_bins$Watermass)) + 1.5))+
  rotate()+theme(legend.position="none")


bincont_ridge

bincomp_ridge<-ggplot(All_bins, aes(x=Completeness, y=Watermass, fill=Watermass))+
  geom_density_ridges(jittered_points = TRUE,alpha = 0.5) +
  scale_y_discrete(drop=FALSE)+
  scale_color_viridis(discrete=TRUE)+
  theme_light()+
  theme(axis.title.x = element_text(face="bold",size=20),
        axis.text.x = element_text(colour = "black", size=22), 
        axis.text.y = element_text(colour = "black", size=18),
        axis.title.y = element_text(face="bold",size=24),
        plot.title = element_text(size = 24),
        legend.title =element_text(face="bold",size = 14),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1, "cm"),
        strip.text.x = element_text(size=22, face="bold"),
        strip.text.y = element_text(size=22, face="bold"),
        panel.border = element_rect(fill = NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  xlab("Completeness") + ylab("Water Mass")+ 
  labs(tags = "C") +
  expand_limits(y= c(1.5, length(levels(All_bins$Watermass)) + 1.5))+
  theme(legend.position="none") +
  xlim(40,110)

bincomp_ridge

blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
  )
library(gridExtra)
pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/bin_quality.pdf", width = 12, height = 10) # Open a new pdf file
grid.arrange(bincomp_ridge, legend, qual_plot, bincont_ridge, 
             ncol=2, nrow=2, widths=c(4, 2.5), heights=c(3, 4))



dev.off()

raw_read_count <- read.csv("/Users/yugibeast/Library/CloudStorage/OneDrive-UniversityofOtago/ARF/Manuscripts/MOTS_Global_MAG/raw_read_count_watermass.csv", header = TRUE)

#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

conflicts_prefer(plyr::rename)
df2 <- data_summary(raw_read_count , varname="reads",
                    groupnames = c("Watermas"))

# Convert dose to a factor variable
#df2$dose=as.factor(df2$dose)
head(df2)

df2$reads <- as.numeric(df2$reads)
df2$sd <- as.numeric(df2$sd)

ggplot(df2, aes(x = factor(Watermas, level=c('Neritic','STW','Front','SAW','Deep')), y = reads)) +
  geom_bar(stat = "summary", fun = "mean") +
  geom_errorbar((aes(ymin=reads-sd(reads), ymax=reads+sd(reads)))) +
  labs(x = "Watermass", y = "Average read count") +
  My_Theme

## ANOVA


library(dplyr)
conflicts_prefer(summarise::dplyr)
group_by(raw_read_count, Watermas) %>%
  summarise(
    count = n(),
    mean = mean(weight, na.rm = TRUE),
    sd = sd(weight, na.rm = TRUE)
  )

library("ggpubr")
ggboxplot(my_data, x = "group", y = "weight", 
          color = "group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("ctrl", "trt1", "trt2"),
          ylab = "Weight", xlab = "Treatment")

# Compute the analysis of variance
res.aov <- aov(reads ~ Watermas, data = raw_read_count)
# Summary of the analysis
summary(res.aov)

kruskal.test(reads ~ Watermas, data = raw_read_count)

kruskal.test(reads ~ Watermas, data = raw_read_count)

library(car)
leveneTest(reads ~ Watermas, data = raw_read_count)


geom_errorbar(aes(ymin=len-sd, ymax=len+sd), width=.2,
              position=position_dodge(.9))

##Phyloseq
#OTU table
mots.otu.csv <- read.csv(file="/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/read_recruitment_to_bins/import/mots_otu.csv", header = TRUE, na.strings = "unclassified", fileEncoding="UTF-8-BOM", check.names = FALSE, row.names = 1)
mots.otu.matrix <- as.matrix(mots.otu.csv)

mots.otu <- otu_table(mots.otu.matrix, taxa_are_rows = TRUE)


#Taxa table
mots.taxa.csv <- read.csv(file="/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/read_recruitment_to_bins/import/taxa_table.csv", header = TRUE, na.strings = "unclassified", fileEncoding="UTF-8-BOM", check.names = FALSE, row.names = 1)



# #####
# #To allow all ranks to be shown I recoded the taxa_table file. This is how it was done. 
# #update taxonomy
# temp_meta_file <- read.csv(file="/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/bin_taxonomy.csv", header = TRUE, na.strings = "unclassified", fileEncoding="UTF-8-BOM", check.names = FALSE, row.names = 1)
# 
# library(stringr)
# #split 'player' column using '_' as the separator
# temp_meta_file[c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')] <- str_split_fixed(temp_meta_file$Classification, ';', 7)
# 
# temp_meta_file[, c(2:8)] <- as.data.frame(apply(temp_meta_file[,c(c(2:8))], 2,
#                                          function(temp_meta_file) {
#                                            gsub(".__", "", temp_meta_file)}
# ))
# 
# temp_meta_file<-temp_meta_file[ -c(1) ]
# 
# #remove taxonomy columns
# mots.taxa.csv <- mots.taxa.csv[ -c(1:5) ]
# 
# 
# #merge taxa_table
# new_taxa_table<-merge(mots.taxa.csv,temp_meta_file,by=0)
# library(tibble)
# new_taxa_table<-column_to_rownames(new_taxa_table, var = "Row.names")
# 
# #reorder
# new_taxa_table<-new_taxa_table %>% 
#    # dplyr::relocate(disp) %>% ## simply make disp the first column
#    relocate(c("Domain":"Species"))
# 
# 
# #save as csv and check file for errors
# write.csv(new_taxa_table, "/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/read_recruitment_to_bins/import/taxa_table.csv", row.names = T)




mots.taxa.matrix <- as.matrix(mots.taxa.csv)

mots.taxa <- tax_table(mots.taxa.matrix)


#Meta-data
mots.sample.meta <- read.csv(file="/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/read_recruitment_to_bins/import/mots_meta.csv", header = TRUE, fileEncoding="UTF-8-BOM", check.names = FALSE, row.names = 1)
mots.meta.df <- as.data.frame(mots.sample.meta)

mots.meta <- sample_data(mots.meta.df)


#Fix samples names to link w/ OTU

sample_names(mots.meta) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59", "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70", "71", "72", "73", "74", "75", "76", "77", "78", "79", "80", "81", "82", "83") # don't do this manually 
# sample_names(mots.meta) <-as.character(seq(1:83)) 


#Create phyloseq object
mots.phylo <- phyloseq(mots.otu, mots.taxa, mots.meta)


abundance <- psmelt(mots.phylo)

abundance_bin_count <- subset(abundance, select = c(OTU, Watermass, Phylum))

conflicts_prefer(dplyr::count)
abundance_bin_count <- dplyr::rename(count(abundance_bin_count, OTU, Watermass, Phylum), Freq = n)

#reorder Phyla by highest number of bins and replot
abundance_bin_count$Phylum<-factor(abundance_bin_count$Phylum,levels =c( "Proteobacteria","Bacteroidota","Thermoplasmatota","Verrucomicrobiota","Actinobacteriota","Cyanobacteria","Chloroflexota","Crenarchaeota","Planctomycetota","SAR324","Gemmatimonadota","Acidobacteriota","Nitrospinota","Myxococcota","Binatota","Dadabacteria","Fibrobacterota","Marinisomatota") )

abundance_bin_count$Watermass<-factor(abundance_bin_count$Watermass,levels =c("Neritic", "STW", "Front", "SAW", "Deep") )

abundance_bin_count_stacked_plot <- ggplot(abundance_bin_count, aes(x=Watermass, y = Freq, fill=Phylum)) +
  geom_bar(stat = "summary", color = "black") +
  scale_fill_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  labs(y = "Absolute abundance", x = "Water Mass", tag = "E") +
  theme(axis.title.x = element_text(face="bold",size=24),
        axis.text.x = element_text(colour = "black", size=22, angle = 45, hjust = 1), 
        axis.text.y = element_text(colour = "black", size=22),
        axis.title.y = element_text(face="bold",size=24),
        plot.title = element_text(size = 24),
        legend.title =element_text(face="bold",size = 14),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1, "cm"),
        strip.text.x = element_text(size=22, face="bold"),
        strip.text.y = element_text(size=22, face="bold"),
        panel.border = element_rect(fill = NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  guides(fill=guide_legend(ncol=2)) +
  theme(legend.position = "none")


abundance_bin_count_stacked_plot


#reorder Phyla by highest number of bins and replot
All_bins$Phylum<-factor(All_bins$Phylum,levels =c( "Proteobacteria","Bacteroidota","Thermoplasmatota","Verrucomicrobiota","Actinobacteriota","Cyanobacteria","Chloroflexota","Crenarchaeota","Planctomycetota","SAR324","Gemmatimonadota","Acidobacteriota","Nitrospinota","Myxococcota","Binatota","Dadabacteria","Fibrobacterota","Marinisomatota") )
All_bins$Watermass<-factor(All_bins$Watermass,levels =c("Neritic", "STW", "Front", "SAW", "Deep") )


taxanomic_breakdown <- ggplot(All_bins, aes(x=Watermass, fill=Phylum)) +
  geom_bar(stat="count", color = "black") +
  scale_fill_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  labs(y = "Absolute abundance", x = "Water Mass", tag = "E") +
  theme(axis.title.x = element_text(face="bold",size=24),
        axis.text.x = element_text(colour = "black", size=22, angle = 45, hjust = 1), 
        axis.text.y = element_text(colour = "black", size=22),
        axis.title.y = element_text(face="bold",size=24),
        plot.title = element_text(size = 24),
        legend.title =element_text(face="bold",size = 14),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1, "cm"),
        strip.text.x = element_text(size=22, face="bold"),
        strip.text.y = element_text(size=22, face="bold"),
        panel.border = element_rect(fill = NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  guides(fill=guide_legend(ncol=2)) +
  theme(legend.position = "none")


taxanomic_breakdown


  

abundance$Phylum<-factor(abundance$Phylum,levels =c( "Proteobacteria","Bacteroidota","Thermoplasmatota","Verrucomicrobiota","Actinobacteriota","Cyanobacteria","Chloroflexota","Crenarchaeota","Planctomycetota","SAR324","Gemmatimonadota","Acidobacteriota","Nitrospinota","Myxococcota","Binatota","Dadabacteria","Fibrobacterota","Marinisomatota"))

abundance$Watermass<-factor(abundance$Watermass,levels =c("Neritic","STW","Front","SAW","Deep"))

abundace_taxonomic_breakdown <- ggplot(abundance, aes(x=Watermass, fill=Phylum)) +
  geom_bar(stat="count", color = "black") +
  scale_fill_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#C3B4D1", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  labs(y = "Absolute abundance", x = "Water Mass") +
  My_Theme +
  guides(fill=guide_legend(ncol=2)) +
  theme(legend.position = "none")

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/for_phylum_legend.pdf",width=7,height=5) # Open a new pdf file

#view
abundace_taxonomic_breakdown

dev.off()

conflicts_prefer(ggpubr::get_legend)
legend <- get_legend(abundace_taxonomic_breakdown)



#75F5FD

ord <- ordinate(mots.phylo, "NMDS", "bray")
p1 = plot_ordination(mots.phylo, ord, type="taxa", color="Phylum")
p1 <- p1 +
  geom_point(size = 3) +
  My_Theme +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#C3B4D1", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey"))

p1

p1$data$NMDS1
p1$d

p2 = plot_ordination(mots.phylo, ord, type="sample", color="Watermass")
p2 <- p2 +
  geom_point(size = 3) +
  My_Theme #+
  scale_color_manual(values = c("Neritic" = "#0C0881", "STW" = "#7316A2", "Front" = "#BD5077", "SAW" = "#EA9953", "Deep" = "#F2F958"))

p2

nmds_taxa <- as.data.frame(p1$data)

ggplot(nmds_taxa, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(color = Phylum)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  theme(axis.title.x = element_text(face="bold",size=24),
        axis.text.x = element_text(colour = "black", size=22), 
        axis.text.y = element_text(colour = "black", size=22),
        axis.title.y = element_text(face="bold",size=24),
        plot.title = element_text(size = 24),
        legend.title =element_text(face="bold",size = 14),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1, "cm"),
        strip.text.x = element_text(size=22, face="bold"),
        strip.text.y = element_text(size=22, face="bold"),
        panel.border = element_rect(fill = NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  labs(x = "NMDS1", colour = "Phylum.y", y = "NMDS2", shape = "Type")

###16S Phyloseq object
Munida_OG

###FROM MERGE_SAMPLES
Munida_OG



MOG = Munida_OG
MOG = prune_taxa(taxa_sums(Munida_OG) > 0, Munida_OG)

S16S_MAGS_Sample_data <- read.csv("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/raw_data/16S_MAGS_Sample_Names.csv")

MAGS = S16S_MAGS_Sample_data$X16S_Samples
MAGS <- data.frame(MAGS)

sample_data(MOG)$MAGsamples <- get_variable(MOG, "SampleID") %in% MAGS

mergedMOG = merge_samples(MOG, "SampleID")
SD = merge_samples(sample_data(MOG), "SampleID")
print(SD[, "SampleID"])

S16_Matrix_merged <- distance(mergedMOG, method = "bray")

S16_Matrix_merged <- as.matrix(S16_Matrix_merged)

#S16_Matrix_merged$name <- row.names(S16_Matrix_merged)



##My ownn stuff

Munida_melted <- psmelt(Munida_OG)

Munida_melted <- subset(Munida_melted, select = c(Abundance, Phylum, SampleID, Water_mass_PCA, Season))

library(data.table)
setDT(Munida_melted)
Munida_melted1  = Munida_melted [ , .(counted = sum(Abundance)), by = .(Phylum, SampleID, Water_mass_PCA, Season)]

Munida_16S_matrix <- distance(Munida_OG, method = "bray", type = "Sample")
Munida_16S_matrix <- as.matrix(Munida_16S_matrix)

mots.phylo_matrix <- distance(mots.phylo, method = "bray", type = "Sample")
mots.phylo_matrix <- as.matrix(mots.phylo_matrix)
mots.phylo_matrix <- read.csv("~/Downloads/mots_phylo_matrix.csv", check.names=FALSE)
rownames(mots.phylo_matrix) <- mots.phylo_matrix$X
mots.phylo_matrix$X <- NULL

T.fm <- format(mots.phylo_matrix)
T.fm[row(T.fm) < col(T.fm)] <- ""
print(T.fm, quote=F)

T.ffm <- format(Munida_16S_matrix)
T.ffm[row(T.ffm) < col(T.ffm)] <- ""
print(T.fm, quote=F)



T.fm <- mots.phylo_matrix[order(row.names(x = mots.phylo_matrix)), order(colnames(x = mots.phylo_matrix))]
T.ffm <- Munida_16S_matrix[order(row.names(x = Munida_16S_matrix)), order(colnames(x = Munida_16S_matrix))]

# Remove rows/column 16S
T.ffm <- data.frame(T.ffm, check.names = FALSE)
T.ffm <- T.ffm[ , c("1A0315500", "1A0415", "1A0714", "1A1214", "1B0315", "1B0614", "1B0614500", "6A0315", "6A0415", "6B0714", "7B0614", "MunidaAug2015ST1a0m", "MunidaAug2015ST1a500m", "MunidaAug2015ST2a", "MunidaAug2015ST3a", "MunidaAug2015ST4a", "MunidaAug2015ST5a", "MunidaAug2015ST6a", "MunidaAug2015ST7a", "MunidaAug2015ST8b", "MunidaAug2016ST3a", "MunidaAug2016ST5a", "MunidaAug2016ST8b0m", "MunidaAug2016ST8b500m", "MunidaDec2015ST1a0m", "MunidaDec2015ST1a500m", "MunidaDec2015ST2a", "MunidaDec2015ST3a", "MunidaDec2015ST4a", "MunidaDec2015ST5a", "MunidaDec2015ST6a", "MunidaDec2015ST7a", "MunidaFeb2017ST2a", "MunidaFeb2017ST3a", "MunidaFeb2017ST8a0m", "MunidaFeb2017ST8a500m", "MunidaFeb22017ST4a", "MunidaFeb22017ST6a", "MunidaFeb22017ST8a0m", "MunidaFeb22017ST8b500m", "MunidaJan2016ST1a0m", "MunidaJan2016ST1a500m", "MunidaJan2016ST7b","MunidaMar2016ST1a0m", "MunidaMar2016ST1a500m", "MunidaMar2016ST2b", "MunidaMar2016ST3a", "MunidaMar2016ST4a", "MunidaMar2016ST5a", "MunidaMar2016ST6b", "MunidaMar2016ST7a", "MunidaMar2016ST8a", "MunidaMar2017ST2a", "MunidaMar2017ST3a", "MunidaMar2017ST8a0m", "MunidaNov2015ST1a0m", "MunidaNov2015ST1a500m", "MunidaNov2015ST7a") ]

# getting rows  
rows <- c("1A0315500", "1A0415", "1A0714", "1A1214", "1B0315", "1B0614", "1B0614500", "6A0315", "6A0415", "6B0714", "7B0614", "MunidaAug2015ST1a0m", "MunidaAug2015ST1a500m", "MunidaAug2015ST2a", "MunidaAug2015ST3a", "MunidaAug2015ST4a", "MunidaAug2015ST5a", "MunidaAug2015ST6a", "MunidaAug2015ST7a", "MunidaAug2015ST8b", "MunidaAug2016ST3a", "MunidaAug2016ST5a", "MunidaAug2016ST8b0m", "MunidaAug2016ST8b500m", "MunidaDec2015ST1a0m", "MunidaDec2015ST1a500m", "MunidaDec2015ST2a", "MunidaDec2015ST3a", "MunidaDec2015ST4a", "MunidaDec2015ST5a", "MunidaDec2015ST6a", "MunidaDec2015ST7a", "MunidaFeb2017ST2a", "MunidaFeb2017ST3a", "MunidaFeb2017ST8a0m", "MunidaFeb2017ST8a500m", "MunidaFeb22017ST4a", "MunidaFeb22017ST6a", "MunidaFeb22017ST8a0m", "MunidaFeb22017ST8b500m", "MunidaJan2016ST1a0m", "MunidaJan2016ST1a500m", "MunidaJan2016ST7b","MunidaMar2016ST1a0m", "MunidaMar2016ST1a500m", "MunidaMar2016ST2b", "MunidaMar2016ST3a", "MunidaMar2016ST4a", "MunidaMar2016ST5a", "MunidaMar2016ST6b", "MunidaMar2016ST7a", "MunidaMar2016ST8a", "MunidaMar2017ST2a", "MunidaMar2017ST3a", "MunidaMar2017ST8a0m", "MunidaNov2015ST1a0m", "MunidaNov2015ST1a500m", "MunidaNov2015ST7a") 

# extracting data frame rows 
T.ffm <- T.ffm [rownames(T.ffm) %in% rows, ] 


# Remove rows/column MAGS

T.fm <- data.frame(T.fm, check.names = FALSE)
T.fm <- T.fm[ , c("1A_03/15_500", "1A_04/15_0", "1A_07/14_0", "1A_12/14_0", "1B_03/15_0", "1B_06/14_0", "1B_06/14_500", "6A_03/15_0", "6A_04/15_0", "6B_07/14_0", "7B_06/14_0", "Munida_Aug2015_ST1a_0m", "Munida_Aug2015_ST1a_500m", "Munida_Aug2015_ST2a", "Munida_Aug2015_ST3a", "Munida_Aug2015_ST4a", "Munida_Aug2015_ST5a", "Munida_Aug2015_ST6a", "Munida_Aug2015_ST7a", "Munida_Aug2015_ST8b", "Munida_Aug2016_ST3a", "Munida_Aug2016_ST5", "Munida_Aug2016_ST8a_0m", "Munida_Aug2016_ST8b_500m", "Munida_Dec2015_ST1a_0m", "Munida_Dec2015_ST1a_500m", "Munida_Dec2015_ST2a", "Munida_Dec2015_ST3a", "Munida_Dec2015_ST4a", "Munida_Dec2015_ST5a", "Munida_Dec2015_ST6a", "Munida_Dec2015_ST7a", "Munida_Feb2017_ST2a", "Munida_Feb2017_ST3", "Munida_Feb2017_ST8a_0m", "Munida_Feb2017_ST8a_500m", "Munida_Feb22017_ST4a", "Munida_Feb22017_ST6", "Munida_Feb22017_ST8a_0m", "Munida_Feb22017_ST8b_500m", "Munida_Jan2016_ST1a_0m", "Munida_Jan2016_ST1a_500m", "Munida_Jan2016_ST7b", "Munida_Mar2016_ST1a_0m", "Munida_Mar2016_ST1a_500m", "Munida_Mar2016_ST2b", "Munida_Mar2016_ST3a", "Munida_Mar2016_ST4a", "Munida_Mar2016_ST5a", "Munida_Mar2016_ST6b", "Munida_Mar2016_ST7a", "Munida_Mar2016_ST8a", "Munida_Mar2017_ST2a", "Munida_Mar2017_ST3", "Munida_Mar2017_ST8a_0m", "Munida_Nov2015_ST1a_0m", "Munida_Nov2015_ST1a_500m", "Munida_Nov2015_ST7a") ]


rows <- c("1A_03/15_500", "1A_04/15_0", "1A_07/14_0", "1A_12/14_0", "1B_03/15_0", "1B_06/14_0", "1B_06/14_500", "6A_03/15_0", "6A_04/15_0", "6B_07/14_0", "7B_06/14_0", "Munida_Aug2015_ST1a_0m", "Munida_Aug2015_ST1a_500m", "Munida_Aug2015_ST2a", "Munida_Aug2015_ST3a", "Munida_Aug2015_ST4a", "Munida_Aug2015_ST5a", "Munida_Aug2015_ST6a", "Munida_Aug2015_ST7a", "Munida_Aug2015_ST8b", "Munida_Aug2016_ST3a", "Munida_Aug2016_ST5", "Munida_Aug2016_ST8a_0m", "Munida_Aug2016_ST8b_500m", "Munida_Dec2015_ST1a_0m", "Munida_Dec2015_ST1a_500m", "Munida_Dec2015_ST2a", "Munida_Dec2015_ST3a", "Munida_Dec2015_ST4a", "Munida_Dec2015_ST5a", "Munida_Dec2015_ST6a", "Munida_Dec2015_ST7a", "Munida_Feb2017_ST2a", "Munida_Feb2017_ST3", "Munida_Feb2017_ST8a_0m", "Munida_Feb2017_ST8a_500m", "Munida_Feb22017_ST4a", "Munida_Feb22017_ST6", "Munida_Feb22017_ST8a_0m", "Munida_Feb22017_ST8b_500m", "Munida_Jan2016_ST1a_0m", "Munida_Jan2016_ST1a_500m", "Munida_Jan2016_ST7b", "Munida_Mar2016_ST1a_0m", "Munida_Mar2016_ST1a_500m", "Munida_Mar2016_ST2b", "Munida_Mar2016_ST3a", "Munida_Mar2016_ST4a", "Munida_Mar2016_ST5a", "Munida_Mar2016_ST6b", "Munida_Mar2016_ST7a", "Munida_Mar2016_ST8a", "Munida_Mar2017_ST2a", "Munida_Mar2017_ST3", "Munida_Mar2017_ST8a_0m", "Munida_Nov2015_ST1a_0m", "Munida_Nov2015_ST1a_500m", "Munida_Nov2015_ST7a")
T.fm <- T.fm [rownames(T.fm) %in% rows, ] 


T.fm <- as.matrix(T.fm)
T.ffm <- as.matrix(T.ffm)

mantel(T.fm, T.ffm, method = "spearman", permutations = 9999, na.rm = TRUE)

m_com <- T.ffm


tot <- rowSums(m_com)
m_com <- m_com[tot > 0, ]

m_com <- na.omit(m_com)

set.seed(123)
nmds = metaMDS(m_com, distance = "bray")
nmds

data.scores_16S = as.data.frame(scores(nmds))

data.scores_16S$SampleID <- row.names(data.scores_16S)

S16_sample_data <- sample_data(Munida_OG)
S16_sample_data <- data.frame(S16_sample_data)

data.scores_16S <- merge(data.scores_16S, S16_sample_data, by = "SampleID")


data.scores_16S$Water_mass_PCA<-factor(data.scores_16S$Water_mass_PCA,levels =c("NW", "STW", "FRONT", "SAW", "DEEP") )


NMDS_16S <- ggplot(data.scores_16S, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(color = Water_mass_PCA)) +
  theme(axis.title.x = element_text(face="bold",size=24),
        axis.text.x = element_text(colour = "black", size=22), 
        axis.text.y = element_text(colour = "black", size=22),
        axis.title.y = element_text(face="bold",size=24),
        plot.title = element_text(size = 24),
        legend.title =element_text(face="bold",size = 14),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1, "cm"),
        strip.text.x = element_text(size=22, face="bold"),
        strip.text.y = element_text(size=22, face="bold"),
        panel.border = element_rect(fill = NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_color_viridis(discrete=TRUE, option = "C") +
  theme(legend.position="none") +
  labs(x = "16S NMDS1", colour = "Water_mass_PCA", y = "16S NMDS2", shape = "Type") +
  xlim(-0.06, 0.06)

NMDS_16S

#MAG
m_com <- T.fm


tot <- rowSums(m_com)
m_com <- m_com[tot > 0, ]

m_com <- na.omit(m_com)

set.seed(123)
nmds = metaMDS(m_com, distance = "bray")
nmds

data.scores_MAG = as.data.frame(scores(nmds))

data.scores_MAG$MOTS_ID <- row.names(data.scores_MAG)

MAG_sample_data <- sample_data(mots.phylo)
MAG_sample_data <- data.frame(MAG_sample_data)

data.scores_MAG <- merge(data.scores_MAG, MAG_sample_data, by = "MOTS_ID")


data.scores_MAG$Water_mass_PCA<-factor(data.scores_MAG$Water_mass_PCA,levels =c("NW", "STW", "FRONT", "SAW", "DEEP") )


NMDS_MAG <- ggplot(data.scores_MAG, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(color = Watermass)) +
  theme(axis.title.x = element_text(face="bold",size=24),
        axis.text.x = element_text(colour = "black", size=22), 
        axis.text.y = element_text(colour = "black", size=22),
        axis.title.y = element_text(face="bold",size=24),
        plot.title = element_text(size = 24),
        legend.title =element_text(face="bold",size = 14),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1, "cm"),
        strip.text.x = element_text(size=22, face="bold"),
        strip.text.y = element_text(size=22, face="bold"),
        panel.border = element_rect(fill = NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_color_viridis(discrete=TRUE, option = "C") +
  theme(legend.position="none") +
  labs(x = "MAG NMDS1", colour = "Watermass", y = "MAG NMDS2", shape = "Type")

NMDS_MAG


pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/16S_vs_MAG_16S.pdf",width=6,height=9) # Open a new pdf file

grid.arrange(NMDS_16S, NMDS_MAG)

dev.off()
#mantel


###mantel
df <- t(mots.otu.csv)

temp <- mots.meta.df$`Temp (oC)`

dist.abun <- vegdist(df, method = "bray")
dist.temp <- dist(temp, method = "euclidean")

abund_temp = mantel(dist.abun, dist.temp, method = "spearman", permutations = 9999, na.rm = TRUE)
abund_temp

#Melt phyloseq object

rank_df <- mots.taxa.csv[ -c(8:15) ]

#rank_df[is.na(rank_df)] <- Unclassified

#calculate unique levels by rank
rank_summary <-as.data.frame(apply(rank_df, 2, function(x) length(unique(x))))
rank_summary

library(tibble)
rank_summary<-rank_df %>% 
  #group_by(Watermass) %>%
  summarise_all(n_distinct)%>% 
  t()%>%
  data.frame()%>%
  rownames_to_column(var = "Rank")%>%
  rename_with(.cols = 2, ~"Count")

rank_summary

#calculate the number of BINS classified at each rank
rank_summary$Count <-as.numeric(rank_summary$Count)
rank_summary$Unclassified_count <-colSums(is.na(rank_df) | rank_df == "")
rank_summary$Classified_count<-1027-rank_summary$Unclassified_count
rank_summary$percent_bins_classified<- ((rank_summary$Classified_count)/nrow(rank_df))*100

#reorder Ranks
rank_summary$Rank <- factor(rank_summary$Rank, levels=c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Watermass"))

rank_summary <- rank_summary[-8,]

library(viridis)
percent_classified_MAGs <- ggplot(rank_summary,aes(x= Rank, y= Count, size=percent_bins_classified, label=Count)) + 
  geom_point()+
  scale_size_continuous(limits = c(27, 100),breaks=c(30, 90), labels=expression("<30%",">90%"))+
  labs(y="Number of Distinct Taxa", x="Taxonomic Rank", size = "Percent of classified MAGS", tag = "E")+
  #My_Theme +
  theme(axis.title.x = element_text(face="bold",size=24),
        axis.text.x = element_text(angle=45, colour = "black", vjust=1, hjust = 1, size=22), 
        axis.text.y = element_text(colour = "black", size=25),
        axis.title.y = element_text(face="bold",size=24),
        plot.title = element_text(size = 24),
        legend.title =element_text(face="bold",size = 14),
        legend.text = element_text(size = 14),
        legend.position="right",
        legend.key.size = unit(1, "cm"),
        panel.background = element_blank(),
        strip.text.x = element_text(size=22, color="black"),
        strip.text.y = element_text(size=22, face="bold"),
        panel.border = element_rect(fill = NA, colour = "black")) +
  theme(legend.position.inside=c(0.3, 0.7)) +
  scale_y_continuous(limits = c(0, 150),breaks=c(0,25,50,75, 100, 125, 150))+geom_text(aes(label = Count),size = 6, color="black", vjust = -1.5)


pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/percent_classified_MAGs.pdf",width=7,height=5) # Open a new pdf file

#view
percent_classified_MAGs

dev.off()

blank <- grid.rect(gp=gpar(col="white"))
#Save as PDF

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/Figure_1.pdf", width = 16, height = 18) # Open a new pdf file
grid.arrange(qual_plot, bincont_ridge, bincomp_ridge, percent_classified_MAGs, taxanomic_breakdown,
             ncol=2, nrow=3, widths=c(3, 3), heights=c(4, 4, 4))

dev.off()

##Figure 2

mots.data <- read.csv("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/raw_data/mots_mag_fnct.csv")


#Rename column lablled n to Counts
names(mots.data)[names(mots.data) == "n"] <- "Counts"
names(mots.data)

#Remove count data that are showing up as NA
mots.data2=na.omit(mots.data)

#Order phylum ascending
mots.data2$Phylum=factor(mots.data2$Phylum, levels=c("Dadabacteria","Binatota","Fibrobacterota","Marinisomatota","Myxococcota","Nitrospinota","Acidobacteriota","Gemmatimonadota","Planctomycetota","SAR324",                                            "Crenarchaeota","Chloroflexota","Cyanobacteria","Verrucomicrobiota","Actinobacteriota","Bacteroidota",
                                                     "Thermoplasmatota","Proteobacteria"))


library(tidyverse)
library(grid)
library(cowplot)



p2 <- ggplot(data = mots.data2, aes(x = Phylum,y = Function, size=Counts))+
  geom_point(aes(fill=Completeness),pch=21)+
  scale_fill_viridis_c(option="plasma", limits = c(75, 90), oob = scales::squish)+ scale_size_continuous(limits=c(0,400), range=c(0.5,30), breaks=c(1, 50, 100, 200 ,300))+
  facet_grid("Category", scales = "free_y", space = "free",labeller=labeller(Category=label_wrap_gen(10)))+
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid = element_line(color = "grey35"),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.title =element_text(face="bold",size = 24),
    legend.text = element_text(size = 22),
    axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1, face="bold",size=22),
    axis.text.y = element_text(size=22),
    axis.title.x = element_text(color = "black",size=22),
    axis.title.y = element_text(color = "black",size=22),
    strip.background = element_rect(fill = "grey20"),
    strip.text = element_text(colour = "white",size=5,lineheight=0.01),
    legend.position="none"
  )+
  labs(y="Gene", x="Phylum")



p2=p2+theme(panel.spacing =unit(.05, "lines"),
            panel.border = element_rect(color = "black", fill = NA, size = .5), 
            strip.background = element_rect(color = "white", size = .5))

p2


g2 <- ggplot_gtable(ggplot_build(p2))
stripr <- which(grepl('strip-r', g2$layout$name))
fills <- c("cadetblue3","aquamarine2","darkgoldenrod1","lightcoral","thistle3","rosybrown4","lightsteelblue3","yellow")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g2$grobs[[i]]$grobs[[1]]$childrenOrder))
  g2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

grid.draw(g2)



func_col <- c("Alt e acceptor"= "cadetblue3",
              "Alt e donor"="aquamarine2",
              "Carbon fixation"="darkgoldenrod1",
              "Nitrogen cycle"="lightcoral",
              "Phototrophy"="thistle3",
              "Respiration"="rosybrown4",
              "Sulfur cycle"="lightsteelblue3",
              "Trace gas metabolism
"="yellow")


legend_func<-ggplot(mots.data2, aes(x=Category, fill=Phylum)) +
  geom_bar(stat="count")+
  scale_fill_manual(values = func_col)+
  My_Theme+
  theme(legend.position="bottom",
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.5, "cm"))+
  guides(fill=guide_legend(nrow=2))


#Save the legend  
func_legend <- get_legend(legend_func)
grid.newpage()
grid.draw(func_legend)

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/diamond_plot.pdf", width = 8, height = 10) # Open a new pdf file
grid.arrange(g2, func_legend, ncol=1, nrow=2, widths=c(8), heights=c(10, 1))
dev.off()

##Figure 3
##NMDS plot using total functional data DRAM
##plot with watermass as colours
MISC1 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/1_genome_summaries/metabolism_summary.xlsx", sheet = "MISC")

MISC1 <-MISC1 %>%
  mutate(MISC1,"function_description" = "Misc", .after = gene_description)

carbon_utilization1 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/1_genome_summaries/metabolism_summary.xlsx", sheet = "carbon utilization")

carbon_utilization1 <-carbon_utilization1 %>%
  mutate(carbon_utilization1,"function_description" = "carbon utilization", .after = gene_description)

Transporters1 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/1_genome_summaries/metabolism_summary.xlsx", sheet = "Transporters")

Transporters1 <-Transporters1 %>%
  mutate(Transporters1,"function_description" = "Transporters", .after = gene_description)

Energy1 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/1_genome_summaries/metabolism_summary.xlsx", sheet = "Energy")

Energy1 <-Energy1 %>%
  mutate(Energy1,"function_description" = "Energy", .after = gene_description)

organic_nitrogen1 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/1_genome_summaries/metabolism_summary.xlsx", sheet = "Organic Nitrogen")

organic_nitrogen1 <-organic_nitrogen1 %>%
  mutate(organic_nitrogen1,"function_description" = "Organic Nitrogen", .after = gene_description)

Woodcroft1 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/1_genome_summaries/metabolism_summary.xlsx", sheet = "carbon utilization (Woodcroft)")

Woodcroft1 <-Woodcroft1 %>%
  mutate(Woodcroft1,"function_description" = "carbon utilization (Woodcroft)", .after = gene_description)

rRNA1 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/1_genome_summaries/metabolism_summary.xlsx", sheet = "rRNA")

rRNA1 <-rRNA1 %>%
  mutate(rRNA1,"function_description" = "rRNA", .after = gene_description)

tRNA1 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/1_genome_summaries/metabolism_summary.xlsx", sheet = "tRNA")

tRNA1 <-tRNA1 %>%
  mutate(tRNA1,"function_description" = "tRNA", .after = gene_description)



##2_genome_summaries
MISC2 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/2_genome_summaries/metabolism_summary.xlsx", sheet = "MISC")

MISC2 <-MISC2 %>%
  mutate(MISC2,"function_description" = "Misc", .after = gene_description)

MISC2 <- subset(MISC2, select = -c(gene_id, gene_description, function_description, module, header, subheader))

carbon_utilization2 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/2_genome_summaries/metabolism_summary.xlsx", sheet = "carbon utilization")

carbon_utilization2 <-carbon_utilization2 %>%
  mutate(carbon_utilization2,"function_description" = "carbon utilization", .after = gene_description)

carbon_utilization2 <- subset(carbon_utilization2, select = -c(gene_id, gene_description, function_description, module, header, subheader))

Transporters2 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/2_genome_summaries/metabolism_summary.xlsx", sheet = "Transporters")

Transporters2 <-Transporters2 %>%
  mutate(Transporters2,"function_description" = "Transporters", .after = gene_description)

Transporters2 <- subset(Transporters2, select = -c(gene_id, gene_description, function_description, module, header, subheader))

Energy2 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/2_genome_summaries/metabolism_summary.xlsx", sheet = "Energy")

Energy2 <-Energy2 %>%
  mutate(Energy2,"function_description" = "Energy", .after = gene_description)

Energy2 <- subset(Energy2, select = -c(gene_id, gene_description, function_description, module, header, subheader))

organic_nitrogen2 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/2_genome_summaries/metabolism_summary.xlsx", sheet = "Organic Nitrogen")

organic_nitrogen2 <-organic_nitrogen2 %>%
  mutate(organic_nitrogen2,"function_description" = "Organic Nitrogen", .after = gene_description)

organic_nitrogen2 <- subset(organic_nitrogen2, select = -c(gene_id, gene_description, function_description, module, header, subheader))

Woodcroft2 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/2_genome_summaries/metabolism_summary.xlsx", sheet = "carbon utilization (Woodcroft)")

Woodcroft2 <-Woodcroft2 %>%
  mutate(Woodcroft2,"function_description" = "carbon utilization (Woodcroft)", .after = gene_description)

Woodcroft2 <- subset(Woodcroft2, select = -c(gene_id, gene_description, function_description, module, header, subheader))

rRNA2 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/2_genome_summaries/metabolism_summary.xlsx", sheet = "rRNA")

rRNA2 <-rRNA2 %>%
  mutate(rRNA2,"function_description" = "rRNA", .after = gene_description)

rRNA2 <- subset(rRNA2, select = -c(gene_id, gene_description, function_description, module, header, subheader))

tRNA2 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/2_genome_summaries/metabolism_summary.xlsx", sheet = "tRNA")

tRNA2 <-tRNA2 %>%
  mutate(tRNA2,"function_description" = "tRNA", .after = gene_description)

tRNA2 <- subset(tRNA2, select = -c(gene_id, gene_description, function_description, module, header, subheader))

##3_genome_summaries
MISC3 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/3_genome_summaries/metabolism_summary.xlsx", sheet = "MISC")

MISC3 <-MISC3 %>%
  mutate(MISC3,"function_description" = "Misc", .after = gene_description)

MISC3 <- subset(MISC3, select = -c(gene_id, gene_description, function_description, module, header, subheader))

carbon_utilization3 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/3_genome_summaries/metabolism_summary.xlsx", sheet = "carbon utilization")

carbon_utilization3 <-carbon_utilization3 %>%
  mutate(carbon_utilization3,"function_description" = "carbon utilization", .after = gene_description)

carbon_utilization3 <- subset(carbon_utilization3, select = -c(gene_id, gene_description, function_description, module, header, subheader))

Transporters3 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/3_genome_summaries/metabolism_summary.xlsx", sheet = "Transporters")

Transporters3 <-Transporters3 %>%
  mutate(Transporters3,"function_description" = "Transporters", .after = gene_description)

Transporters3 <- subset(Transporters3, select = -c(gene_id, gene_description, function_description, module, header, subheader))

Energy3 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/3_genome_summaries/metabolism_summary.xlsx", sheet = "Energy")

Energy3 <-Energy3 %>%
  mutate(Energy3,"function_description" = "Energy", .after = gene_description)

Energy3 <- subset(Energy3, select = -c(gene_id, gene_description, function_description, module, header, subheader))

organic_nitrogen3 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/3_genome_summaries/metabolism_summary.xlsx", sheet = "Organic Nitrogen")

organic_nitrogen3 <-organic_nitrogen3 %>%
  mutate(organic_nitrogen3,"function_description" = "Organic Nitrogen", .after = gene_description)

organic_nitrogen3 <- subset(organic_nitrogen3, select = -c(gene_id, gene_description, function_description, module, header, subheader))

Woodcroft3 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/3_genome_summaries/metabolism_summary.xlsx", sheet = "carbon utilization (Woodcroft)")

Woodcroft3 <-Woodcroft3 %>%
  mutate(Woodcroft3,"function_description" = "carbon utilization (Woodcroft)", .after = gene_description)

Woodcroft3 <- subset(Woodcroft3, select = -c(gene_id, gene_description, function_description, module, header, subheader))

rRNA3 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/3_genome_summaries/metabolism_summary.xlsx", sheet = "rRNA")

rRNA3 <-rRNA3 %>%
  mutate(rRNA3,"function_description" = "rRNA", .after = gene_description)

rRNA3 <- subset(rRNA3, select = -c(gene_id, gene_description, function_description, module, header, subheader))

tRNA3 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/3_genome_summaries/metabolism_summary.xlsx", sheet = "tRNA")

tRNA3 <-tRNA3 %>%
  mutate(tRNA3,"function_description" = "tRNA", .after = gene_description)

tRNA3 <- subset(tRNA3, select = -c(gene_id, gene_description, function_description, module, header, subheader))

##4_genome_summaries
MISC4 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/4_genome_summaries/metabolism_summary.xlsx", sheet = "MISC")

MISC4 <-MISC4 %>%
  mutate(MISC4,"function_description" = "Misc", .after = gene_description)

MISC4 <- subset(MISC4, select = -c(gene_id, gene_description, function_description, module, header, subheader))

carbon_utilization4 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/4_genome_summaries/metabolism_summary.xlsx", sheet = "carbon utilization")

carbon_utilization4 <-carbon_utilization4 %>%
  mutate(carbon_utilization4,"function_description" = "carbon utilization", .after = gene_description)

carbon_utilization4 <- subset(carbon_utilization4, select = -c(gene_id, gene_description, function_description, module, header, subheader))

Transporters4 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/4_genome_summaries/metabolism_summary.xlsx", sheet = "Transporters")

Transporters4 <-Transporters4 %>%
  mutate(Transporters4,"function_description" = "Transporters", .after = gene_description)

Transporters4 <- subset(Transporters4, select = -c(gene_id, gene_description, function_description, module, header, subheader))

Energy4 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/4_genome_summaries/metabolism_summary.xlsx", sheet = "Energy")

Energy4 <-Energy4 %>%
  mutate(Energy4,"function_description" = "Energy", .after = gene_description)

Energy4 <- subset(Energy4, select = -c(gene_id, gene_description, function_description, module, header, subheader))

organic_nitrogen4 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/4_genome_summaries/metabolism_summary.xlsx", sheet = "Organic Nitrogen")

organic_nitrogen4 <-organic_nitrogen4 %>%
  mutate(organic_nitrogen4,"function_description" = "Organic Nitrogen", .after = gene_description)

organic_nitrogen4 <- subset(organic_nitrogen4, select = -c(gene_id, gene_description, function_description, module, header, subheader))

Woodcroft4 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/4_genome_summaries/metabolism_summary.xlsx", sheet = "carbon utilization (Woodcroft)")

Woodcroft4 <-Woodcroft4 %>%
  mutate(Woodcroft4,"function_description" = "carbon utilization (Woodcroft)", .after = gene_description)

Woodcroft4 <- subset(Woodcroft4, select = -c(gene_id, gene_description, function_description, module, header, subheader))

rRNA4 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/4_genome_summaries/metabolism_summary.xlsx", sheet = "rRNA")

rRNA4 <-rRNA4 %>%
  mutate(rRNA4,"function_description" = "rRNA", .after = gene_description)

rRNA4 <- subset(rRNA4, select = -c(gene_id, gene_description, function_description, module, header, subheader))

tRNA4 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/4_genome_summaries/metabolism_summary.xlsx", sheet = "tRNA")

tRNA4 <-tRNA4 %>%
  mutate(tRNA4,"function_description" = "tRNA", .after = gene_description)

tRNA4 <- subset(tRNA4, select = -c(gene_id, gene_description, function_description, module, header, subheader))


##5_genome_summaries
MISC5 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/5_genome_summaries/metabolism_summary.xlsx", sheet = "MISC")

MISC5 <-MISC5 %>%
  mutate(MISC5,"function_description" = "Misc", .after = gene_description)

MISC5 <- subset(MISC5, select = -c(gene_id, gene_description, function_description, module, header, subheader))

carbon_utilization5 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/5_genome_summaries/metabolism_summary.xlsx", sheet = "carbon utilization")

carbon_utilization5 <-carbon_utilization5 %>%
  mutate(carbon_utilization5,"function_description" = "carbon utilization", .after = gene_description)

carbon_utilization5 <- subset(carbon_utilization5, select = -c(gene_id, gene_description, function_description, module, header, subheader))

Transporters5 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/5_genome_summaries/metabolism_summary.xlsx", sheet = "Transporters")

Transporters5 <-Transporters5 %>%
  mutate(Transporters5,"function_description" = "Transporters", .after = gene_description)

Transporters5 <- subset(Transporters5, select = -c(gene_id, gene_description, function_description, module, header, subheader))

Energy5 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/5_genome_summaries/metabolism_summary.xlsx", sheet = "Energy")

Energy5 <-Energy5 %>%
  mutate(Energy5,"function_description" = "Energy", .after = gene_description)

Energy5 <- subset(Energy5, select = -c(gene_id, gene_description, function_description, module, header, subheader))

organic_nitrogen5 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/5_genome_summaries/metabolism_summary.xlsx", sheet = "Organic Nitrogen")

organic_nitrogen5 <-organic_nitrogen5 %>%
  mutate(organic_nitrogen5,"function_description" = "Organic Nitrogen", .after = gene_description)

organic_nitrogen5 <- subset(organic_nitrogen5, select = -c(gene_id, gene_description, function_description, module, header, subheader))

Woodcroft5 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/5_genome_summaries/metabolism_summary.xlsx", sheet = "carbon utilization (Woodcroft)")

Woodcroft5 <-Woodcroft5 %>%
  mutate(Woodcroft5,"function_description" = "carbon utilization (Woodcroft)", .after = gene_description)

Woodcroft5 <- subset(Woodcroft5, select = -c(gene_id, gene_description, function_description, module, header, subheader))

rRNA5 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/5_genome_summaries/metabolism_summary.xlsx", sheet = "rRNA")

rRNA5 <-rRNA5 %>%
  mutate(rRNA5,"function_description" = "rRNA", .after = gene_description)

rRNA5 <- subset(rRNA5, select = -c(gene_id, gene_description, function_description, module, header, subheader))

tRNA5 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/5_genome_summaries/metabolism_summary.xlsx", sheet = "tRNA")

tRNA5 <-tRNA5 %>%
  mutate(tRNA5,"function_description" = "tRNA", .after = gene_description)

tRNA5 <- subset(tRNA5, select = -c(gene_id, gene_description, function_description, module, header, subheader))

##6_genome_summaries
MISC6 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/6_genome_summaries/metabolism_summary.xlsx", sheet = "MISC")

MISC6 <-MISC6 %>%
  mutate(MISC6,"function_description" = "Misc", .after = gene_description)

MISC6 <- subset(MISC6, select = -c(gene_id, gene_description, function_description, module, header, subheader))

carbon_utilization6 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/6_genome_summaries/metabolism_summary.xlsx", sheet = "carbon utilization")

carbon_utilization6 <-carbon_utilization6 %>%
  mutate(carbon_utilization6,"function_description" = "carbon utilization", .after = gene_description)

carbon_utilization6 <- subset(carbon_utilization6, select = -c(gene_id, gene_description, function_description, module, header, subheader))

Transporters6 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/6_genome_summaries/metabolism_summary.xlsx", sheet = "Transporters")

Transporters6 <-Transporters6 %>%
  mutate(Transporters6,"function_description" = "Transporters", .after = gene_description)

Transporters6 <- subset(Transporters6, select = -c(gene_id, gene_description, function_description, module, header, subheader))

Energy6 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/6_genome_summaries/metabolism_summary.xlsx", sheet = "Energy")

Energy6 <-Energy6 %>%
  mutate(Energy6,"function_description" = "Energy", .after = gene_description)

Energy6 <- subset(Energy6, select = -c(gene_id, gene_description, function_description, module, header, subheader))

organic_nitrogen6 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/6_genome_summaries/metabolism_summary.xlsx", sheet = "Organic Nitrogen")

organic_nitrogen6 <-organic_nitrogen6 %>%
  mutate(organic_nitrogen6,"function_description" = "Organic Nitrogen", .after = gene_description)

organic_nitrogen6 <- subset(organic_nitrogen6, select = -c(gene_id, gene_description, function_description, module, header, subheader))

Woodcroft6 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/6_genome_summaries/metabolism_summary.xlsx", sheet = "carbon utilization (Woodcroft)")

Woodcroft6 <-Woodcroft6 %>%
  mutate(Woodcroft6,"function_description" = "carbon utilization (Woodcroft)", .after = gene_description)

Woodcroft6 <- subset(Woodcroft6, select = -c(gene_id, gene_description, function_description, module, header, subheader))

rRNA6 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/6_genome_summaries/metabolism_summary.xlsx", sheet = "rRNA")

rRNA6 <-rRNA6 %>%
  mutate(rRNA6,"function_description" = "rRNA", .after = gene_description)

rRNA6 <- subset(rRNA6, select = -c(gene_id, gene_description, function_description, module, header, subheader))

tRNA6 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/6_genome_summaries/metabolism_summary.xlsx", sheet = "tRNA")

tRNA6 <-tRNA6 %>%
  mutate(tRNA6,"function_description" = "tRNA", .after = gene_description)

tRNA6 <- subset(tRNA6, select = -c(gene_id, gene_description, function_description, module, header, subheader))


##7_genome_summaries
MISC7 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/7_genome_summaries/metabolism_summary.xlsx", sheet = "MISC")

MISC7 <-MISC7 %>%
  mutate(MISC7,"function_description" = "Misc", .after = gene_description)

MISC7 <- subset(MISC7, select = -c(gene_id, gene_description, function_description, module, header, subheader))

carbon_utilization7 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/7_genome_summaries/metabolism_summary.xlsx", sheet = "carbon utilization")

carbon_utilization7 <-carbon_utilization7 %>%
  mutate(carbon_utilization7,"function_description" = "carbon utilization", .after = gene_description)

carbon_utilization7 <- subset(carbon_utilization7, select = -c(gene_id, gene_description, function_description, module, header, subheader))

Transporters7 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/7_genome_summaries/metabolism_summary.xlsx", sheet = "Transporters")

Transporters7 <-Transporters7 %>%
  mutate(Transporters7,"function_description" = "Transporters", .after = gene_description)

Transporters7 <- subset(Transporters7, select = -c(gene_id, gene_description, function_description, module, header, subheader))

Energy7 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/7_genome_summaries/metabolism_summary.xlsx", sheet = "Energy")

Energy7 <-Energy7 %>%
  mutate(Energy7,"function_description" = "Energy", .after = gene_description)

Energy7 <- subset(Energy7, select = -c(gene_id, gene_description, function_description, module, header, subheader))

organic_nitrogen7 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/7_genome_summaries/metabolism_summary.xlsx", sheet = "Organic Nitrogen")

organic_nitrogen7 <-organic_nitrogen7 %>%
  mutate(organic_nitrogen7,"function_description" = "Organic Nitrogen", .after = gene_description)

organic_nitrogen7 <- subset(organic_nitrogen7, select = -c(gene_id, gene_description, function_description, module, header, subheader))

Woodcroft7 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/7_genome_summaries/metabolism_summary.xlsx", sheet = "carbon utilization (Woodcroft)")

Woodcroft7 <-Woodcroft7 %>%
  mutate(Woodcroft7,"function_description" = "carbon utilization (Woodcroft)", .after = gene_description)

Woodcroft7 <- subset(Woodcroft7, select = -c(gene_id, gene_description, function_description, module, header, subheader))

rRNA7 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/7_genome_summaries/metabolism_summary.xlsx", sheet = "rRNA")

rRNA7 <-rRNA7 %>%
  mutate(rRNA7,"function_description" = "rRNA", .after = gene_description)

rRNA7 <- subset(rRNA7, select = -c(gene_id, gene_description, function_description, module, header, subheader))

tRNA7 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/7_genome_summaries/metabolism_summary.xlsx", sheet = "tRNA")

tRNA7 <-tRNA7 %>%
  mutate(tRNA7,"function_description" = "tRNA", .after = gene_description)

tRNA7 <- subset(tRNA7, select = -c(gene_id, gene_description, function_description, module, header, subheader))

##8_genome_summaries
MISC8 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/8_genome_summaries/metabolism_summary.xlsx", sheet = "MISC")

MISC8 <-MISC8 %>%
  mutate(MISC8,"function_description" = "Misc", .after = gene_description)

MISC8 <- subset(MISC8, select = -c(gene_id, gene_description, function_description, module, header, subheader))

carbon_utilization8 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/8_genome_summaries/metabolism_summary.xlsx", sheet = "carbon utilization")

carbon_utilization8 <-carbon_utilization8 %>%
  mutate(carbon_utilization8,"function_description" = "carbon utilization", .after = gene_description)

carbon_utilization8 <- subset(carbon_utilization8, select = -c(gene_id, gene_description, function_description, module, header, subheader))

Transporters8 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/8_genome_summaries/metabolism_summary.xlsx", sheet = "Transporters")

Transporters8 <-Transporters8 %>%
  mutate(Transporters8,"function_description" = "Transporters", .after = gene_description)

Transporters8 <- subset(Transporters8, select = -c(gene_id, gene_description, function_description, module, header, subheader))

Energy8 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/8_genome_summaries/metabolism_summary.xlsx", sheet = "Energy")

Energy8 <-Energy8 %>%
  mutate(Energy8,"function_description" = "Energy", .after = gene_description)

Energy8 <- subset(Energy8, select = -c(gene_id, gene_description, function_description, module, header, subheader))

organic_nitrogen8 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/8_genome_summaries/metabolism_summary.xlsx", sheet = "Organic Nitrogen")

organic_nitrogen8 <-organic_nitrogen8 %>%
  mutate(organic_nitrogen8,"function_description" = "Organic Nitrogen", .after = gene_description)

organic_nitrogen8 <- subset(organic_nitrogen8, select = -c(gene_id, gene_description, function_description, module, header, subheader))

Woodcroft8 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/8_genome_summaries/metabolism_summary.xlsx", sheet = "carbon utilization (Woodcroft)")

Woodcroft8 <-Woodcroft8 %>%
  mutate(Woodcroft8,"function_description" = "carbon utilization (Woodcroft)", .after = gene_description)

Woodcroft8 <- subset(Woodcroft8, select = -c(gene_id, gene_description, function_description, module, header, subheader))

rRNA8 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/8_genome_summaries/metabolism_summary.xlsx", sheet = "rRNA")

rRNA8 <-rRNA8 %>%
  mutate(rRNA8,"function_description" = "rRNA", .after = gene_description)

rRNA8 <- subset(rRNA8, select = -c(gene_id, gene_description, function_description, module, header, subheader))

tRNA8 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/8_genome_summaries/metabolism_summary.xlsx", sheet = "tRNA")

tRNA8 <-tRNA8 %>%
  mutate(tRNA8,"function_description" = "tRNA", .after = gene_description)

tRNA8 <- subset(tRNA8, select = -c(gene_id, gene_description, function_description, module, header, subheader))

##9_genome_summaries
MISC9 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/9_genome_summaries/metabolism_summary.xlsx", sheet = "MISC")

MISC9 <-MISC9 %>%
  mutate(MISC9,"function_description" = "Misc", .after = gene_description)

MISC9 <- subset(MISC9, select = -c(gene_id, gene_description, function_description, module, header, subheader))

carbon_utilization9 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/9_genome_summaries/metabolism_summary.xlsx", sheet = "carbon utilization")

carbon_utilization9 <-carbon_utilization9 %>%
  mutate(carbon_utilization9,"function_description" = "carbon utilization", .after = gene_description)

carbon_utilization9 <- subset(carbon_utilization9, select = -c(gene_id, gene_description, function_description, module, header, subheader))

Transporters9 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/9_genome_summaries/metabolism_summary.xlsx", sheet = "Transporters")

Transporters9 <-Transporters9 %>%
  mutate(Transporters9,"function_description" = "Transporters", .after = gene_description)

Transporters9 <- subset(Transporters9, select = -c(gene_id, gene_description, function_description, module, header, subheader))

Energy9 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/9_genome_summaries/metabolism_summary.xlsx", sheet = "Energy")

Energy9 <-Energy9 %>%
  mutate(Energy9,"function_description" = "Energy", .after = gene_description)

Energy9 <- subset(Energy9, select = -c(gene_id, gene_description, function_description, module, header, subheader))

organic_nitrogen9 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/9_genome_summaries/metabolism_summary.xlsx", sheet = "Organic Nitrogen")

organic_nitrogen9 <-organic_nitrogen9 %>%
  mutate(organic_nitrogen9,"function_description" = "Organic Nitrogen", .after = gene_description)

organic_nitrogen9 <- subset(organic_nitrogen9, select = -c(gene_id, gene_description, function_description, module, header, subheader))

Woodcroft9 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/9_genome_summaries/metabolism_summary.xlsx", sheet = "carbon utilization (Woodcroft)")

Woodcroft9 <-Woodcroft9 %>%
  mutate(Woodcroft9,"function_description" = "carbon utilization (Woodcroft)", .after = gene_description)

Woodcroft9 <- subset(Woodcroft9, select = -c(gene_id, gene_description, function_description, module, header, subheader))

rRNA9 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/9_genome_summaries/metabolism_summary.xlsx", sheet = "rRNA")

rRNA9 <-rRNA9 %>%
  mutate(rRNA9,"function_description" = "rRNA", .after = gene_description)

rRNA9 <- subset(rRNA9, select = -c(gene_id, gene_description, function_description, module, header, subheader))

tRNA9 <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/9_genome_summaries/metabolism_summary.xlsx", sheet = "tRNA")

tRNA9 <-tRNA9 %>%
  mutate(tRNA9,"function_description" = "tRNA", .after = gene_description)

tRNA9 <- subset(tRNA9, select = -c(gene_id, gene_description, function_description, module, header, subheader))


MISC_combined <- dplyr::bind_cols(MISC1,MISC2,MISC3,MISC4,MISC5,MISC6,MISC7,MISC8,MISC9)

carbon_utilization_combined <- dplyr::bind_cols(carbon_utilization1,carbon_utilization2,carbon_utilization3,carbon_utilization4,carbon_utilization5,carbon_utilization6,carbon_utilization7,carbon_utilization8,carbon_utilization9)

Transporters_combined <- dplyr::bind_cols(Transporters1,Transporters2,Transporters3,Transporters4,Transporters5,Transporters6,Transporters7,Transporters8,Transporters9)

Energy_combined <- dplyr::bind_cols(Energy1,Energy2,Energy3,Energy4,Energy5,Energy6,Energy7,Energy8,Energy9)

organic_nitrogen_combined <- dplyr::bind_cols(organic_nitrogen1,organic_nitrogen2,organic_nitrogen3,organic_nitrogen4,organic_nitrogen5,organic_nitrogen6,organic_nitrogen7,organic_nitrogen8,organic_nitrogen9)

Woodcroft_combined <- dplyr::bind_cols(Woodcroft1,Woodcroft2,Woodcroft3,Woodcroft4,Woodcroft5,Woodcroft6,Woodcroft7,Woodcroft8,Woodcroft9)

rRNA_combined <- dplyr::bind_cols(rRNA1,rRNA2,rRNA3,rRNA4,rRNA5,rRNA6,rRNA7,rRNA8,rRNA9)

#tRNA_combined <- dplyr::bind_cols(tRNA1,tRNA2,tRNA3,tRNA4,tRNA5,tRNA6,tRNA7,tRNA8,tRNA9)

library(vegan)
remove(metabolism_summary)

metabolism_summary<-dplyr::bind_rows(MISC_combined,carbon_utilization_combined,Transporters_combined,Energy_combined,organic_nitrogen_combined,Woodcroft_combined,rRNA_combined)

names(metabolism_summary)[1] <- "gene_id"
names(metabolism_summary)[2] <- "gene_description"
names(metabolism_summary)[3] <- "function_description"
names(metabolism_summary)[4] <- "module"
names(metabolism_summary)[5] <- "header"
names(metabolism_summary)[6] <- "subheader"


metabolism_summary$unique_id <- paste(metabolism_summary$function_description, metabolism_summary$gene_id,  metabolism_summary$module, sep="_")

#metabolism_summary<-dplyr::bind_rows(MISC1,carbon_utilization1,Transporters1,Energy1,Organic_Nitrogen1,Woodcroft1,rRNA1,tRNA1)



metabolism_summary$unique_id <- paste(metabolism_summary$function_description, metabolism_summary$gene_id,  metabolism_summary$module, sep="_")
metabolism_summary<-select(metabolism_summary,unique_id, everything())

metabolism_summary_for_dotplot <- gather(metabolism_summary, bin, count, "10_bin.1":"9_bin.9")

metabolism_summary_for_dotplot <- subset(metabolism_summary_for_dotplot, select = c(header, count))

conflicts_prefer(dplyr::filter)
metabolism_summary_for_dotplot <- filter(metabolism_summary_for_dotplot, count > 0)

metabolism_summary_for_dotplot$count <- as.numeric(metabolism_summary_for_dotplot$count)

library(data.table)
setDT(metabolism_summary_for_dotplot)
metabolism_summary_for_dotplot1  = metabolism_summary_for_dotplot [ , .(counted = sum(count)), by = .(header)]

##change header names

metabolism_summary_for_dotplot1$header <- str_replace(metabolism_summary_for_dotplot1$header,"central carbon", "Central carbon")
metabolism_summary_for_dotplot1$header <- str_replace(metabolism_summary_for_dotplot1$header,"Electron transport Chain", "Electron transport chain")
metabolism_summary_for_dotplot1$header <- str_replace(metabolism_summary_for_dotplot1$header,"Antibiotic Resistance", "Antibiotic resistance")
metabolism_summary_for_dotplot1$header <- str_replace(metabolism_summary_for_dotplot1$header,"pyruvate metabolism", "Pyruvate metabolism")
metabolism_summary_for_dotplot1$header <- str_replace(metabolism_summary_for_dotplot1$header,"Flagella Structure", "Flagella structure")
metabolism_summary_for_dotplot1$header <- str_replace(metabolism_summary_for_dotplot1$header,"hydrocarbon degradation", "Hydrocarbon degradation")
metabolism_summary_for_dotplot1$header <- str_replace(metabolism_summary_for_dotplot1$header,"sugar utilization (woodcroft)", "Sugar utilization (Woodcroft)")
metabolism_summary_for_dotplot1$header <- str_replace(metabolism_summary_for_dotplot1$header,"aerobic corrin ring synthesis", "Aerobic corrin ring synthesis")
metabolism_summary_for_dotplot1$header <- str_replace(metabolism_summary_for_dotplot1$header,"Metal Reduction", "Metal reduction")

DRAM_overview <- ggplot(metabolism_summary_for_dotplot1, aes(x = reorder(header, -counted), y = counted)) +
  geom_point(size = 3) +
  scale_y_continuous(trans='log10') +
  labs(x = "Gene category", y = "Total gene counts (log10)") +
  theme(axis.title.x = element_text(face="bold",size=24),
        axis.text.x = element_text(angle=90, colour = "black", vjust=0.5, hjust = 1, size=22), 
        axis.text.y = element_text(colour = "black", size=22),
        axis.title.y = element_text(face="bold",size=24),
        plot.title = element_text(size = 22),
        legend.title =element_text(face="bold",size = 24),
        legend.text = element_text(size = 22),
        legend.key.size = unit(1, "cm"),
        strip.text.x = element_text(size=18, face="bold"),
        strip.text.y = element_text(size=18, face="bold"),
        panel.border = element_rect(fill = NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

DRAM_overview
dev.off()

library(ragg)



pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/Figure_2.pdf", width = 14, height = 22) # Open a new pdf file
grid.arrange(g2, DRAM_overview,
             ncol=1, nrow=2, widths=c(10), heights=c(12, 6))

dev.off()

###DRAM_overview per watermass

metabolism_summary_for_dotplot <- gather(metabolism_summary, bin, count, "10_bin.1":"9_bin.9")
metabolism_summary_for_dotplot_watermass <- merge(metabolism_summary_for_dotplot, NMDS_metadata, by.x = "bin", by.y = "bin_ID")
metabolism_summary_for_dotplot_watermass <- subset(metabolism_summary_for_dotplot_watermass, select = -c(gene_description, gene_id, function_description, module, subheader, cluster, Sample, Phylum))

conflicts_prefer(dplyr::filter)
metabolism_summary_for_dotplot_watermass <- filter(metabolism_summary_for_dotplot_watermass, count > 0)

library(data.table)
setDT(metabolism_summary_for_dotplot_watermass)
metabolism_summary_for_dotplot_watermass1  = metabolism_summary_for_dotplot_watermass [ , .(counted = sum(count)), by = .(header, Watermass, bin)]


metabolism_summary_for_dotplot_watermass1$Watermass <- factor(metabolism_summary_for_dotplot_watermass1$Watermass,levels = c("Neritic", "STW", "Front", "SAW", "Deep"))

metabolism_summary_for_dotplot_watermass$Watermass <- factor(metabolism_summary_for_dotplot_watermass$Watermass,levels = c("Neritic", "STW", "Front", "SAW", "Deep"))

metabolism_summary_for_dotplot_watermass1$header <- str_replace(metabolism_summary_for_dotplot_watermass1$header,"central carbon", "Central carbon")
metabolism_summary_for_dotplot_watermass1$header <- str_replace(metabolism_summary_for_dotplot_watermass1$header,"Electron transport Chain", "Electron transport chain")
metabolism_summary_for_dotplot_watermass1$header <- str_replace(metabolism_summary_for_dotplot_watermass1$header,"Antibiotic Resistance", "Antibiotic resistance")
metabolism_summary_for_dotplot_watermass1$header <- str_replace(metabolism_summary_for_dotplot_watermass1$header,"pyruvate metabolism", "Pyruvate metabolism")
metabolism_summary_for_dotplot_watermass1$header <- str_replace(metabolism_summary_for_dotplot_watermass1$header,"Flagella Structure", "Flagella structure")
metabolism_summary_for_dotplot_watermass1$header <- str_replace(metabolism_summary_for_dotplot_watermass1$header,"hydrocarbon degradation", "Hydrocarbon degradation")
metabolism_summary_for_dotplot_watermass1$header <- str_replace(metabolism_summary_for_dotplot_watermass1$header,"sugar utilization (woodcroft)", "Sugar utilization (Woodcroft)")
metabolism_summary_for_dotplot_watermass1$header <- str_replace(metabolism_summary_for_dotplot_watermass1$header,"aerobic corrin ring synthesis", "Aerobic corrin ring synthesis")
metabolism_summary_for_dotplot_watermass1$header <- str_replace(metabolism_summary_for_dotplot_watermass1$header,"Metal Reduction", "Metal reduction")

DRAM_overview_per_watermass <- ggplot(metabolism_summary_for_dotplot_watermass1, aes(x = reorder(header, -counted), y = counted)) +
  scale_y_continuous(trans='log10') +
  geom_boxplot() +
  facet_wrap(~Watermass, ncol = 1)+
  labs(x = "Gene category", y = "Average gene counts / bin / water mass (Log)") +
  theme(axis.title.x = element_text(face="bold",size=24),
        axis.text.x = element_text(angle=90, colour = "black", vjust=1, hjust = 1, size=22), 
        axis.text.y = element_text(colour = "black", size=22),
        axis.title.y = element_text(face="bold",size=24),
        plot.title = element_text(size = 22),
        legend.title =element_text(face="bold",size = 24),
        legend.text = element_text(size = 22),
        legend.key.size = unit(1, "cm"),
        strip.text.x = element_text(size=18, face="bold"),
        strip.text.y = element_text(size=18, face="bold"),
        panel.border = element_rect(fill = NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())


DRAM_overview_per_watermass

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/DRAM_overview_by_watermass.pdf",width=10,height=18) # Open a new pdf file

#view
DRAM_overview_per_watermass

dev.off()

metabolism_summary_for_dotplot_watermass <- merge(metabolism_summary_for_dotplot, NMDS_metadata, by.x = "bin", by.y = "bin_ID")
metabolism_summary_for_dotplot_watermass <- subset(metabolism_summary_for_dotplot_watermass, select = -c(gene_description, gene_id, function_description, module, subheader, cluster, Sample, bin, Phylum))

conflicts_prefer(dplyr::filter)
metabolism_summary_for_dotplot_watermass <- filter(metabolism_summary_for_dotplot_watermass, count > 0)

library(data.table)
setDT(metabolism_summary_for_dotplot_watermass)
metabolism_summary_for_dotplot_watermass1  = metabolism_summary_for_dotplot_watermass [ , .(counted = sum(count)), by = .(header, Watermass)]


metabolism_summary_for_dotplot_watermass1$Watermass <- factor(metabolism_summary_for_dotplot_watermass1$Watermass,levels = c("Neritic", "STW", "Front", "SAW", "Deep"))

metabolism_summary_for_dotplot_watermass$Watermass <- factor(metabolism_summary_for_dotplot_watermass$Watermass,levels = c("Neritic", "STW", "Front", "SAW", "Deep"))

metabolism_summary_for_dotplot_watermass1$header <- str_replace(metabolism_summary_for_dotplot_watermass1$header,"central carbon", "Central carbon")
metabolism_summary_for_dotplot_watermass1$header <- str_replace(metabolism_summary_for_dotplot_watermass1$header,"Electron transport Chain", "Electron transport chain")
metabolism_summary_for_dotplot_watermass1$header <- str_replace(metabolism_summary_for_dotplot_watermass1$header,"Antibiotic Resistance", "Antibiotic resistance")
metabolism_summary_for_dotplot_watermass1$header <- str_replace(metabolism_summary_for_dotplot_watermass1$header,"pyruvate metabolism", "Pyruvate metabolism")
metabolism_summary_for_dotplot_watermass1$header <- str_replace(metabolism_summary_for_dotplot_watermass1$header,"Flagella Structure", "Flagella structure")
metabolism_summary_for_dotplot_watermass1$header <- str_replace(metabolism_summary_for_dotplot_watermass1$header,"hydrocarbon degradation", "Hydrocarbon degradation")
metabolism_summary_for_dotplot_watermass1$header <- str_replace(metabolism_summary_for_dotplot_watermass1$header,"sugar utilization (woodcroft)", "Sugar utilization (Woodcroft)")
metabolism_summary_for_dotplot_watermass1$header <- str_replace(metabolism_summary_for_dotplot_watermass1$header,"aerobic corrin ring synthesis", "Aerobic corrin ring synthesis")
metabolism_summary_for_dotplot_watermass1$header <- str_replace(metabolism_summary_for_dotplot_watermass1$header,"Metal Reduction", "Metal reduction")

for_watermass_bin <- subset(All_bins, select = c(Watermass))
conflicts_prefer(dplyr::count)
for_watermass_bin_count <- dplyr::rename(count(for_watermass_bin, Watermass), Freq = n)
for_watermass_bin_count$Watermass <- str_replace(for_watermass_bin_count$Watermass,"NW", "Neritic")

metabolism_summary_for_dotplot_watermass_norm <- merge(metabolism_summary_for_dotplot_watermass1, for_watermass_bin_count, by = "Watermass")
metabolism_summary_for_dotplot_watermass_norm$norm <- metabolism_summary_for_dotplot_watermass_norm$counted/metabolism_summary_for_dotplot_watermass_norm$Freq


DRAM_overview_between_watermass <- ggplot(metabolism_summary_for_dotplot_watermass1, aes(x = reorder(header, -counted), y = counted), group = counted, color = Watermass) +
  scale_color_viridis(discrete=TRUE, option = "C") +
  geom_point(size = 3, aes(color = factor(Watermass))) +
  labs(x = "Gene category", y = "Total gene counts per water mass (Log)", color = "Watermass") +
  theme(axis.title.x = element_text(face="bold",size=24),
        axis.text.x = element_text(angle=90, colour = "black", vjust=0.5, hjust = 1, size=22), 
        axis.text.y = element_text(colour = "black", size=22),
        axis.title.y = element_text(face="bold",size=24),
        plot.title = element_text(size = 22),
        legend.title =element_text(face="bold",size = 24),
        legend.text = element_text(size = 22),
        legend.key.size = unit(1, "cm"),
        strip.text.x = element_text(size=18, face="bold"),
        strip.text.y = element_text(size=18, face="bold"),
        panel.border = element_rect(fill = NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_y_continuous(trans='log10')

DRAM_overview_between_watermass_norm <- ggplot(metabolism_summary_for_dotplot_watermass_norm, aes(x = reorder(header, -norm), y = norm), group = norm, color = Watermass) +
  scale_color_viridis(discrete=TRUE, option = "C") +
  geom_point(size = 3, aes(color = factor(Watermass))) +
  labs(x = "Gene category", y = "Normalized total gene counts per water mass (Log)", color = "Watermass") +
  theme(axis.title.x = element_text(face="bold",size=24),
        axis.text.x = element_text(angle=90, colour = "black", vjust=0.5, hjust = 1, size=22), 
        axis.text.y = element_text(colour = "black", size=22),
        axis.title.y = element_text(face="bold",size=24),
        plot.title = element_text(size = 22),
        legend.title =element_text(face="bold",size = 24),
        legend.text = element_text(size = 22),
        legend.key.size = unit(1, "cm"),
        strip.text.x = element_text(size=18, face="bold"),
        strip.text.y = element_text(size=18, face="bold"),
        panel.border = element_rect(fill = NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_y_continuous(trans='log10')

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/DRAM_overview_between_watermass.pdf",width=10,height=18) # Open a new pdf file

#view
DRAM_overview_between_watermass

dev.off()

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/Normalized_DRAM_overview_between_watermass.pdf",width=10,height=18) # Open a new pdf file

#view
DRAM_overview_between_watermass_norm

dev.off()

library(vegan)
remove(metabolism_summary)

metabolism_summary<-dplyr::bind_rows(MISC_combined,carbon_utilization_combined,Transporters_combined,Energy_combined,organic_nitrogen_combined,Woodcroft_combined)

t_metabolism_summary <- t(metabolism_summary)

com = metabolism_summary[,7:ncol(metabolism_summary)]

com <- t(com)

m_com <- as.matrix(com)


tot <- rowSums(m_com)
m_com <- m_com[tot > 0, ]

m_com <- na.omit(m_com)

set.seed(123)
nmds = metaMDS(m_com, distance = "bray")
nmds

data.scores = as.data.frame(scores(nmds)$sites)

data.scores$bin_ID <- row.names(data.scores)

#add columns to data frame 
NMDS_metadata <- read.csv("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/read_recruitment_to_bins/Yugen_working_files/for_NMDS_metadata.csv")

with_season_metadata <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/read_recruitment_to_bins/Yugen_working_files/Nova_metadata_season.xlsx")

with_season_meta <- merge(x=NMDS_metadata, y=with_season_metadata, by.x='Sample', by.y='Sample')

mots.taxa.csv$bin_ID <- rownames(mots.taxa.csv)

mots_taxa <- mots.taxa.csv
mots_taxa <- subset(mots_taxa, select = c(bin_ID, Phylum, Watermass))

with_season_meta <- merge(x=with_season_meta, y = mots_taxa, by = "bin_ID")

with_season_meta <- subset(with_season_meta, select = c(bin_ID, Phylum.y, Watermass.y, Sample, Season, cluster, Year))

data.scores <- merge(x=with_season_meta, y=data.scores, by.x='bin_ID', by.y='bin_ID')

#data.scores <- cbind(with_season_meta, data.scores)

head(data.scores)

#data.scores <- subset(data.scores, select = c(bin_ID, Watermass.y, Phylum.y, Season, cluster.x, NMDS1, NMDS2))

rownames(data.scores ) <- NULL


#data.scores$Phylum.y<-factor(data.scores$Phylum.y,levels =c( "Proteobacteria","Bacteroidota","Thermoplasmatota","Verrucomicrobiota","Actinobacteriota","Cyanobacteria","Chloroflexota","Crenarchaeota","Planctomycetota","SAR324","Gemmatimonadota","Acidobacteriota","Nitrospinota","Myxococcota","Binatota","Dadabacteria","Fibrobacterota","Marinisomatota") )

data.scores$Phylum<-factor(data.scores$Phylum,levels =c( "Proteobacteria","Bacteroidota","Thermoplasmatota","Verrucomicrobiota","Actinobacteriota","Cyanobacteria","Chloroflexota","Crenarchaeota","Planctomycetota","SAR324","Gemmatimonadota","Acidobacteriota","Nitrospinota","Myxococcota","Binatota","Dadabacteria","Fibrobacterota","Marinisomatota") )

data.scores$cluster <- str_replace(data.scores$cluster, "clust1", "Cluster 1")
data.scores$cluster <- str_replace(data.scores$cluster, "clust2", "Cluster 2")
data.scores$cluster <- str_replace(data.scores$cluster, "clust3", "Cluster 3")

functional_NMDS_Phylum <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(color = Phylum.y)) +
  geom_mark_ellipse(aes(group = cluster, label = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  theme(axis.title.x = element_text(face="bold",size=24),
        axis.text.x = element_text(colour = "black", size=22), 
        axis.text.y = element_text(colour = "black", size=22),
        axis.title.y = element_text(face="bold",size=24),
        plot.title = element_text(size = 24),
        legend.title =element_text(face="bold",size = 14),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1, "cm"),
        strip.text.x = element_text(size=22, face="bold"),
        strip.text.y = element_text(size=22, face="bold"),
        panel.border = element_rect(fill = NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  labs(x = "NMDS1", colour = "Phylum.y", y = "NMDS2", shape = "Type", tag = "A") +
  ylim(-1.5, 1.1)


#view
functional_NMDS_Phylum

###create phyloseq for anosim and adonis






#OTU table
otu.function <- t(m_com)


mots.otu.matrix <- as.matrix(otu.function)

mots.function.otu <- otu_table(mots.otu.matrix, taxa_are_rows = TRUE)


#Taxa table
taxa.function <- metabolism_summary%>%
  select(gene_id,gene_description,function_description,module,header)

rownames(otu.function) <- taxa.function$gene_id

mots.taxa.matrix <- as.matrix(taxa.function)

mots.function.taxa <- tax_table(mots.taxa.matrix)


#Meta-data
#bins <- as.data.frame(rownames(m_com))
#names(bins)[names(bins) == 'rownames(m_com)'] <- 'bin_ID'
#to_combine_taxa <- tibble::rownames_to_column(mots.taxa.csv, "bin_ID")
meta.function <- subset(data.scores, select = -c(NMDS1, NMDS2))



meta.function.df <- as.data.frame(meta.function)

mots.function.meta <- sample_data(meta.function.df)
rownames(mots.function.meta) = mots.function.meta$bin_ID

#Fix samples names to link w/ OTU

#sample_names(mots.meta) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59", "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70", "71", "72", "73", "74", "75", "76", "77", "78", "79", "80", "81", "82", "83") # don't do this manually 
# sample_names(mots.meta) <-as.character(seq(1:83)) 


#Create phyloseq object
mots.function.phylo <- phyloseq(mots.function.otu, mots.function.taxa, mots.function.meta)

cluster_group = get_variable(mots.function.phylo, "cluster")
cluster_group_ano = anosim(phyloseq::distance(mots.function.phylo, "bray"), cluster_group)
cluster_group_ano$signif
cluster_group_ano$statistic

#anosim(m_com, data.scores$Phylum, distance = "bray", permutations = 9999)

#Create a data frame using your sample_data
df_ado_cluster = as(sample_data(mots.function.phylo), "data.frame")
#Calculate your Bray distance matrix
cluster_ado_horse = phyloseq::distance(mots.function.phylo, "bray")
#Perform your ADONIS test
cluster_ado_function_stats <- adonis2(cluster_ado_horse ~ cluster, df_ado_cluster)
cluster_ado_function_stats

#adonis2(m_com ~ cluster, data=data.scores, perm=999)

#NMDS2
for_bin_size <- read.csv("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/read_recruitment_to_bins/Yugen_working_files/for_bin_size_bin_metadata.csv")


#nmds_bin_size <- cbind(for_bin_size, data.scores_clust)
#for_bin_size %>%
#join(data.scores_clust, by = c("Bin_ID", "bin_ID"))

nmds_bin_size <- merge(x=for_bin_size, y=data.scores, by.x='Bin_ID', by.y='bin_ID')

NMDS2_binsize <- ggplot(nmds_bin_size, aes(x = NMDS2, y = Bin_Size.Mbp.)) + 
  geom_point(size = 4, aes(color = Phylum.y)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  theme(axis.title.x = element_text(face="bold",size=24),
        axis.text.x = element_text(colour = "black", size=22), 
        axis.text.y = element_text(colour = "black", size=22),
        axis.title.y = element_text(face="bold",size=24),
        plot.title = element_text(size = 24),
        legend.title =element_text(face="bold",size = 14),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1, "cm"),
        strip.text.x = element_text(size=22, face="bold"),
        strip.text.y = element_text(size=22, face="bold"),
        panel.border = element_rect(fill = NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position="none") +
  labs(x = "NMDS2", colour = "Phylum", y = "Bin size (Mbp)", shape = "Type", tag = "B")


NMDS2_binsize

library('rstatix')
nmds_bin_size$cluster <- str_replace(nmds_bin_size$cluster, "clust1", "Cluster 1")
nmds_bin_size$cluster <- str_replace(nmds_bin_size$cluster, "clust2", "Cluster 2")
nmds_bin_size$cluster <- str_replace(nmds_bin_size$cluster, "clust3", "Cluster 3")

cluster_comparisons <- list( c("Cluster 1", "Cluster 2"),  c("Cluster 2", "Cluster 3"), c("Cluster 1", "Cluster 3") )

 

Bin_Size_clust<- ggplot(nmds_bin_size, aes(x = cluster, y = Bin_Size.Mbp.)) +
  geom_boxplot() +
  labs(y = "Bin size (Mbp)", x = "Cluster", tag = "C") +
  theme(axis.title.x = element_text(face="bold",size=24),
        axis.text.x = element_text(colour = "black", size=22, angle = 90),
        axis.text.y = element_text(colour = "black", size=22),
        axis.title.y = element_text(face="bold",size=24),
        plot.title = element_text(size = 24),
        legend.title =element_text(face="bold",size = 14),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1, "cm"),
        strip.text.x = element_text(size=22, face="bold"),
        strip.text.y = element_text(size=22, face="bold"),
        panel.border = element_rect(fill = NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  stat_pvalue_manual(nmds_bin_size %>% 
                       wilcox_test(Bin_Size.Mbp. ~ cluster, exact = FALSE) %>% 
                       add_xy_position())

  #stat_compare_means(comparisons = cluster_comparisons, label = "p.signif") #+
 
Bin_Size_clust

nmds_bin_size$ORF_binsize <- nmds_bin_size$ORF_count/nmds_bin_size$Bin_Size.Mbp.

ORF_binsize_boxplot <- ggplot(nmds_bin_size, aes(x = cluster, y = ORF_binsize)) +
  geom_boxplot() +
  stat_pvalue_manual(nmds_bin_size %>% 
                       wilcox_test(ORF_binsize ~ cluster, exact = NULL) %>% 
                       add_xy_position())+
  theme(axis.title.x = element_text(face="bold",size=24),
        axis.text.x = element_text(colour = "black", size=22, angle = 90, hjust = 1), 
        axis.text.y = element_text(colour = "black", size=22),
        axis.title.y = element_text(face="bold",size=24),
        plot.title = element_text(size = 24),
        legend.title =element_text(face="bold",size = 14),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1, "cm"),
        strip.text.x = element_text(size=22, face="bold"),
        strip.text.y = element_text(size=22, face="bold"),
        panel.border = element_rect(fill = NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  labs(y = "Gene density", x = "Cluster", tag = "D")


ORF_binsize_boxplot

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/Figure_3.pdf", width = 15, height = 13)

grid.arrange(functional_NMDS_Phylum, NMDS2_binsize, Bin_Size_clust, ORF_binsize_boxplot,
             widths = c(3, 1, 1),
             layout_matrix = rbind(c(1, 2, 2),
                                   c(1, 3, 4)))

dev.off()



###Figure 4
#Cluster 1&2 vs Cluster 3 = Clust 1 vs Clust 2
metabolism_summary_for_dotplot <- gather(metabolism_summary, bin, count, "10_bin.1":"9_bin.9")
metabolism_summary_for_dotplot_phylum <- subset(metabolism_summary_for_dotplot, select = -c(gene_description, gene_id, function_description, module, subheader))

conflicts_prefer(dplyr::filter)
metabolism_summary_for_dotplot_phylum <- filter(metabolism_summary_for_dotplot_phylum, count > 0)

metabolism_summary_for_dotplot_phylum <- merge(metabolism_summary_for_dotplot_phylum, mots_taxa, by.x ="bin", by.y = "bin_ID")

metabolism_summary_for_dotplot_Crenarchaeota <- metabolism_summary_for_dotplot_phylum %>% filter(grepl('Crenarchaeota', Phylum))

metabolism_summary_for_dotplot_Crenarchaeota <- subset(metabolism_summary_for_dotplot_Crenarchaeota, select = -c(Watermass))

library(data.table)
setDT(metabolism_summary_for_dotplot_Crenarchaeota)
metabolism_summary_for_dotplot_Crenarchaeota1  = metabolism_summary_for_dotplot_Crenarchaeota [ , .(counted = sum(count)), by = .(header, bin, Phylum)]

metabolism_summary_for_dotplot_Thermoplasmatota <- metabolism_summary_for_dotplot_phylum %>% filter(grepl('Thermoplasmatota', Phylum))

metabolism_summary_for_dotplot_Thermoplasmatota <- subset(metabolism_summary_for_dotplot_Thermoplasmatota, select = -c(Watermass))

library(data.table)
setDT(metabolism_summary_for_dotplot_Thermoplasmatota)
metabolism_summary_for_dotplot_Thermoplasmatota1  = metabolism_summary_for_dotplot_Thermoplasmatota [ , .(counted = sum(count)), by = .(header, bin, Phylum)]


to_merge_gene_counts_Cren <- metabolism_summary_for_dotplot_Crenarchaeota1 %>%
  add_column(add_column = "Cluster 1")
to_merge_gene_counts_Therm <- metabolism_summary_for_dotplot_Thermoplasmatota1 %>%
  add_column(add_column = "Cluster 1")

Cluster1_Cren_Therm_gene <- rbind(to_merge_gene_counts_Cren, to_merge_gene_counts_Therm)

metabolism_summary_for_dotplot <- gather(metabolism_summary, bin, count, "10_bin.1":"9_bin.9")

#mots_taxa <- tibble::rownames_to_column(mots_taxa, "bin_ID")

metabolism_summary_for_dotplot <- merge(metabolism_summary_for_dotplot, mots_taxa, by.x ="bin", by.y = "bin_ID")

metabolism_summary_for_dotplot <- subset(metabolism_summary_for_dotplot, select = -c(gene_description, gene_id, function_description, subheader, Watermass))

library(data.table)
setDT(metabolism_summary_for_dotplot)
metabolism_summary_for_dotplot1  = metabolism_summary_for_dotplot [ , .(counted = sum(count)), by = .(header, bin, Phylum)]

Cluster3_gene <- metabolism_summary_for_dotplot1 %>%
  add_column(add_column = "Cluster 2")

Cluster1_2_Cren_Therm_gene_count <- rbind(Cluster1_Cren_Therm_gene, Cluster3_gene, use.names=FALSE)


library(data.table)
setDT(Cluster1_2_Cren_Therm_gene_count)
Cluster1_2_Cren_Therm_gene_count1  = Cluster1_2_Cren_Therm_gene_count [ , .(counted = sum(counted)), by = .(header, bin, Phylum, counted, add_column)]

names(Cluster1_2_Cren_Therm_gene_count1)[names(Cluster1_2_Cren_Therm_gene_count1) == 'add_column'] <- 'Cluster'

kruskal.test(counted ~ Cluster, data = Cluster1_2_Cren_Therm_gene_count1)

pairwise.wilcox.test(Cluster1_2_Cren_Therm_gene_count1$counted, Cluster1_2_Cren_Therm_gene_count1$Cluster)

Cluster1_2_Cren_Therm_gene_count$header <- factor(Cluster1_2_Cren_Therm_gene_count$header, levels=c("Information systems", "Peptidase", "Amino Acid", "central carbon", "CAZY", "C1", "Electron transport Chain", "MISC", "Flagella Structure", "Oxygen", "Antibiotic Resistance", "pyruvate metabolism", "ADO-CBL synthesis", "Photosynthesis", "Sulfur", "hydrocarbon degradation", "sugar utilization (woodcroft)", "C1-methane", "aerobic corrin ring synthesis", "Flagellar cytoplasmic chaperone", "Metal Reduction", "Nitrogen", "SCFA and alcohol conversions", "CRISPR", "Vitamin B12 transport system", "Hydrogenases", "NA"))

Cluster1_2_Cren_Therm_gene_count$add_column <- str_replace(Cluster1_2_Cren_Therm_gene_count$add_column, "Cluster 1", "Cluster 1+2")

Cluster1_2_Cren_Therm_gene_count$add_column <- str_replace(Cluster1_2_Cren_Therm_gene_count$add_column, "Cluster 2", "Cluster 3")

ggboxplot(Cluster1_2_Cren_Therm_gene_count, x = "header", y = "counted",
          facet.by = "add_column", width = .5) +stat_compare_means(label = "p.signif", label.y = c(45, 55, 55)) +
  labs(x = "Gene category", y = "Total gene count (log)") +
  scale_y_continuous(trans='log10', breaks = c(1,10,100, 1000), limits = c(1,1000)) +
  My_Theme 


##rename headers
Cluster1_2_Cren_Therm_gene_count$header <- str_replace(Cluster1_2_Cren_Therm_gene_count$header,"Amino Acid", "Amino acid")
Cluster1_2_Cren_Therm_gene_count$header <- str_replace(Cluster1_2_Cren_Therm_gene_count$header,"central carbon", "Central carbon")
Cluster1_2_Cren_Therm_gene_count$header <- str_replace(Cluster1_2_Cren_Therm_gene_count$header,"Electron transport Chain", "Electron transport chain")
Cluster1_2_Cren_Therm_gene_count$header <- str_replace(Cluster1_2_Cren_Therm_gene_count$header,"Antibiotic Resistance", "Antibiotic resistance")
Cluster1_2_Cren_Therm_gene_count$header <- str_replace(Cluster1_2_Cren_Therm_gene_count$header,"pyruvate metabolism", "Pyruvate metabolism")
Cluster1_2_Cren_Therm_gene_count$header <- str_replace(Cluster1_2_Cren_Therm_gene_count$header,"Flagella Structure", "Flagella structure")
Cluster1_2_Cren_Therm_gene_count$header <- str_replace(Cluster1_2_Cren_Therm_gene_count$header,"hydrocarbon degradation", "Hydrocarbon degradation")
Cluster1_2_Cren_Therm_gene_count$header <- str_replace(Cluster1_2_Cren_Therm_gene_count$header,"sugar utilization (woodcroft)", "Sugar utilization (Woodcroft)")
Cluster1_2_Cren_Therm_gene_count$header <- str_replace(Cluster1_2_Cren_Therm_gene_count$header,"aerobic corrin ring synthesis", "Aerobic corrin ring synthesis")
Cluster1_2_Cren_Therm_gene_count$header <- str_replace(Cluster1_2_Cren_Therm_gene_count$header,"Metal Reduction", "Metal reduction")


Cluster1_2_Cren_Therm_gene_count_2_panel <- ggplot(Cluster1_2_Cren_Therm_gene_count, aes(x = reorder(header, -counted), y = counted, fill = header)) +
  geom_boxplot() +
   scale_fill_manual(values = c("white", "white", "#E31A1C", "white", "#E31A1C", "white", "white", "#E31A1C", "white", "#E31A1C", "#E31A1C", "white", "white", "white", "#E31A1C", "#E31A1C","#E31A1C", "#E31A1C", "#E31A1C","#E31A1C", "white", "white", "white", "white", "white", "white")) +
  facet_wrap(~add_column, ncol = 1) +
  theme(axis.title.x = element_text(face="bold",size=24),
        axis.text.x = element_text(angle=90, colour = "black", size=22, vjust=0.5, hjust = 1), 
        axis.text.y = element_text(colour = "black", size=22),
        axis.title.y = element_text(face="bold",size=24),
        plot.title = element_text(size = 24),
        legend.title =element_text(face="bold",size = 14),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1, "cm"),
        strip.text.x = element_text(size=22, face="bold"),
        strip.text.y = element_text(size=22, face="bold"),
        panel.border = element_rect(fill = NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none") +
  labs(x = "Functional category", y = "Gene count/cluster (log10 scaling)", tags = 
         "A") +
  scale_y_continuous(trans='log10')
#scale_fill_manual(values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#1ff8ff", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#4b6a53", "#b249d5", "#7edc45", "#5c47b8", "#cfd251"))
Cluster1_2_Cren_Therm_gene_count_2_panel

add_cluster <- subset(with_season_meta, select = c("bin_ID", "cluster"))

msummary_for_dotplot_phylum_clust <- merge(metabolism_summary_for_dotplot_phylum, add_cluster, by.x="bin", by.y="bin_ID")

msummary_for_dotplot_phylum_clust <- filter(msummary_for_dotplot_phylum_clust, count > 0)

msummary_for_dotplot_phylum_clust <- merge(msummary_for_dotplot_phylum_clust, mots_taxa, by.x ="bin", by.y = "bin_ID")

#add_cluster <- subset(with_season_meta, select = c("bin_ID", "cluster"))
msummary_for_dotplot_phylum_clust_boxplot <- merge(msummary_for_dotplot_phylum_clust, add_cluster, by.x="bin", by.y="bin_ID")

msummary_for_dotplot_phylum_clust <- subset(msummary_for_dotplot_phylum_clust_boxplot, select = -c(bin, Watermass.x,Watermass.y, Phylum.y))

library(data.table)
setDT(msummary_for_dotplot_phylum_clust)
msummary_for_dotplot_phylum_clust1  = msummary_for_dotplot_phylum_clust [ , .(counted = sum(count)), by = .(header, cluster.x, Phylum.x)]

phylum_stats <- compare_means(counted ~ cluster.x,  data = msummary_for_dotplot_phylum_clust1,
                              group.by = "header")


###Metabolism by Cren, Thermo and the rest
metabolism_summary_for_dotplot <- gather(metabolism_summary, bin, count, "10_bin.1":"9_bin.9")
metabolism_summary_for_dotplot_phylum <- subset(metabolism_summary_for_dotplot, select = -c(gene_description, gene_id, function_description, module, subheader))



metabolism_summary_for_dotplot_phylum <- filter(metabolism_summary_for_dotplot_phylum, count > 0)

#metabolism_summary_for_dotplot_phylum <- merge(metabolism_summary_for_dotplot_phylum, mots_taxa, by.x ="bin", by.y = "bin_ID")

metabolism_summary_a_p_n <- metabolism_summary_for_dotplot_phylum %>% filter(grepl('Amino Acid|Peptidase|Nitrogen', header))
add_cluster <- subset(with_season_meta, select = c("bin_ID", "cluster"))
metabolism_summary_apn_boxplot <- merge(metabolism_summary_a_p_n, add_cluster, by.x="bin", by.y="bin_ID")

#metabolism_summary_apn_boxplot <- subset(metabolism_summary_apn_boxplot, select = -c(bin, Watermass.x,Phylum.x))

metabolism_summary_apn_boxplot$cluster <- str_replace(metabolism_summary_apn_boxplot$cluster, "clust1", "Cluster 1")
metabolism_summary_apn_boxplot$cluster <- str_replace(metabolism_summary_apn_boxplot$cluster, "clust2", "Cluster 2")
metabolism_summary_apn_boxplot$cluster <- str_replace(metabolism_summary_apn_boxplot$cluster, "clust3", "Cluster 3")

metabolism_summary_apn_boxplot$header<-factor(metabolism_summary_apn_boxplot$header,levels =c("Amino Acid", "Peptidase", "Nitrogen"))

metabolism_summary_apn_boxplot$cluster<-factor(metabolism_summary_apn_boxplot$cluster,levels =c("Cluster 1", "Cluster 2", "Cluster 3"))

library(data.table)
setDT(metabolism_summary_apn_boxplot)
metabolism_summary_apn_boxplot1  = metabolism_summary_apn_boxplot [ , .(counted = sum(count)), by = .(header, cluster, bin)]

###TRY NORMALIZE BY BIN COUNT
metabolism_summary_apn_boxplot1$bin_count <- metabolism_summary_apn_boxplot1$cluster

rep_str = c('Cluster 1'="26",'Cluster 2'="122",'Cluster 3'="857")
metabolism_summary_apn_boxplot1$bin_count <- str_replace_all(metabolism_summary_apn_boxplot1$bin_count, rep_str)

metabolism_summary_apn_boxplot1$bin_count[metabolism_summary_apn_boxplot1$bin_count == "Cluster 1"] <- 26
metabolism_summary_apn_boxplot1$bin_count[metabolism_summary_apn_boxplot1$bin_count == "Cluster 2"] <- 122
metabolism_summary_apn_boxplot1$bin_count[metabolism_summary_apn_boxplot1$bin_count == "Cluster 3"] <- 857

metabolism_summary_apn_boxplot1$bin_count <- as.numeric(as.character(metabolism_summary_apn_boxplot1$bin_count))

metabolism_summary_apn_boxplot1$norm <- metabolism_summary_apn_boxplot1$counted/metabolism_summary_apn_boxplot1$bin_count

unique_functions_count <- as.numeric(unique_functions_del$count)
unique_functions_clustnum <- as.numeric(unique_functions_del$clustnum)
unique_functions_norm <- unique_functions_count/unique_functions_clustnum


metabolism_summary_a.a_n_p <- ggplot(metabolism_summary_apn_boxplot1, aes(x = header, y = counted)) +
  geom_boxplot() +
  facet_wrap(~cluster) +
  labs(x = "Gene", y = "Total gene count per cluster (log10)", tag = "B") +
  theme(axis.title.x = element_text(face="bold",size=24),
axis.text.x = element_text(angle=90, colour = "black", size=22, vjust=0.5, hjust = 1), 
axis.text.y = element_text(colour = "black", size=22),
axis.title.y = element_text(face="bold",size=24),
plot.title = element_text(size = 24),
legend.title =element_text(face="bold",size = 14),
legend.text = element_text(size = 14),
legend.key.size = unit(1, "cm"),
strip.text.x = element_text(size=22, face="bold"),
strip.text.y = element_text(size=22, face="bold"),
panel.border = element_rect(fill = NA, colour = "black"),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.background = element_blank(),
legend.position = "none") +
  scale_y_continuous(trans = 'log10 scaling')

metabolism_summary_a.a_n_p

metabolism_summary_a.a_n_p_norm <- ggplot(metabolism_summary_apn_boxplot1, aes(x = header, y = norm)) +
  geom_boxplot() +
  facet_wrap(~cluster) +
  labs(x = "Gene", y = "Total gene count per cluster (log10)", tag = "B") +
  theme(axis.title.x = element_text(face="bold",size=24),
        axis.text.x = element_text(angle=90, colour = "black", size=22, vjust=0.5, hjust = 1), 
        axis.text.y = element_text(colour = "black", size=22),
        axis.title.y = element_text(face="bold",size=24),
        plot.title = element_text(size = 24),
        legend.title =element_text(face="bold",size = 14),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1, "cm"),
        strip.text.x = element_text(size=22, face="bold"),
        strip.text.y = element_text(size=22, face="bold"),
        panel.border = element_rect(fill = NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none") +
  scale_y_continuous(trans = 'log10 scaling')

metabolism_summary_a.a_n_p_norm

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/Figure_4.pdf",15, height = 13)

grid.arrange(Cluster1_2_Cren_Therm_gene_count_2_panel, metabolism_summary_a.a_n_p,
             widths = c(2, 2))

#grid.arrange(Cluster1_Cren_Therm_gene_plot, Cluster3_gene_plot, metabolism_summary_a.a_n_p,
#layout_matrix = rbind(c(1, 2),
# c(3)))

dev.off()


###Figure 5
###TRYING TO COMPARE STATS BETWEEN ALL HEADER
add_cluster <- subset(with_season_meta, select = c("bin_ID", "cluster"))
msummary_for_dotplot_clust <- merge(metabolism_summary_for_dotplot_phylum, add_cluster, by.x="bin", by.y="bin_ID")

#msummary_for_dotplot_clust <- subset(msummary_for_dotplot_clust, select = c(bin, count, cluster))

library(data.table)
setDT(msummary_for_dotplot_clust)
msummary_for_dotplot_clust1  = msummary_for_dotplot_clust [ , .(counted = sum(count)), by = .(header,cluster,bin)]

msummary_for_dotplot_clust1$cluster <- str_replace(msummary_for_dotplot_clust1$cluster, "clust1", "Cluster 1")
msummary_for_dotplot_clust1$cluster <- str_replace(msummary_for_dotplot_clust1$cluster, "clust2", "Cluster 2")
msummary_for_dotplot_clust1$cluster <- str_replace(msummary_for_dotplot_clust1$cluster, "clust3", "Cluster 3")

msummary_for_dotplot_clust$cluster <- str_replace(msummary_for_dotplot_clust$cluster, "clust1", "Cluster 1")
msummary_for_dotplot_clust$cluster <- str_replace(msummary_for_dotplot_clust$cluster, "clust2", "Cluster 2")
msummary_for_dotplot_clust$cluster <- str_replace(msummary_for_dotplot_clust$cluster, "clust3", "Cluster 3")

ggplot(msummary_for_dotplot_clust, aes(x = header, y = count)) +
  geom_boxplot() +
  facet_wrap(~cluster) +
  My_Theme

ggplot(msummary_for_dotplot_clust1, aes(x = header, y = counted)) +
  geom_boxplot(aes(fill = cluster)) +
  facet_wrap(~cluster) +
  My_Theme +
  labs(x = "Energy source", y = "Total gene count (log)") +
  stat_compare_means(paired = TRUE) +
  scale_y_continuous(trans='log10')

cluster_stats <- compare_means(count ~ cluster,  data = msummary_for_dotplot_clust,
                               group.by = "header")

my_comparisons <- list( c("Cluster 1", "Cluster 2"), c("Cluster 2", "Cluster 3"), c("Custer 1", "Cluster 3") )

msummary_for_dotplot_clust$cluster <- factor(msummary_for_dotplot_clust$cluster, levels=c("Cluster 1", "Custer 2", "Cluster 3"))



msummary_for_dotplot_clust1$header <- str_replace(msummary_for_dotplot_clust1$header, "ADO-CBL synthesis", "ADO-CBL syn")
msummary_for_dotplot_clust1$header <- str_replace(msummary_for_dotplot_clust1$header, "anAerobic corrin ring synthesis", "Anaerobic corrin ring synthesis")
msummary_for_dotplot_clust1$header <- str_replace(msummary_for_dotplot_clust1$header,"Amino Acid", "Amino acid")
msummary_for_dotplot_clust1$header <- str_replace(msummary_for_dotplot_clust1$header,"central carbon", "Central carbon")
msummary_for_dotplot_clust1$header <- str_replace(msummary_for_dotplot_clust1$header,"Electron transport Chain", "Electron transport chain")
msummary_for_dotplot_clust1$header <- str_replace(msummary_for_dotplot_clust1$header,"Antibiotic Resistance", "Antibiotic resistance")
msummary_for_dotplot_clust1$header <- str_replace(msummary_for_dotplot_clust1$header,"pyruvate metabolism", "Pyruvate metabolism")
msummary_for_dotplot_clust1$header <- str_replace(msummary_for_dotplot_clust1$header,"Flagella Structure", "Flagella structure")
msummary_for_dotplot_clust1$header <- str_replace(msummary_for_dotplot_clust1$header,"hydrocarbon degradation", "Hydrocarbon degradation")
msummary_for_dotplot_clust1$header <- str_replace(msummary_for_dotplot_clust1$header,"sugar utilization (woodcroft)", "Sugar utilization (Woodcroft)")
msummary_for_dotplot_clust1$header <- str_replace(msummary_for_dotplot_clust1$header,"aerobic corrin ring synthesis", "Aerobic corrin ring synthesis")
msummary_for_dotplot_clust1$header <- str_replace(msummary_for_dotplot_clust1$header,"Metal Reduction", "Metal reduction")


msummary_phylum_clust_pallete = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#1ff8ff", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#4b6a53", "#b249d5", "#7edc45", "#5c47b8", "#cfd251", "#ff69b4", "#69c86c")

msummary_for_dotplot_clust1 = msummary_for_dotplot_clust1[!grepl("Aerobic corrin ring synthesis", msummary_for_dotplot_clust1$header),]
msummary_for_dotplot_clust1 = msummary_for_dotplot_clust1[!grepl("Anaerobic corrin ring synthesis", msummary_for_dotplot_clust1$header),]
msummary_for_dotplot_clust1 = msummary_for_dotplot_clust1[!grepl("CRISPR", msummary_for_dotplot_clust1$header),]
msummary_for_dotplot_clust1 = msummary_for_dotplot_clust1[!grepl("Flagella structure", msummary_for_dotplot_clust1$header),]
msummary_for_dotplot_clust1 = msummary_for_dotplot_clust1[!grepl("Flagellar cytoplasmic chaperone", msummary_for_dotplot_clust1$header),]
msummary_for_dotplot_clust1 = msummary_for_dotplot_clust1[!grepl("Hydrocarbon degradation", msummary_for_dotplot_clust1$header),]
msummary_for_dotplot_clust1 = msummary_for_dotplot_clust1[!grepl("Hydrogenases", msummary_for_dotplot_clust1$header),]
msummary_for_dotplot_clust1 = msummary_for_dotplot_clust1[!grepl("Metal reduction", msummary_for_dotplot_clust1$header),]
msummary_for_dotplot_clust1 = msummary_for_dotplot_clust1[!grepl("Nitrogen", msummary_for_dotplot_clust1$header),]
msummary_for_dotplot_clust1 = msummary_for_dotplot_clust1[!grepl("Pyruvate metabolism", msummary_for_dotplot_clust1$header),]
msummary_for_dotplot_clust1 = msummary_for_dotplot_clust1[!grepl("SCFA  and alcohol conversions", msummary_for_dotplot_clust1$header),]
msummary_for_dotplot_clust1 = msummary_for_dotplot_clust1[!grepl('sugar', msummary_for_dotplot_clust1$header),]
msummary_for_dotplot_clust1 = msummary_for_dotplot_clust1[!grepl("Vitamin B12 transport system", msummary_for_dotplot_clust1$header),]
msummary_for_dotplot_clust1 = msummary_for_dotplot_clust1[!grepl("C1", msummary_for_dotplot_clust1$header),]

msummary_for_dotplot_clust1$header = factor(msummary_for_dotplot_clust1$header, levels=c('Information systems','Amino acid','Peptidase','Central carbon','CAZY','C1','Electron transport chain','MISC','Sulfur','Oxygen','Antibiotic resistance','Pyruvate metabolism','Photosynthesis','ADO_CBL syn','Sulfur','Flagella structure','Hydrocarbon degradation','Sugar utilization (woodcroft)', 'C1-methane', 'Aerobic corrin ring synthesis', 'Nitrogen', 'Flagellar cytoplasmic chaperone','Metal reduction', 'CRISPR', 'Vitamin B12 transport system','Hydrogenases', 'NA'))

msummary_for_dotplot_clust1$header = factor(msummary_for_dotplot_clust1$header, levels=c('Information systems','Amino acid','Peptidase','Central carbon','CAZY','C1','Electron transport chain','MISC','Sulfur','Oxygen','Antibiotic resistance','C1-methane','Photosynthesis','ADO_CBL synthesis','NA'))
msummary_for_dotplot_clust1

cluster_comparisons <- list( c("Cluster 1", "Cluster 2"),  c("Cluster 2", "Cluster 3"), c("Cluster 1", "Cluster 3") )

msummary_for_dotplot_clust1$bin_count <- msummary_for_dotplot_clust1$cluster
rep_str = c('Cluster 1'="26",'Cluster 2'="122",'Cluster 3'="857")
msummary_for_dotplot_clust1$bin_count <- str_replace_all(msummary_for_dotplot_clust1$bin_count, rep_str)

msummary_for_dotplot_clust1$bin_count[msummary_for_dotplot_clust1$bin_count == "Cluster 1"] <- 26
msummary_for_dotplot_clust1$bin_count[msummary_for_dotplot_clust1$bin_count == "Cluster 2"] <- 122
msummary_for_dotplot_clust1$bin_count[msummary_for_dotplot_clust1$bin_count == "Cluster 3"] <- 857

msummary_for_dotplot_clust1$bin_count <- as.numeric(as.character(msummary_for_dotplot_clust1$bin_count))

msummary_for_dotplot_clust1$norm <- msummary_for_dotplot_clust1$counted/msummary_for_dotplot_clust1$bin_count


msummary_for_dotplot_clust1$norm <- as.numeric(msummary_for_dotplot_clust1$norm)


msummary_for_dotplot_clust_plot <- ggplot(msummary_for_dotplot_clust1, aes(x = cluster, y = counted, fill = cluster)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#CAB2D6", "lightgrey", "#FFFFB0")) +
  facet_wrap(~header, scales = "free") +
  stat_compare_means(comparisons = cluster_comparisons, hide.ns = TRUE, label = "p.signif") +
  theme(axis.title.x = element_text(face="bold",size=24),
        axis.text.x = element_text(angle=90, colour = "black", size=22, vjust=0.5, hjust = 1), 
        axis.text.y = element_text(colour = "black", size=22),
        axis.title.y = element_text(face="bold",size=24),
        plot.title = element_text(size = 24),
        legend.title =element_text(face="bold",size = 14),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1, "cm"),
        strip.text.x = element_text(size=18, face="bold"),
        strip.text.y = element_text(size=18, face="bold"),
        panel.border = element_rect(fill = NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none") +
  labs(x = "Cluster", y = "Normlized gene count/cluster (log10 scaling)") +
  scale_y_continuous(trans='log10')

msummary_for_dotplot_clust_plot

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/Figure_5_free_scale.pdf",15, height = 14)

msummary_for_dotplot_clust_plot

dev.off()

###Figure 6
Munida_mapping <- ps_melt(mots.phylo)




Munida_mapping$Phylum <- factor(Munida_mapping$Phylum, levels=c("SAR324", "Verrucomicrobiota", "Proteobacteria", "Cyanobacteria", "Planctomycetota", "Gemmatimonadota", "Chloroflexota", "Nitrospinota", "Myxococcota", "Marinisomatota", "Dadabacteria","Fibrobacterota","Bacteroidota", "Binatota","Crenarchaeota"))

Munida_mapping$Watermass <- factor(Munida_mapping$Watermass, levels=c('Neritic','STW','Front','SAW','Deep'))

Munida_mapping$Depth = Munida_mapping$Watermass

Munida_mapping$Depth <- str_replace(Munida_mapping$Depth, "Neritic", "Surface")
Munida_mapping$Depth <- str_replace(Munida_mapping$Depth, "STW", "Surface")
Munida_mapping$Depth <- str_replace(Munida_mapping$Depth, "Front", "Surface")
Munida_mapping$Depth <- str_replace(Munida_mapping$Depth, "SAW", "Surface")

Munida_mapping$Depth <- factor(Munida_mapping$Depth, levels=c('Surface', 'Deep'))


Munida_mapping$Phylum <- factor(Munida_mapping$Phylum, levels=c('Marinisomatota','SAR324','Crenarchaeota','Cyanobacteria','Dadabacteria','Nitrospinota','Bacteroidota','Proteobacteria','Gemmatimonadota','Thermoplasmatota','Actinobacteriota','Acidobacteriota','Chloroflexota','Verrucomicrobiota','Planctomycetota','Binatota','Myxococcota','Fibrobacterota'))

Munida_Tara_Comparison <- ggplot(Munida_mapping, aes(x= Watermass, y=Abundance, group = Phylum, fill=Phylum)) + 
  facet_wrap(~Depth , scales = "free")+ 
  geom_bar(stat="summary", color = "black")+ 
  labs(x = "Location", y = "Reads per million (RPM)", tag = "A") +
  scale_fill_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#C3B4D1", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#4B9B7A", "Thermoplasmatota" = "lightgrey")) +
  theme_light()+
  theme(axis.title.x = element_text(face="bold",size=12),
        axis.text.x = element_text(angle=45, colour = "black", size=10,vjust = 0.5, hjust=0.5), 
        axis.text.y = element_text(colour = "black", size=10),
        axis.title.y = element_text(face="bold",size=12),
        plot.title = element_text(size = 12),
        legend.title =element_text(face="bold",size = 12),
        legend.text = element_text(size = 10),
        legend.position="none",
        legend.key.size = unit(1, "cm"),
        panel.background = element_blank(),
        strip.text.x = element_text(size=12, color="black"),
        panel.border = element_rect(fill = NA, colour = "black"))

###FOR TARA
#OTU table
tara.mg.otu.csv <- read.csv(file="/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/read_recruitment_to_bins/import/tara_mg_otu.csv", header = TRUE, na.strings = "unclassified", fileEncoding="UTF-8-BOM", check.names = FALSE, row.names = 1)
tara.mg.otu.matrix <- as.matrix(tara.mg.otu.csv)

tara.mg.otu <- otu_table(tara.mg.otu.matrix, taxa_are_rows = TRUE)


#Taxa table
tara.mg.taxa.csv <- read.csv(file="/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/read_recruitment_to_bins/import/taxa_table.csv", header = TRUE, na.strings = "unclassified", fileEncoding="UTF-8-BOM", check.names = FALSE, row.names = 1)
tara.mg.taxa.matrix <- as.matrix(tara.mg.taxa.csv)

tara.mg.taxa <- tax_table(tara.mg.taxa.matrix)


#Meta-data
tara.mg.sample.meta <- read.csv(file="/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/read_recruitment_to_bins/import/tara_mg_meta.csv", header = TRUE, fileEncoding="UTF-8-BOM", check.names = FALSE, row.names = 1)
tara.mg.meta.df <- as.data.frame(tara.mg.sample.meta)

tara.mg.meta <- sample_data(tara.mg.meta.df)


#Create phyloseq object
tara.mg.phylo <- phyloseq(tara.mg.otu, tara.mg.taxa, tara.mg.meta)


#Melt phyloseq object for figure-making df
tara.mg.mapping <- psmelt(tara.mg.phylo)


##Fix column names
tara.mg.mapping<-tara.mg.mapping %>% rename(c('X16S' = '16S',
                                              'Community....' = 'Community_%',
                                              'Bin_Size..Mbp.' = 'Bin Size (Mbp)',
                                              'Completeness....' = 'Completeness (%)',
                                              'Contamination....' = 'Contamination (%)',
                                              'Strain_Heterogeneity....' = 'Strain_Heterogeneity (%)',
                                              'OTU' = 'Bin_ID'))


colnames(tara.mg.mapping)



tara_mg <- filter(tara.mg.mapping)

#import thermo data...remove in cleanup since used in other figures
#thermo_meta <- read.csv("/Users/morse47p/Library/CloudStorage/OneDrive-UniversityofOtago/Dropbox/Alida_Thermo_bins/figure_code/raw_data/Thermoplasmatota_base_data.csv",fill = TRUE, header = TRUE, sep = ",")

#Bin_ID<-c(tara.mg.mapping$Bin_ID)

#tara_mg<-filter(tara.mg.mapping, Bin_ID)

#save table with data
#write.csv(tara_mg, "/Users/morse47p/Library/CloudStorage/OneDrive-UniversityofOtago/Dropbox/Alida_Thermo_bins/figure_code/raw_data/tara_mg.csv", row.names = T)

library(scales)

#reorder locations
tara_mg$Location <- factor(tara_mg$Location, levels=c("PON", "PSE", "PSW", "ANE", "ANW", "ASE","ASW","ION","IOS","SOC", "MED","RED"))
#reorder depths
tara_mg$Depth <- factor(tara_mg$Depth, levels=c("SUR", "DCM", "MES"))
tara_mg$mOTU<-as.character(tara_mg$mOTU)
#remove leading space
tara_mg$mOTU <- trimws(tara_mg$mOTU, which = c("left"))

library(microViz)
brewerPlus <- distinct_palette()

tara_mg$Phylum <- reorder(tara_mg$Phylum, tara_mg$Abundance)
tara_mg$Phylum <- factor(tara_mg$Phylum, levels=rev(levels(tara_mg$Phylum)))

Tara_Comparison <- ggplot(tara_mg, aes(x= Location, y=Abundance, group = Phylum, fill=Phylum)) + 
  facet_wrap(~Depth ,
             labeller = labeller(Habitat = label_wrap_gen(width = 20, multi_line=TRUE)))+ 
  geom_bar(stat="summary", color = "black")+ 
  ylab("reads per million (RPM)")+ 
  xlab("Location")+
  theme_light()+
  theme(axis.title.x = element_text(face="bold",size=12),
        axis.text.x = element_text(angle=45, colour = "black", size=10,vjust = 0.5, hjust=0.5), 
        axis.text.y = element_text(colour = "black", size=10),
        axis.title.y = element_text(face="bold",size=12),
        plot.title = element_text(size = 12),
        legend.title =element_text(face="bold",size = 12),
        legend.text = element_text(size = 10),
        legend.position="right",
        legend.key.size = unit(1, "cm"),
        panel.background = element_blank(),
        strip.text.x = element_text(size=12, color="black"),
        panel.border = element_rect(fill = NA, colour = "black")) +
  scale_fill_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#C3B4D1", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#4B9B7A", "Thermoplasmatota" = "lightgrey"))

legend <- get_legend(Tara_Comparison)

Tara_Comparison <- ggplot(tara_mg, aes(x= Location, y=Abundance, group = Phylum, fill=Phylum)) + 
  facet_wrap(~Depth ,
             labeller = labeller(Habitat = label_wrap_gen(width = 20, multi_line=TRUE)))+ 
  geom_bar(stat="summary", color = "black")+ 
  labs(x = "Location", y = "Reads per million (RPM)", tag = "B") +
  theme_light()+
  theme(axis.title.x = element_text(face="bold",size=12),
        axis.text.x = element_text(angle=45, colour = "black", size=10,vjust = 0.5, hjust=0.5), 
        axis.text.y = element_text(colour = "black", size=10),
        axis.title.y = element_text(face="bold",size=12),
        plot.title = element_text(size = 12),
        legend.title =element_text(face="bold",size = 12),
        legend.text = element_text(size = 10),
        legend.position="none",
        legend.key.size = unit(1, "cm"),
        panel.background = element_blank(),
        strip.text.x = element_text(size=12, color="black"),
        panel.border = element_rect(fill = NA, colour = "black")) +
  scale_fill_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#C3B4D1", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#4B9B7A", "Thermoplasmatota" = "lightgrey"))

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/Figure_6.pdf",15, height = 13)

grid.arrange(Munida_Tara_Comparison, Tara_Comparison, legend,
             heights = c(1,2),
             layout_matrix = rbind(c(1, 3),
                                   c(2, 3)))

#grid.arrange(Cluster1_Cren_Therm_gene_plot, Cluster3_gene_plot, metabolism_summary_a.a_n_p,
#layout_matrix = rbind(c(1, 2),
# c(3)))

dev.off()










###HEATMAP FOR Central Carbon Metabolism and ETC
remove(metabolism_summary)
metabolism_summary<-dplyr::bind_rows(MISC_combined,carbon_utilization_combined,Transporters_combined,Energy_combined,organic_nitrogen_combined,Woodcroft_combined,rRNA_combined)

names(metabolism_summary)[1] <- "gene_id"
names(metabolism_summary)[2] <- "gene_description"
names(metabolism_summary)[3] <- "function_description"
names(metabolism_summary)[4] <- "module"
names(metabolism_summary)[5] <- "header"
names(metabolism_summary)[6] <- "subheader"


metabolism_summary$unique_id <- paste(metabolism_summary$function_description, metabolism_summary$gene_id,  metabolism_summary$module, sep="_")

#metabolism_summary<-dplyr::bind_rows(MISC1,carbon_utilization1,Transporters1,Energy1,Organic_Nitrogen1,Woodcroft1,rRNA1,tRNA1)



metabolism_summary$unique_id <- paste(metabolism_summary$function_description, metabolism_summary$gene_id,  metabolism_summary$module, sep="_")
metabolism_summary<-select(metabolism_summary,unique_id, everything())





conflicts_prefer(rstatix::filter)
metabolism_summary_for_central_ETC <- metabolism_summary_for_central_ETC %>%
  filter(grepl('central carbon|Electron transport Chain', header))

metabolism_summary_for_central_ETC <- subset(metabolism_summary_for_central_ETC, select = c(gene_description, bin, count))

library(data.table)
setDT(metabolism_summary_for_central_ETC)
metabolism_summary_for_central_ETC_counted  = metabolism_summary_for_central_ETC [ , .(counted = sum(count)), by = .(gene_description, bin)]

###HEATMAP FUNCTION
met_summary_central_ETC_matrix <- reshape(metabolism_summary_for_central_ETC_counted, idvar = "gene_description", timevar = "bin", direction = "wide")

met_summary_central_ETC_matrix <- met_summary_central_ETC_matrix %>% 
  column_to_rownames(var = "gene_description")

met_summary_central_ETC_matrix <- as.matrix(met_summary_central_ETC_matrix)

central_carbon_ETC_heatmap <- pheatmap(met_summary_central_ETC_matrix, scale = 'column',
        annotation_row = col_annotation,
        annotation_color = header_color,
        show_colnames=TRUE,
        drop_levels = FALSE,
        annotation_legend= TRUE,
        annotation_names_row = TRUE,
        legend_breaks = c(1, 0),
        fontsize = 10, 
        fontsize_row = 6, 
        fontsize_col =3)

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/central_carbon_ETC_heatmap.pdf",30, height = 20)
central_carbon_ETC_heatmap
dev.off()

heatmap(met_summary_central_ETC_matrix, scale = 'column',
        Colv = col_annotation,
        ColSideColors = header_color)

metabolism_summary_for_central_ETC <- gather(metabolism_summary, bin, count, "10_bin.1":"9_bin.9")

conflicts_prefer(rstatix::filter)
metabolism_summary_for_central_ETC <- metabolism_summary_for_central_ETC %>%
  filter(grepl('central carbon|Electron transport Chain', header))

metabolism_summary_for_central_ETC <- subset(metabolism_summary_for_central_ETC, select = c(gene_description, header, count))


library(data.table)
setDT(metabolism_summary_for_central_ETC)
metabolism_summary_for_central_ETC_counted  = metabolism_summary_for_central_ETC [ , .(counted = sum(count)), by = .(gene_description, header)]

metabolism_summary_for_central_ETC_annotation <- subset(metabolism_summary_for_central_ETC_counted, select = c(gene_description, header))
col_annotation <- metabolism_summary_for_central_ETC_annotation %>%
  column_to_rownames(var = "gene_description")

header_color <- list(
  header = c("central carbon"="#296A8F", "Electron transport Chain"="#5BFB72"))


metabolism_summary_for_central_ETC <- gather(metabolism_summary, bin, count, "10_bin.1":"9_bin.9")
conflicts_prefer(rstatix::filter)
metabolism_summary_for_central_ETC <- metabolism_summary_for_central_ETC %>%
  filter(grepl('central carbon|Electron transport Chain', header))

metabolism_summary_for_central_ETC_Phylum <- merge(metabolism_summary_for_central_ETC, with_season_meta, by.y = "bin_ID", by.x = "bin")

metabolism_summary_for_central_ETC_Phylum <- subset(metabolism_summary_for_central_ETC_Phylum, select = c(gene_description, header, Phylum, count))


library(data.table)
setDT(metabolism_summary_for_central_ETC_Phylum)
metabolism_summary_for_central_ETC_Phylum_counted  = metabolism_summary_for_central_ETC_Phylum [ , .(counted = sum(count)), by = .(gene_description, Phylum)]

metabolism_summary_for_central_ETC_Phylum_annotation <- subset(metabolism_summary_for_central_ETC_Phylum_counted, select = c(gene_description))
col_annotation <- metabolism_summary_for_central_ETC_annotation %>%
  column_to_rownames(var = "gene_description")

met_summary_central_ETC_Phylum_matrix <- reshape(metabolism_summary_for_central_ETC_Phylum_counted, idvar = "gene_description", timevar = "Phylum", direction = "wide")

met_summary_central_ETC_Phylum_matrix <- met_summary_central_ETC_Phylum_matrix %>% 
  column_to_rownames(var = "gene_description")

met_summary_central_ETC_Phylum_matrix <- as.matrix(met_summary_central_ETC_Phylum_matrix)

central_carbon_ETC_Phylum_heatmap <- pheatmap(met_summary_central_ETC_Phylum_matrix, scale = 'column',
                                              filename = "~/Desktop/heatmap6.pdf")
central_carbon_ETC_Phylum_heatmap


###FROM NITROSO PAPER
allbins_metabolism_summary<-metabolism_summary

conflicts_prefer(rstatix::filter)
allbins_metabolism_summary <- allbins_metabolism_summary %>%
  filter(grepl('central carbon|Electron transport Chain', header))

#combine into single file
#metabolism_summary1<-dplyr::bind_rows(MISC1,carbon_utilization1,Transporters1,Energy1,Organic_Nitrogen1,Woodcroft1,rRNA1,tRNA1)
#metabolism_summary2<-dplyr::bind_rows(MISC2,carbon_utilization2,Transporters2,Energy2,Organic_Nitrogen2,Woodcroft2,rRNA2,tRNA2)

#combine into all bin results into a single file
#allbins_metabolism_summary<-merge.data.frame(metabolism_summary2,metabolism_summary1, all.x = T, all.y=T)


###X IS BY unique_ID
allbins_metabolism_summary$unique_id <- paste(allbins_metabolism_summary$function_description, allbins_metabolism_summary$gene_id,  allbins_metabolism_summary$module, sep="_")
allbins_metabolism_summary<-select(allbins_metabolism_summary,unique_id, everything())

allbins_metabolism_summary_nr<-allbins_metabolism_summary %>% 
  distinct(unique_id,.keep_all = TRUE)

#allbins_metabolism_summary_nr <- allbins_metabolism_summary_nr %>%
  #filter(grepl('central carbon|Electron transport Chain', header))

#retain only presence/absence data
allbins_metabolism_summary_nr_ap <- allbins_metabolism_summary_nr[ -c(2:7) ]
allbins_metabolism_summary[is.na(allbins_metabolism_summary)] <- 0 



#make functional data into matrix
allbins_metabolism_summary_nr_ap<-column_to_rownames(allbins_metabolism_summary_nr_ap, var = "unique_id")
transposed_metabolism_summary<-t(allbins_metabolism_summary_nr_ap)
transposed_metabolism_summary[is.na(transposed_metabolism_summary)] <- 0

transposed_metabolism_summary<-data.frame(transposed_metabolism_summary)
metabolism_summary_onlypresent<-transposed_metabolism_summary[, colSums(transposed_metabolism_summary) > 0,]





matrix_abundance_metabolism <- as.matrix(metabolism_summary_onlypresent)

Abund_functions_matrix <- mapply(matrix_abundance_metabolism, FUN=as.numeric)
Abund_functions_matrix <- matrix(data=Abund_functions_matrix, ncol=length(colnames(matrix_abundance_metabolism)), nrow=length(row.names(matrix_abundance_metabolism)))
row.names(Abund_functions_matrix) <- row.names(matrix_abundance_metabolism)
colnames(Abund_functions_matrix) <- colnames(matrix_abundance_metabolism)



#run this to check if you have any NAs
is.na(Abund_functions_matrix) %>% table()

dim(Abund_functions_matrix)

#heatmap

pheatmap(Abund_functions_matrix,
         filename ="~/Desktop/heatmap1.pdf")

#convert to presence absence and replot wihtout labels
library(vegan)
PA_functions_matrix<-decostand(Abund_functions_matrix, method="pa")
pheatmap(PA_functions_matrix,show_colnames=FALSE,
         filename ="~/Desktop/heatmap2.pdf")

###X IS BY HEADER

allbins_metabolism_summary$unique_id <- paste(allbins_metabolism_summary$function_description, allbins_metabolism_summary$gene_id,  allbins_metabolism_summary$module, sep="_")
allbins_metabolism_summary<-select(allbins_metabolism_summary,unique_id, everything())

allbins_metabolism_summary_nr<-allbins_metabolism_summary %>% 
  distinct(unique_id,.keep_all = TRUE)

#retain only presence/absence data
allbins_metabolism_summary_nr_ap <- allbins_metabolism_summary_nr[ -c(1:5, 7) ]
allbins_metabolism_summary[is.na(allbins_metabolism_summary)] <- 0 

library(plyr)
allbins_metabolism_summary_nr_ap <- ddply(allbins_metabolism_summary_nr_ap,"header",numcolwise(sum))

#make functional data into matrix
allbins_metabolism_summary_nr_ap$header <- gsub("-", "", allbins_metabolism_summary_nr_ap$header)
allbins_metabolism_summary_nr_ap$header <- gsub("(", "", allbins_metabolism_summary_nr_ap$header)
allbins_metabolism_summary_nr_ap$header <- gsub(")", "", allbins_metabolism_summary_nr_ap$header)

allbins_metabolism_summary_nr_ap$header[is.na(allbins_metabolism_summary_nr_ap$header)] <- "NA"

allbins_metabolism_summary_nr_ap<-column_to_rownames(allbins_metabolism_summary_nr_ap, var = "header")

transposed_metabolism_summary<-t(allbins_metabolism_summary_nr_ap)
transposed_metabolism_summary[is.na(transposed_metabolism_summary)] <- 0

transposed_metabolism_summary<-data.frame(transposed_metabolism_summary)
metabolism_summary_onlypresent<-transposed_metabolism_summary[, colSums(transposed_metabolism_summary) > 0,]





matrix_abundance_metabolism <- as.matrix(metabolism_summary_onlypresent)

Abund_functions_matrix <- mapply(matrix_abundance_metabolism, FUN=as.numeric)
Abund_functions_matrix <- matrix(data=Abund_functions_matrix, ncol=length(colnames(matrix_abundance_metabolism)), nrow=length(row.names(matrix_abundance_metabolism)))
row.names(Abund_functions_matrix) <- row.names(matrix_abundance_metabolism)
colnames(Abund_functions_matrix) <- colnames(matrix_abundance_metabolism)



#run this to check if you have any NAs
is.na(Abund_functions_matrix) %>% table()

dim(Abund_functions_matrix)

#heatmap

pheatmap(Abund_functions_matrix,
         #drop_levels = FALSE,
         #annotation_legend= TRUE,
         #annotation_names_row = TRUE,
         legend_breaks = c(1, 0),
         fontsize = 10, 
         fontsize_row = 10, 
         fontsize_col =10,
         cellheight = 10,
         angle_col = 45,
         filename ="~/Desktop/heatmap3.pdf")

pheatmap(Abund_functions_matrix,
         filename ="~/Desktop/heatmap34.pdf",
         fontfamily = "mono",
         cellheight=10)

#convert to presence absence and replot wihtout labels
library(vegan)
PA_functions_matrix<-decostand(Abund_functions_matrix, method="pa")
pheatmap(PA_functions_matrix,show_colnames=FALSE,
         filename ="~/Desktop/heatmap4.pdf")


#build annotated heatmap

#create a document with gene annotations using the non redundant dataframe
all_functions <- allbins_metabolism_summary_nr[ -c(8:1015) ]
all_functional_categories<-distinct(all_functions,unique_id, .keep_all= TRUE)
all_functional_categories$unique_id<-gsub(" ", ".", all_functional_categories$unique_id)

all_functional_categories$unique_id<-gsub(",", ".", all_functional_categories$unique_id)
all_functional_categories$unique_id<-gsub("(", ".", all_functional_categories$unique_id,fixed=TRUE)
all_functional_categories$unique_id<-gsub(")", ".", all_functional_categories$unique_id,fixed=TRUE)
all_functional_categories$unique_id<-gsub("=>", "..", all_functional_categories$unique_id)
all_functional_categories$unique_id<-gsub("-", ".", all_functional_categories$unique_id)
all_functional_categories$unique_id<-gsub(":", ".", all_functional_categories$unique_id)
all_functional_categories$unique_id<-gsub("/", ".", all_functional_categories$unique_id)


rownames(all_functional_categories) <- all_functional_categories$unique_id
all_functional_categories <- all_functional_categories[ -c(1:3,5,7) ]
conflicts_prefer(plyr::rename)

all_functional_categories<-rename(all_functional_categories, c("header" = "Specific Function" , "function_description" = "Broad Function"))
all_functional_categories_annotation <- all_functional_categories[ -c(2) ]



pheatmap(PA_functions_matrix, #annotation_row = Row_clusters, #annotation_col = all_functional_categories,
         show_colnames=FALSE,
         filename ="~/Desktop/heatmap5.pdf") #labels_row = Row_clusters$mOTU)



all_functions_colour = list(
  mOTU = c("1"="#000000", "2"="#444444", "3"="#888888", "4"="#F25F14", "5"="#F49027", "6"="#F7C13A", "7"="#30123B","8"="#296A8F","9"="#23C3E4","10"="#23C3E4", "11"="#5BFB72", "12"="#3FDFAB", "13"="#AFFA37", "14"="#CDE538", "15"="#EBD239", "16"="#332288", "17"="#6E3390", "18"="#AA4499", "19"="#970E01", "20"="#B72203", "21"="#D83706", "22"="#117733", "23"="#558833", "24"="#999933", "25"="#CC79A7", "26"="#E5A3C4", "27"="#FFCFE2", "28"="#00C2F9", "29"="#3EE0F9", "30"="#7CFFFA", "31"="#999999", "32"="#FFFFFF"),
  "Broad Function" = c("Carbon Utilization"="#CC79A7", "rRNA"="red","Organic Nitrogen"="sky blue","tRNA"="#1D91C0", "Energy"="yellow", "Misc"="orange","carbon utilization (Woodcroft)"="black","Transporters"="blue"),
  "Specific Function" = c("hydrocarbon degradation"="#E5F5F9","tRNA"="#1D91C0","Amino Acid"="#67001F","C1"="#F7FCFD","Photosynthesis"="#CB181D","Hydrogenases"="#78C679","anaerobic corrin ring synthesis"="#F46D43","CRISPR"="#A6CEE3","NA"="#FD8D3C","Metal Reduction"="#A6D854","sugar utilization (woodcroft)"="#D4B9DA","SCFA  and alcohol conversions"="#6A51A3","Nitrogen"="#7F0000","Antibiotic Resistance"="#D9D9D9","Flagellar cytoplasmic chaperone"="#FFF7BC","Peptidase"="#000000","pyruvate metabolism"="#F0F0F0","Information systems"="#C7EAE5","C1-methane"="#003C30","Sulfur"="#F16913","aerobic corrin ring synthesis"="#FFF7FB","Flagella Structure"="#8C6BB1","CAZY"="#C7E9B4","MISC"="#762A83","central carbon"="#FC9272","Electron transport Chain"="#AE017E","Oxygen"="#F7F7F7","ADO-CBL synthesis"="#DF65B0","Vitamin B12 transport system"="#EF3B2C"))

##to add colors for functions
#"Specific Function" = c("hydrocarbon degradation"="#543005",
# "tRNA"="#67001F",
# "Amino Acid"="#A6CEE3",
# "C1"="#F7FBFF",
# "Photosynthesis"="#FFF5EB",
# "Hydrogenases"="#FFF7F3",
# "anaerobic corrin ring synthesis"="#8C510A",
# "CRISPR"="#D73027",
# "NA"="#B3CDE3",
# "Metal Reduction"="#E5F5F9",
# "sugar utilization (woodcroft)"="#FEE8C8",
# "SCFA  and alcohol conversions"="#FEE0D2",
# "Nitrogen"="#DE77AE",
# "Antibiotic Resistance"="#F46D43",
# "Flagellar cytoplasmic chaperone"="#CBD5E8",
# "Peptidase"="#BFD3E6",
# "pyruvate metabolism"="#D0D1E6",
# "Information systems"="#D9F0A3",
# "C1-methane"="#C2A5CF",
# "Sulfur"="#FDAE61",
# "aerobic corrin ring synthesis"="#984EA3",
# "Flagella Structure"="#A8DDB5",
# "CAZY"="#A6BDDB",
# "MISC"="#7FCDBB",
# "central carbon"="#FEE0B6",
# "Electron transport Chain"="#386CB0",
# "Oxygen"="#A6D854",
# "ADO-CBL synthesis"="#74C476",
# "Vitamin B12 transport system"="#DF65B0")



pheatmap(PA_functions_matrix,
         color=c("red", "blue"),
         annotation_colors = all_functions_colour,
         #annotation_row = Row_clusters, 
         #annotation_col = all_functional_categories_annotation,
         labels_col = all_functional_categories$'Specific Function',
         show_colnames=TRUE,
         drop_levels = FALSE,
         annotation_legend= TRUE,
         annotation_names_row = TRUE,
         legend_breaks = c(1, 0),
         fontsize = 10, 
         fontsize_row = 6, 
         fontsize_col =3,
         cellheight = 5,
         angle_col = 45,
         cutree_col = 2,
         cutree_row = 4,
         filename ="/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/ETC_central_carbon_heatmap.pdf")

###Thermo Heatmaps cazy and ETC
#Functional Summary DRAM
#load in data for summary heatmap
product1 <- read.csv("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/1_genome_summaries/product.tsv",fill = TRUE, header = TRUE, sep = "\t",check.names=FALSE)
product2 <- read.csv("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/2_genome_summaries/product.tsv",fill = TRUE, header = TRUE, sep = "\t",check.names=FALSE)
product3 <- read.csv("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/3_genome_summaries/product.tsv",fill = TRUE, header = TRUE, sep = "\t",check.names=FALSE)
product4 <- read.csv("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/4_genome_summaries/product.tsv",fill = TRUE, header = TRUE, sep = "\t",check.names=FALSE)
product5 <- read.csv("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/5_genome_summaries/product.tsv",fill = TRUE, header = TRUE, sep = "\t",check.names=FALSE)
product6 <- read.csv("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/6_genome_summaries/product.tsv",fill = TRUE, header = TRUE, sep = "\t",check.names=FALSE)
product7 <- read.csv("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/7_genome_summaries/product.tsv",fill = TRUE, header = TRUE, sep = "\t",check.names=FALSE)
product8 <- read.csv("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/8_genome_summaries/product.tsv",fill = TRUE, header = TRUE, sep = "\t",check.names=FALSE)
product9 <- read.csv("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/9_genome_summaries/product.tsv",fill = TRUE, header = TRUE, sep = "\t",check.names=FALSE)

remove(product_combined)
product_combined <- dplyr::bind_rows(product1, product2, product3, product4, product5, product6, product7, product8, product9)

remove(mdata)

Row_clusters<-NMDS_metadata[ -c(2,5)]
Row_clusters<-distinct(Row_clusters,bin_ID,.keep_all= TRUE)
rownames(Row_clusters) <- Row_clusters$bin_ID
Row_clusters$Phylum<-as.character(Row_clusters$Phylum)
Row_clusters <- Row_clusters[ -c(1) ]

#Row_clusters<-rename(Row_clusters, c("ANI_Cluster" = "mOTU"))

#keep only bins being analysed
conflicts_prefer(rstatix::filter)
product_combined<-filter(product_combined, genome %in% All_bins$Bin_ID)


#split into presence/absence vs pathway completeness results
all_bin_pcompleteness <- product_combined[ -c(34:99) ]

all_bin_cazy <- product_combined[ -c(2:33) ]

#rename 'genome' to 'Bin_ID' for consistency
conflicts_prefer(plyr::rename)
all_bin_pcompleteness<-rename(all_bin_pcompleteness, c("genome" = "Bin_ID" ))
all_bin_cazy<-rename(all_bin_cazy, c("genome" = "Bin_ID"  ))

#remove duplicate rows/bins
all_bin_pcompleteness<-unique(all_bin_pcompleteness, by = "Bin_ID")
all_bin_cazy<-unique(all_bin_cazy, by = "Bin_ID")

#replace True/False for 1/0
all_bin_cazy[all_bin_cazy == "False"] <- "0"
all_bin_cazy[all_bin_cazy == "True"] <- "1"

#make functional data into matrix

rownames(all_bin_pcompleteness) <- all_bin_pcompleteness$Bin_ID

#identify and delete genes with no hits
all_bin_pcompleteness<-all_bin_pcompleteness[, colSums(all_bin_pcompleteness != 0) > 0]
pcompleteness_matrix <- as.matrix(all_bin_pcompleteness[,-1])

rownames(pcompleteness_matrix) <- all_bin_pcompleteness[,1]


completeness_matrix <- mapply(pcompleteness_matrix, FUN=as.numeric)
completeness_matrix <- matrix(data=completeness_matrix, ncol=length(colnames(pcompleteness_matrix)), nrow=length(row.names(pcompleteness_matrix)))
row.names(completeness_matrix) <- row.names(pcompleteness_matrix)
colnames(completeness_matrix) <- colnames(pcompleteness_matrix)

library(pheatmap)
pheatmap(completeness_matrix)


#import modified annotations
functional_categories <- read.csv("/Volumes/micro-shared$/MoralesLab/Projects/Syaliny_Soils/MAGs/data_analysis_by_paper/Gemm/raw_data/DRAM/genome_summaries/functional_catgories_rev.csv",fill = TRUE, header = TRUE, sep = ",")

#subset by those detected
func_list <- colnames(pcompleteness_matrix)
functional_categories<-filter(functional_categories, Function %in% func_list)

#Keep only two columns used for plots
rownames(functional_categories) <- functional_categories$Function
functional_categories <- functional_categories[ -c(1) ]
functional_categories<-rename(functional_categories, c("Renamed_function" = "Function"))
functional_categories_annotation <- functional_categories[ -c(2) ]



#test heatmap
pheatmap(completeness_matrix, annotation_row = Row_clusters, annotation_col = functional_categories,show_colnames=FALSE, labels_row = Row_clusters$Watermass)

func_category_colour = list(
  Phylum = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey"), "Category" = c("Central C metabolism"="#000000", "ETC complex I"="#990F0F", "ETC complex II"="#99540F", "ETC complex III"="#6B990F", "ETC complex IV High affinity"="#0F6B99", "ETC complex IV Low affinity"="#260F99", "ETC complex V"="#CC79A7"))
#set color palette

func_category_colour = list(
  Phylum = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey"),
  Watermass = c("Neritic" = "#0C0881", "STW" = "#7316A2", "Front" = "#BD5077", "SAW" = "#EA9953", "Deep" = "#F2F958"))

library(viridis)
pheatmap(completeness_matrix,
         color=turbo(10),
         annotation_colors = func_category_colour,
         annotation_row = Row_clusters, 
         annotation_col = functional_categories_annotation,
         labels_col = functional_categories$Function,
         show_colnames=TRUE,
         drop_levels = TRUE,
         annotation_legend= TRUE,
         annotation_names_row = TRUE,
         fontsize = 10, 
         fontsize_row = 6, 
         fontsize_col = 6,
         cellwidth = 8,
         angle_col = 45,
         cellheight = 5,
         cutree_rows = 6,
         filename ="/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/DRAM_MOTS_summary_heatmap.pdf")

dev.off()


#CAZY Summary DRAM
rownames(all_bin_cazy) <- all_bin_cazy$Bin_ID

#identify and delete genes with no hits
all_bin_cazy<-all_bin_cazy[, colSums(all_bin_cazy != 0) > 0]
cazy_matrix <- as.matrix(all_bin_cazy[,-1])

rownames(cazy_matrix) <- all_bin_cazy[,1]


bincazy_matrix <- mapply(cazy_matrix, FUN=as.numeric)
bincazy_matrix <- matrix(data=bincazy_matrix, ncol=length(colnames(cazy_matrix)), nrow=length(row.names(cazy_matrix)))
row.names(bincazy_matrix) <- row.names(cazy_matrix)
colnames(bincazy_matrix) <- colnames(cazy_matrix)

pheatmap(bincazy_matrix)

#create a document with functional annotations
cazy_categories<-as.data.frame(t(cazy_matrix))
cazy_categories <- tibble::rownames_to_column(cazy_categories, "Cazy_Function") # Apply rownames_to_column
#save as csv and modify file to create category annotations
#write.csv(cazy_categories, "/Volumes/micro-shared$/MoralesLab/Projects/Syaliny_Soils/MAGs/data_analysis_by_paper/Global_MAG_paper/raw_data_for_figs/cazy_categories.csv", row.names = T)

#re-import modified annotations
cazy_categories <- read.csv("/Volumes/micro-shared$/MoralesLab/Projects/Syaliny_Soils/MAGs/data_analysis_by_paper/Global_MAG_paper/raw_data_for_figs/cazy_categories_rev.csv",fill = TRUE, header = TRUE, sep = ",")

#remove leading space
cazy_categories$Renamed_function <- trimws(cazy_categories$Renamed_function, which = c("left"))

#subset by those detected
cazy_list <- colnames(cazy_matrix)
cazy_categories<-filter(cazy_categories, Cazy_Function %in% cazy_list)

#Keep only two columns used for plots
rownames(cazy_categories) <- cazy_categories$Cazy_Function
cazy_categories <- cazy_categories[ -c(1) ]
cazy_categories<-rename(cazy_categories, c("Renamed_function" = "CAZy_Function"))
cazy_categories_annotation <- cazy_categories[ -c(2) ]

#Row_clusters <- Row_clusters[ -c(3:4) ]





#test heatmap
pheatmap(bincazy_matrix, 
         annotation_row = Row_clusters, 
         annotation_col = cazy_categories,
         show_colnames=FALSE)
         #labels_row = Row_clusters$mOTU)


#set color palette
cazy_category_colour = list(
  Phylum = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey"),
  Category = c("CAZy"="#CC79A7", "Methanogenesis and methanotrophy"="#003C30", "Nitrogen metabolism"="#7F0000", "Other Reductases"="black", "SCFA and alcohol conversions"="#6A51A3", "Sulfur metabolism"="#FDAE61"),
  Watermass = c("Neritic" = "#0C0881", "STW" = "#7316A2", "Front" = "#BD5077", "SAW" = "#EA9953", "Deep" = "#F2F958"))

pheatmap(bincazy_matrix,
         color=c("red", "blue"),
         annotation_colors = cazy_category_colour,
         annotation_row = Row_clusters, 
         annotation_col = cazy_categories_annotation,
         labels_col = cazy_categories$CAZy_Function,
         show_colnames=TRUE,
         drop_levels = TRUE,
         annotation_legend= TRUE,
         annotation_names_row = TRUE,
         fontsize = 10, 
         fontsize_row = 6, 
         fontsize_col = 6,
         cellwidth = 9,
         angle_col = 45,
         cellheight = 5,
         cutree_rows = 6,
         filename ="/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/DRAM_MOTS_cazy_heatmap.pdf")

dev.off()

###PER GENUS
product1 <- read.csv("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/1_genome_summaries/product.tsv",fill = TRUE, header = TRUE, sep = "\t",check.names=FALSE)
product2 <- read.csv("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/2_genome_summaries/product.tsv",fill = TRUE, header = TRUE, sep = "\t",check.names=FALSE)
product3 <- read.csv("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/3_genome_summaries/product.tsv",fill = TRUE, header = TRUE, sep = "\t",check.names=FALSE)
product4 <- read.csv("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/4_genome_summaries/product.tsv",fill = TRUE, header = TRUE, sep = "\t",check.names=FALSE)
product5 <- read.csv("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/5_genome_summaries/product.tsv",fill = TRUE, header = TRUE, sep = "\t",check.names=FALSE)
product6 <- read.csv("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/6_genome_summaries/product.tsv",fill = TRUE, header = TRUE, sep = "\t",check.names=FALSE)
product7 <- read.csv("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/7_genome_summaries/product.tsv",fill = TRUE, header = TRUE, sep = "\t",check.names=FALSE)
product8 <- read.csv("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/8_genome_summaries/product.tsv",fill = TRUE, header = TRUE, sep = "\t",check.names=FALSE)
product9 <- read.csv("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/DRAM/9_genome_summaries/product.tsv",fill = TRUE, header = TRUE, sep = "\t",check.names=FALSE)

remove(product_combined)
product_combined <- dplyr::bind_rows(product1, product2, product3, product4, product5, product6, product7, product8, product9)

remove(mdata)

#Row_clusters<-NMDS_metadata[ -c(2,5)]
#Row_clusters<-distinct(Row_clusters,bin_ID,.keep_all= TRUE)
#rownames(Row_clusters) <- Row_clusters$bin_ID
#Row_clusters$Phylum<-as.character(Row_clusters$Phylum)
#Row_clusters <- Row_clusters[ -c(1) ]

mdata <- melt(All_bins, id=c("Bin_ID","Watermass", "Genus", "Phylum"))

##combining for just by genome
All_bins <- read.csv("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/raw_data/MOTS_bin_base_data.csv",fill = TRUE, header = TRUE, sep = ",")

All_bins_taxa <- as.matrix(All_bins[13:19])

All_bins_taxa <- tax_table(All_bins_taxa)


All_bins_taxa <- All_bins_taxa %>%
  tax_fix()

All_bins_taxa <- data.frame(All_bins_taxa)

All_bins <- All_bins %>%
  .[-c(13:19)]

All_bins_taxa$Bin_ID <- All_bins$Bin_ID

All_bins <- merge(All_bins, All_bins_taxa, by = c("Bin_ID"))

#All_bins <- column_to_rownames(All_bins, var = "Bin_ID")

mdata <- melt(All_bins, id=c("Bin_ID","Watermass", "Phylum", "Genus"))



##combining for just by genome
product_combined_genus <- merge(product_combined, mdata, by.x = ("genome"), by.y = ("Bin_ID"))
product_combined_genus <- subset(product_combined_genus, select = -c(genome, Watermass, Phylum, variable, value))

product_combined_genus <- product_combined_genus %>%
  select(Genus, everything())


Row_clusters<-mdata[ -c(5:6)]
Row_clusters<-distinct(Row_clusters,Genus,.keep_all= TRUE)
rownames(Row_clusters) <- Row_clusters$Genus
Row_clusters <- Row_clusters[ -c(1) ]
#Row_clusters$Treatment<-as.character(Row_clusters$Treatment)

Row_clusters <- Row_clusters[ -c(3)]


#keep only bins being analysed
library(conflicted)
conflicts_prefer(rstatix::filter)
product_combined<-filter(product_combined_genus, Genus %in% All_bins$Genus)




#split into presence/absence vs pathway completeness results
all_bin_pcompleteness <- product_combined[ -c(34:99) ]

all_bin_cazy <- product_combined[ -c(2:33) ]

#rename 'genome' to 'Bin_ID' for consistency
conflicts_prefer(plyr::rename)
all_bin_pcompleteness<-rename(all_bin_pcompleteness, c("Genus" = "Bin_ID" ))
all_bin_cazy<-rename(all_bin_cazy, c("Genus" = "Bin_ID"  ))



#remove duplicate rows/bins
all_bin_pcompleteness<-unique(all_bin_pcompleteness, by = "Bin_ID")
all_bin_cazy<-unique(all_bin_cazy, by = "Bin_ID")

#replace True/False for 1/0
all_bin_cazy[all_bin_cazy == "False"] <- "0"
all_bin_cazy[all_bin_cazy == "True"] <- "1"

#make functional data into matrix

all_bin_pcompleteness <- aggregate(. ~ Bin_ID, data = all_bin_pcompleteness, FUN = mean)

rownames(all_bin_pcompleteness) <- all_bin_pcompleteness$Bin_ID

#identify and delete genes with no hits
all_bin_pcompleteness<-all_bin_pcompleteness[, colSums(all_bin_pcompleteness != 0) > 0]
pcompleteness_matrix <- as.matrix(all_bin_pcompleteness[,-1])

rownames(pcompleteness_matrix) <- all_bin_pcompleteness[,1]


completeness_matrix <- mapply(pcompleteness_matrix, FUN=as.numeric)
completeness_matrix <- matrix(data=completeness_matrix, ncol=length(colnames(pcompleteness_matrix)), nrow=length(row.names(pcompleteness_matrix)))
row.names(completeness_matrix) <- row.names(pcompleteness_matrix)
colnames(completeness_matrix) <- colnames(pcompleteness_matrix)

library(pheatmap)
library(viridis)

pheatmap(completeness_matrix)


#import modified annotations
functional_categories <- read.csv("/Volumes/micro-shared$/MoralesLab/Projects/Syaliny_Soils/MAGs/data_analysis_by_paper/Gemm/raw_data/DRAM/genome_summaries/functional_catgories_rev.csv",fill = TRUE, header = TRUE, sep = ",")

#subset by those detected
func_list <- colnames(pcompleteness_matrix)
functional_categories<-filter(functional_categories, Function %in% func_list)

#Keep only two columns used for plots
rownames(functional_categories) <- functional_categories$Function
functional_categories <- functional_categories[ -c(1) ]
functional_categories<-rename(functional_categories, c("Renamed_function" = "Function"))
functional_categories_annotation <- functional_categories[ -c(2) ]



#test heatmap
pheatmap(completeness_matrix, annotation_row = Row_clusters, annotation_col = functional_categories,show_colnames=FALSE, labels_row = Row_clusters$Watermass)

#func_category_colour = list(
  #Phylum = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey"), "Category" = c("Central C metabolism"="#000000", "ETC complex I"="#990F0F", "ETC complex II"="#99540F", "ETC complex III"="#6B990F", "ETC complex IV High affinity"="#0F6B99", "ETC complex IV Low affinity"="#260F99", "ETC complex V"="#CC79A7"))
#set color palette

func_category_colour = list(
  Phylum = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey"),
  Watermass = c("NW" = "#0C0881", "STW" = "#7316A2", Front = "#BD5077", "SAW" = "#EA9953", "Deep" = "#F2F958"))

library(viridis)
library(pheatmap)
pheatmap(completeness_matrix,
         color=turbo(10),
         annotation_colors = func_category_colour,
         annotation_row = Row_clusters, 
         annotation_col = functional_categories_annotation,
         labels_col = functional_categories$Function,
         show_colnames=TRUE,
         drop_levels = TRUE,
         annotation_legend= TRUE,
         annotation_names_row = TRUE,
         fontsize = 10, 
         fontsize_row = 6, 
         fontsize_col = 6,
         cellwidth = 8,
         angle_col = 45,
         cellheight = 5,
         cutree_rows = 6,
         filename ="/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/MOTS_colourful_DRAM_summary_heatmap.pdf")

dev.off()

write.csv(all_bin_cazy, "/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/raw_data/all_bin_cazy.csv")

all_bin_cazy <- read.csv("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/raw_data/all_bin_cazy.csv", header = TRUE)

all_bin_cazy <- all_bin_cazy[-1]

all_bin_cazy <- aggregate(. ~ Bin_ID, all_bin_cazy, sum)

rownames(all_bin_cazy) <- all_bin_cazy$Bin_ID

all_bin_cazy_list <- all_bin_cazy$Bin_ID

#replace True/False for 1/0
all_bin_cazy[all_bin_cazy > 0] <- "1"
all_bin_cazy[all_bin_cazy == "0"] <- "0"



all_bin_cazy <- all_bin_cazy[-1]
all_bin_cazy$Bin_ID <- all_bin_cazy_list

all_bin_cazy = all_bin_cazy %>% dplyr::select("Bin_ID",  
                                              everything())

rownames(all_bin_cazy) <- all_bin_cazy$Bin_ID

#identify and delete genes with no hits
all_bin_cazy<-all_bin_cazy[, colSums(all_bin_cazy != 0) > 0]
cazy_matrix <- as.matrix(all_bin_cazy[,-1])

rownames(cazy_matrix) <- all_bin_cazy[,1]


bincazy_matrix <- mapply(cazy_matrix, FUN=as.numeric)
bincazy_matrix <- matrix(data=bincazy_matrix, ncol=length(colnames(cazy_matrix)), nrow=length(row.names(cazy_matrix)))
row.names(bincazy_matrix) <- row.names(cazy_matrix)
colnames(bincazy_matrix) <- colnames(cazy_matrix)

pheatmap(bincazy_matrix)

#create a document with functional annotations
cazy_categories<-as.data.frame(t(cazy_matrix))
cazy_categories <- tibble::rownames_to_column(cazy_categories, "Cazy_Function") # Apply rownames_to_column
#save as csv and modify file to create category annotations
write.csv(cazy_categories, "/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/raw_data/cazy_categories.csv", row.names = T)



#re-import modified annotations
cazy_categories <- read.csv("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/raw_data/cazy_categories_rev.csv",fill = TRUE, header = TRUE, sep = ",")
#cazy_categories <- cazy_categories[-1]
#remove leading space
cazy_categories$Renamed_function <- trimws(cazy_categories$Renamed_function, which = c("left"))

#subset by those detected


cazy_list <- colnames(all_bin_cazy)
cazy_list <- cazy_list[-1]
cazy_categories<-filter(cazy_categories, Cazy_Function %in% cazy_list)

#Keep only two columns used for plots
rownames(cazy_categories) <- cazy_categories$Cazy_Function
cazy_categories <- cazy_categories[ -c(1) ]
cazy_categories<-rename(cazy_categories, c("Renamed_function" = "CAZy_Function"))
cazy_categories_annotation <- cazy_categories[ -c(2) ]





#test heatmap
pheatmap(bincazy_matrix,
         annotation_col = cazy_categories,
         show_colnames=FALSE, 
         labels_row = Row_clusters$Genus)


#set color palette
cazy_category_colour = list(
  Phylum = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey"),
  Watermass = c("NW" = "#0C0881", "STW" = "#7316A2", Front = "#BD5077", "SAW" = "#EA9953", "Deep" = "#F2F958"),
  Category = c("CAZy"="#CC79A7", "Methanogenesis and methanotrophy"="#003C30", "Nitrogen metabolism"="#7F0000", "Photosynthesis"="#C2DD9B", "Other Reductases"="black", "SCFA and alcohol conversions"="#6A51A3", "Sulfur metabolism"="#FDAE61"))

pheatmap(bincazy_matrix,
         color=c("red", "blue"),
         annotation_colors = cazy_category_colour,
         annotation_row = Row_clusters, 
         annotation_col = cazy_categories_annotation,
         labels_col = cazy_categories$CAZy_Function,
         show_colnames=TRUE,
         drop_levels = TRUE,
         annotation_legend= TRUE,
         annotation_names_row = TRUE,
         fontsize = 10, 
         fontsize_row = 6, 
         fontsize_col = 6,
         cellwidth = 9,
         angle_col = 45,
         cellheight = 5,
         cutree_rows = 6,
         filename ="/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/DRAM_MOTS_cazy_heatmap_genus.pdf")

dev.off()



#CAZY Summary DRAM
rownames(all_bin_cazy) <- all_bin_cazy$Bin_ID

#identify and delete genes with no hits
all_bin_cazy<-all_bin_cazy[, colSums(all_bin_cazy != 0) > 0]
cazy_matrix <- as.matrix(all_bin_cazy[,-1])

rownames(cazy_matrix) <- all_bin_cazy[,1]


bincazy_matrix <- mapply(cazy_matrix, FUN=as.numeric)
bincazy_matrix <- matrix(data=bincazy_matrix, ncol=length(colnames(cazy_matrix)), nrow=length(row.names(cazy_matrix)))
row.names(bincazy_matrix) <- row.names(cazy_matrix)
colnames(bincazy_matrix) <- colnames(cazy_matrix)

pheatmap(bincazy_matrix)

#create a document with functional annotations
cazy_categories<-as.data.frame(t(cazy_matrix))
cazy_categories <- tibble::rownames_to_column(cazy_categories, "Cazy_Function") # Apply rownames_to_column
#save as csv and modify file to create category annotations
#write.csv(cazy_categories, "/Volumes/micro-shared$/MoralesLab/Projects/Syaliny_Soils/MAGs/data_analysis_by_paper/Global_MAG_paper/raw_data_for_figs/cazy_categories.csv", row.names = T)

#re-import modified annotations
cazy_categories <- read.csv("/Volumes/micro-shared$/MoralesLab/Projects/Syaliny_Soils/MAGs/data_analysis_by_paper/Global_MAG_paper/raw_data_for_figs/cazy_categories_rev.csv",fill = TRUE, header = TRUE, sep = ",")

#remove leading space
cazy_categories$Renamed_function <- trimws(cazy_categories$Renamed_function, which = c("left"))

#subset by those detected
cazy_list <- colnames(cazy_matrix)
cazy_categories<-filter(cazy_categories, Cazy_Function %in% cazy_list)

#Keep only two columns used for plots
rownames(cazy_categories) <- cazy_categories$Cazy_Function
cazy_categories <- cazy_categories[ -c(1) ]
cazy_categories<-rename(cazy_categories, c("Renamed_function" = "CAZy_Function"))
cazy_categories_annotation <- cazy_categories[ -c(2) ]

Row_clusters <- Row_clusters[ -c(3:4) ]





#test heatmap
pheatmap(bincazy_matrix, 
         annotation_row = Row_clusters, 
         annotation_col = cazy_categories,
         show_colnames=FALSE)
#labels_row = Row_clusters$mOTU)


#set color palette
cazy_category_colour = list(
  Phylum = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey"),
  "Category" = c("CAZy"="#CC79A7", "Methanogenesis and methanotrophy"="#003C30", "Nitrogen metabolism"="#7F0000", "Other Reductases"="black", "SCFA and alcohol conversions"="#6A51A3", "Sulfur metabolism"="#FDAE61"),
  Watermass = c("NW" = "#0C0881", "STW" = "#7316A2", "Front" = "#BD5077", "SAW" = "#EA9953", "Deep" = "#F2F958"))

pheatmap(bincazy_matrix,
         color=c("red", "blue"),
         annotation_colors = cazy_category_colour,
         annotation_row = Row_clusters, 
         annotation_col = cazy_categories_annotation,
         labels_col = cazy_categories$CAZy_Function,
         show_colnames=TRUE,
         drop_levels = TRUE,
         annotation_legend= TRUE,
         annotation_names_row = TRUE,
         fontsize = 10, 
         fontsize_row = 6, 
         fontsize_col = 6,
         cellwidth = 9,
         angle_col = 45,
         cellheight = 5,
         cutree_rows = 6,
         filename ="/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/DRAM_cazy_heatmap.pdf")

dev.off()

pheatmap(bincazy_matrix,
         color=c("red", "blue"),
         annotation_colors = cazy_category_colour,
         annotation_row = Row_clusters, 
         annotation_col = cazy_categories_annotation,
         labels_col = cazy_categories$CAZy_Function,
         show_colnames=TRUE,
         drop_levels = TRUE,
         annotation_legend= TRUE,
         annotation_names_row = TRUE,
         fontsize = 10, 
         fontsize_row = 6, 
         fontsize_col = 6,
         cellwidth = 9,
         angle_col = 45,
         cellheight = 5,
         cutree_rows = 6,
         filename ="/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/DRAM_cazy_heatmap_Treatment.pdf")

dev.off()


###GGPLOT FUNCTION
metabolism_summary_for_central_ETC_counted

ggplot(metabolism_summary_for_central_ETC_counted, aes(gene_description, bin, fill= counted)) + 
  geom_tile()

