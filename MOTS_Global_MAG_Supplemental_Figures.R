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

#Created by Cecilia Wang. Make function to collect info on classification per rank

#for (i in 1:nrow(mots.taxa.csv)){
#if (mots.taxa.csv[i,2] == ""){
#kingdom <- paste("Kingdom_", mots.taxa.csv[i,1], sep = "")
#mots.taxa.csv[i, 2:7] <- kingdom
#} else if (mots.taxa.csv[i,3] == ""){
#phylum <- paste("Phylum_", mots.taxa.csv[i,2], sep = "")
#mots.taxa.csv[i, 3:7] <- phylum
#} else if (mots.taxa.csv[i,4] == ""){
#class <- paste("Class_", mots.taxa.csv[i,3], sep = "")
#mots.taxa.csv[i, 4:7] <- class
#} else if (mots.taxa.csv[i,5] == ""){
#order <- paste("Order_", mots.taxa.csv[i,4], sep = "")
#mots.taxa.csv[i, 5:7] <- order
#} else if (mots.taxa.csv[i,6] == ""){
#family <- paste("Family_", mots.taxa.csv[i,5], sep = "")
#mots.taxa.csv[i, 6:7] <- family
#} else if (mots.taxa.csv[i,7] == ""){
#mots.taxa.csv$Species[i] <- paste("Genus", mots.taxa.csv$Genus[i], sep = "_")
#}
#}


rank_df <- mots.taxa.csv[ -c(8:15) ]

unique_bin_across_rank<-function(df){
  rank_sum<-NULL # create an empty variable to save results into a table
  for (r in colnames(df)) {
    # print(r) # sanity check
    t<-data.frame(rank_df[r],rank_df[,"Phylum"]) %>% rownames_to_column() # collect info of a certain rank name, the phylum name, and the bin ID.
    temp_df<-t %>% dplyr::group_by(t[,c(2,3)]) %>% dplyr::summarise(bin_count=length(unique(rowname))) # group by phylum and another rank, then count the number of unique bins
    colnames(temp_df)<-c("taxa_name","Phylum","bin_count") # change column names so the temp dataframe can be combined
    temp_df$rank<-r # add a column to record the rank name
    rank_sum<-rbind(temp_df,rank_sum) # combine all temp data
  }
  return(rank_sum) # export the results
}

rank_df_less_watermass <- rank_df[-8]

unique_bin_across_rank_summary<-unique_bin_across_rank(rank_df_less_watermass)

Phylum_rank_bin_summary<-unique_bin_across_rank_summary[unique_bin_across_rank_summary$taxa_name!="",] %>% dplyr::group_by(Phylum,rank) %>% dplyr::summarise(unique_bin_count=sum(bin_count),unique_taxa=length(unique(taxa_name)))

unique_bin_across_rank_summary_phylum_sum<-subset(unique_bin_across_rank_summary, unique_bin_across_rank_summary$rank=="Phylum")

Phylum_rank_bin_summary$total_bin_by_phyla<-unique_bin_across_rank_summary_phylum_sum$bin_count[match(Phylum_rank_bin_summary$Phylum,unique_bin_across_rank_summary_phylum_sum$Phylum)]


conflict_prefer(ggpubr::mutate)
Phylum_rank_bin_summary <- Phylum_rank_bin_summary%>%
  dplyr::rename(Phylum = 1,
                Rank = 2,
                Number_of_classified_bins = 3,
                Number_of_unique_groups = 4)%>%
  ggpubr::mutate(Rank=factor(Rank), Phylum=factor(Phylum)) %>% 
  ggpubr::mutate(Rank=fct_relevel(Rank,c("Domain","Phylum","Class","Order", "Family", "Genus","Species")))
Phylum_rank_bin_summary$classified_bin_percent<-Phylum_rank_bin_summary$Number_of_classified_bins/Phylum_rank_bin_summary$total_bin_by_phyla*100 

library(viridis)
ggplot(Phylum_rank_bin_summary,aes(x= Rank, y= Number_of_unique_groups, size=classified_bin_percent)) + 
  facet_wrap(~Phylum)+
  geom_point(stat="identity")+
  labs(y="Number of Distinct Taxa", x="Taxonomic Rank", size = "Percent of classified MAGS")+
  theme_light()+
  theme(axis.title.x = element_text(face="bold",size=12),
        axis.text.x = element_text(angle=45, colour = "black", size=8,vjust = 0.5, hjust=0.5), 
        axis.text.y = element_text(colour = "black", size=10),
        axis.title.y = element_text(face="bold",size=10),
        plot.title = element_text(size = 12),
        legend.title =element_text(face="bold",size = 12),
        legend.text = element_text(size = 10),
        legend.position="right",
        legend.key.size = unit(1, "cm"),
        panel.background = element_blank(),
        strip.text.x = element_text(size=12, color="black"),
        panel.border = element_rect(fill = NA, colour = "black"))+
  scale_size_continuous(limits = c(0, 100),breaks=c(1,10,100))+
  scale_y_continuous(limits = c(0, 100),breaks=c(1,20,40,60,80,100))


#reorder Phyla by highest number of bins and replot
Phylum_rank_bin_summary$Phylum<-factor(Phylum_rank_bin_summary$Phylum,levels =c( "Proteobacteria","Bacteroidota","Thermoplasmatota","Verrucomicrobiota","Actinobacteriota","Cyanobacteria","Chloroflexota","Crenarchaeota","Planctomycetota","SAR324","Gemmatimonadota","Acidobacteriota","Nitrospinota","Myxococcota","Binatota","Dadabacteria","Fibrobacterota","Marinisomatota") )

ggplot(Phylum_rank_bin_summary,aes(x= Rank, y= Number_of_unique_groups, size=classified_bin_percent)) + 
  facet_wrap(~Phylum)+
  geom_point(stat="identity")+
  labs(y="Number of Distinct Taxa", x="Taxonomic Rank", size = "Percent of classified MAGS") +
  theme(axis.title.x = element_text(face="bold",size=20),
        axis.text.x = element_text(angle=45, colour = "black", size=14,vjust = 0.5, hjust=0.5), 
        axis.text.y = element_text(colour = "black", size=14),
        axis.title.y = element_text(face="bold",size=14),
        plot.title = element_text(size = 12),
        legend.title =element_text(face="bold",size = 14),
        legend.text = element_text(size = 14),
        legend.position="right",
        legend.key.size = unit(1, "cm"),
        panel.background = element_blank(),
        strip.text.x = element_text(size=14, color="black"),
        panel.border = element_rect(fill = NA, colour = "black"))+
  scale_size_continuous(limits = c(0, 100),breaks=c(1,10,100))+
  scale_y_continuous(limits = c(0, 100),breaks=c(1,20,40,60,80,100))


#theme_light()+
#theme(axis.title.x = element_text(face="bold",size=12),
axis.text.x = element_text(angle=45, colour = "black", size=8,vjust = 0.5, hjust=0.5), 
axis.text.y = element_text(colour = "black", size=10),
axis.title.y = element_text(face="bold",size=10),
plot.title = element_text(size = 12),
legend.title =element_text(face="bold",size = 12),
legend.text = element_text(size = 10),
legend.position="right",
legend.key.size = unit(1, "cm"),
panel.background = element_blank(),
strip.text.x = element_text(size=12, color="black"),
panel.border = element_rect(fill = NA, colour = "black"))

#stacked barplot by watermass and phylum 
#Counting phylum

conflicts_prefer(dplyr::count)

rank_df <- read.csv("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/read_recruitment_to_bins/import/rank_df.csv", header = TRUE)
rank_df_phylum_counts <- dplyr::rename(count(rank_df, Phylum, Watermass), Freq = n)
rank_df_phylum_counts <- rank_df_phylum_counts[!grepl("Neritic", rank_df_phylum_counts$Watermass),]
rank_df_phylum_counts <- rank_df_phylum_counts[!grepl("Front", rank_df_phylum_counts$Watermass),]


#counting by order
rank_df_order_counts <- dplyr::rename(count(rank_df, Order, Watermass), Freq = n)
rank_df_order_counts <- rank_df_order_counts[!grepl("Neritic", rank_df_order_counts$Watermass),]
rank_df_order_counts <- rank_df_order_counts[!grepl("Front", rank_df_order_counts$Watermass),]

#counting by family
rank_df_family_counts <- dplyr::rename(count(rank_df, Family, Watermass), Freq = n)
rank_df_family_counts <- rank_df_family_counts[!grepl("Neritic", rank_df_family_counts$Watermass),]
rank_df_family_counts <- rank_df_family_counts[!grepl("Front", rank_df_family_counts$Watermass),]

library("viridis")



##order plot
ggplot(rank_df_order_counts, aes(x=factor(Watermass, c("STW","SAW","Deep")), Watermass, y=Freq, fill = Order)) +
  stat_summary(geom = "bar", fun = "sum", position = "stack") +
  theme_minimal() +
  labs(x = "Watermass", y = "Total number of detected MAGS")
##family plot
ggplot(rank_df_family_counts, aes(x=factor(Watermass, c("STW","SAW","Deep")), Watermass, y=Freq, fill = Family)) +
  stat_summary(geom = "bar", fun = "sum", position = "stack") +
  theme_minimal() +
  labs(x = "Watermass", y = "Total number of detected MAGS")

##FOR STW
MOTS_bin_base_data <- read.csv("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/raw_data/MOTS_bin_func_long.csv", header = TRUE)

conflicts_prefer(dplyr::summarise)

counts <- MOTS_bin_base_data %>%
  group_by(Watermass, Phylum, Gene, Category) %>%
  summarise(across(c(Counts), sum))

completeness <- MOTS_bin_base_data %>%
  group_by(Watermass, Phylum, Gene, Category) %>%
  summarise(across(c(Completeness), mean))

#MOTS_bin_base_data <- merge(counts$Counts, completeness$Completeness)
MOTS_bin_base_data <- cbind(counts, completeness$Completeness)
names(MOTS_bin_base_data)[names(MOTS_bin_base_data) == '...6'] <- 'Completeness'


MOTS_bin_func_STW_long = subset(MOTS_bin_base_data, Watermass == "STW")
MOTS_bin_func_SAW_long = subset(MOTS_bin_base_data, Watermass == "SAW")
MOTS_bin_func_Deep_long = subset(MOTS_bin_base_data, Watermass == "Deep")



#Remove count data that are showing up as NA
MOTS_bin_func_STW_long=na.omit(MOTS_bin_func_STW_long)

#Order phylum ascending
MOTS_bin_func_STW_long$Phylum=factor(MOTS_bin_func_STW_long$Phylum, levels=c("Dadabacteria","Binatota","Fibrobacterota","Marinisomatota","Myxococcota","Nitrospinota","Acidobacteriota","Gemmatimonadota","Planctomycetota","SAR324",                                            "Crenarchaeota","Chloroflexota","Cyanobacteria","Verrucomicrobiota","Actinobacteriota","Bacteroidota",
                                                                             "Thermoplasmatota","Proteobacteria"))


library(tidyverse)
library(grid)
library(cowplot)


p2_STW <- ggplot(data = MOTS_bin_func_STW_long, aes(x = Phylum, y = Gene, size=Counts))+
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

p2_STW

p2_STW=p2_STW+theme(panel.spacing =unit(.05, "lines"),
                    panel.border = element_rect(color = "black", fill = NA, size = .5), 
                    strip.background = element_rect(color = "white", size = .5))

g2_STW <- ggplot_gtable(ggplot_build(p2_STW))
stripr <- which(grepl('strip-r', g2_STW$layout$name))
fills <- c("cadetblue3","aquamarine2","darkgoldenrod1","lightcoral","thistle3","rosybrown4","lightsteelblue3","yellow")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g2_STW$grobs[[i]]$grobs[[1]]$childrenOrder))
  g2_STW$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

grid.draw(g2_STW)



func_col <- c("Alternative e acceptor"= "cadetblue3",
              "Alternative e donor"="aquamarine2",
              "Carbon fixation"="darkgoldenrod1",
              "Nitrogen cycle"="lightcoral",
              "Phototrophy"="thistle3",
              "Respiration"="rosybrown4",
              "Sulfur cycle"="lightsteelblue3",
              "Trace gas metabolism
"="yellow")


legend_func<-ggplot(MOTS_bin_func_STW_long, aes(x=Category, fill=Phylum)) +
  geom_bar(stat="count")+
  scale_fill_manual(values = func_col)+
  My_Theme+
  theme(legend.position="bottom",
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"))+
  guides(fill=guide_legend(nrow=2))


#Save the legend  
func_legend <- get_legend(legend_func)
grid.newpage()
grid.draw(func_legend)




p2_SAW <- ggplot(data = MOTS_bin_func_SAW_long, aes(x = Phylum, y = Gene, size=Counts))+
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

p2_SAW

p2_SAW=p2_SAW+theme(panel.spacing =unit(.05, "lines"),
                    panel.border = element_rect(color = "black", fill = NA, size = .5), 
                    strip.background = element_rect(color = "white", size = .5))

g2_SAW <- ggplot_gtable(ggplot_build(p2_SAW))
stripr <- which(grepl('strip-r', g2_SAW$layout$name))
fills <- c("cadetblue3","aquamarine2","darkgoldenrod1","lightcoral","thistle3","rosybrown4","lightsteelblue3","yellow")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g2_SAW$grobs[[i]]$grobs[[1]]$childrenOrder))
  g2_SAW$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

grid.draw(g2_SAW)



func_col <- c("Alternative e acceptor"= "cadetblue3",
              "Alternative e donor"="aquamarine2",
              "Carbon fixation"="darkgoldenrod1",
              "Nitrogen cycle"="lightcoral",
              "Phototrophy"="thistle3",
              "Respiration"="rosybrown4",
              "Sulfur cycle"="lightsteelblue3",
              "Trace gas metabolism
"="yellow")


legend_func<-ggplot(MOTS_bin_func_STW_long, aes(x=Category, fill=Phylum)) +
  geom_bar(stat="count")+
  scale_fill_manual(values = func_col)+
  My_Theme+
  theme(legend.position="bottom",
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"))+
  guides(fill=guide_legend(nrow=2))


#Save the legend  
func_legend <- get_legend(legend_func)
grid.newpage()
grid.draw(func_legend)



p2_Deep <- ggplot(data = MOTS_bin_func_Deep_long, aes(x = Phylum, y = Gene, size=Counts))+
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

p2_Deep

p2_Deep=p2_Deep + My_Theme
#+theme(panel.spacing =unit(.05, "lines"),
panel.border = element_rect(color = "black", fill = NA, size = .5)
strip.background = element_rect(color = "white", size = .5)

g2_Deep <- ggplot_gtable(ggplot_build(p2_Deep))
stripr <- which(grepl('strip-r', g2_Deep$layout$name))
fills <- c("cadetblue3","aquamarine2","darkgoldenrod1","lightcoral","thistle3","rosybrown4","lightsteelblue3","yellow")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g2_Deep$grobs[[i]]$grobs[[1]]$childrenOrder))
  g2_Deep$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

grid.draw(g2_Deep)



func_col <- c("Alternative e acceptor"= "cadetblue3",
              "Alternative e donor"="aquamarine2",
              "Carbon fixation"="darkgoldenrod1",
              "Nitrogen cycle"="lightcoral",
              "Phototrophy"="thistle3",
              "Respiration"="rosybrown4",
              "Sulfur cycle"="lightsteelblue3",
              "Trace gas metabolism
"="yellow")


legend_func<-ggplot(MOTS_bin_func_STW_long, aes(x=Category, fill=Phylum)) +
  geom_bar(stat="count")+
  scale_fill_manual(values = func_col)+
  My_Theme+
  theme(legend.position="bottom",
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"))+
  guides(fill=guide_legend(nrow=2))


#Save the legend  
func_legend <- get_legend(legend_func)
grid.newpage()
grid.draw(func_legend)


pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/STW_diamond_plot.pdf", width = 8, height = 10) # Open a new pdf file
grid.arrange(g2_STW, func_legend, ncol=1, nrow=2, widths=c(8), heights=c(10, 1))
dev.off()

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/SAW_diamond_plot.pdf", width = 8, height = 10) # Open a new pdf file
grid.arrange(g2_SAW, func_legend, ncol=1, nrow=2, widths=c(8), heights=c(10, 1))
dev.off()

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/Deep_diamond_plot.pdf", width = 8, height = 10) # Open a new pdf file
grid.arrange(g2_Deep, func_legend, ncol=1, nrow=2, widths=c(8), heights=c(10, 1))

dev.off()


##NMDS plot using total functional data DRAM
##plot with watermass as colours
library(vegan)
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

metabolism_summary<-dplyr::bind_rows(MISC_combined,carbon_utilization_combined,Transporters_combined,Energy_combined,organic_nitrogen_combined,Woodcroft_combined,rRNA_combined)
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

with_season_meta <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/read_recruitment_to_bins/Yugen_working_files/Nova_metadata_season.xlsx")

with_season_meta <- merge(x=NMDS_metadata, y=with_season_metadata, by.x='Sample', by.y='Sample')
with_season_meta <- subset(with_season_meta, select = -c(Phylum))

mots.taxa.csv$bin_ID <- rownames(mots.taxa.csv)

mots_taxa <- mots.taxa.csv
mots_taxa <- subset(mots_taxa, select = c(bin_ID, Phylum))

with_season_meta <- merge(x=with_season_meta, y = mots_taxa, by = "bin_ID")



with_season_meta <- subset(with_season_meta, select = c(bin_ID, Phylum, Watermass, Sample, Season, cluster, Year))

data.scores <- merge(x=with_season_meta, y=data.scores, by.x='bin_ID', by.y='bin_ID')

#data.scores <- cbind(with_season_meta, data.scores)

head(data.scores)

#data.scores <- subset(data.scores, select = c(bin_ID, Watermass.y, Phylum.y, Season, cluster.x, NMDS1, NMDS2))

rownames(data.scores ) <- NULL

#data.scores_clust <- read.csv("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/read_recruitment_to_bins/Yugen_working_files/data.scores_clust.csv")

NMDS_watermass <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(color = factor(Watermass, level=c('Neritic','STW','Front','SAW','Deep')))) + 
  labs(x = "NMDS1", colour = "Time", y = "NMDS2", shape = "Type", tags = "A")  + 
  scale_color_viridis(discrete=TRUE, option = "C") +
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
        panel.background = element_blank())

NMDS_year <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(color = factor(Year, level=c('2014','2015','2016','2017', '2018')))) + 
  labs(x = "NMDS1", colour = "Time", y = "NMDS2", shape = "Type", tags = "C")  + 
  scale_color_viridis(discrete=TRUE, option = "C") +
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
        panel.background = element_blank())

NMDS_season <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(color = factor(Season, level=c('Autumn','Winter','Spring','Summer')))) + 
  labs(x = "NMDS1", colour = "Time", y = "NMDS2", shape = "Type" , tags = "B")  + 
  scale_color_viridis(discrete=TRUE, option = "C") +
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
        panel.background = element_blank())

pdf("~/Desktop/NMDS_watermass.pdf",width=12,height=10) # Open a new pdf file

#view
NMDS_watermass

dev.off()

pdf("~/Desktop/NMDS_year.pdf",width=12,height=10) # Open a new pdf file

#view
NMDS_year

dev.off()

pdf("~/Desktop/NMDS_season.pdf",width=12,height=10) # Open a new pdf file

#view
NMDS_season

dev.off()

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/grid_NMDS.pdf", width = 15 , height = 15)
grid.arrange(NMDS_watermass, NMDS_season, NMDS_year,
             ncol = 2)
dev.off()


#data.scores$Phylum.y<-factor(data.scores$Phylum.y,levels =c( "Proteobacteria","Bacteroidota","Thermoplasmatota","Verrucomicrobiota","Actinobacteriota","Cyanobacteria","Chloroflexota","Crenarchaeota","Planctomycetota","SAR324","Gemmatimonadota","Acidobacteriota","Nitrospinota","Myxococcota","Binatota","Dadabacteria","Fibrobacterota","Marinisomatota") )

data.scores$Phylum<-factor(data.scores$Phylum,levels =c( "Proteobacteria","Bacteroidota","Thermoplasmatota","Verrucomicrobiota","Actinobacteriota","Cyanobacteria","Chloroflexota","Crenarchaeota","Planctomycetota","SAR324","Gemmatimonadota","Acidobacteriota","Nitrospinota","Myxococcota","Binatota","Dadabacteria","Fibrobacterota","Marinisomatota") )

data.scores$cluster <- str_replace(data.scores$cluster, "clust1", "Cluster 1")
data.scores$cluster <- str_replace(data.scores$cluster, "clust2", "Cluster 2")
data.scores$cluster <- str_replace(data.scores$cluster, "clust3", "Cluster 3")

functional_NMDS_Phylum <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(color = Phylum)) +
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



pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/functional_NMDS_Phylum.pdf",width=12,height=10) # Open a new pdf file

#view
functional_NMDS_Phylum

dev.off()



anosim(m_com, data.scores$Phylum, distance = "bray", permutations = 9999)

adonis2(m_com ~ cluster, data=data.scores, perm=999)

## Bray-Curtis distances between samples
dis <- vegdist(m_com)

## Calculate multivariate dispersions
mod <- betadisper(dis, Watermass)
mod

#ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
#geom_point(size = 4, aes(color = Phylum)) +
# scale_color_manual(values = c("Chloroflexota" = "#B4CCDF", "Cyanobacteria" = "#4A75AA", "Dadabacteria" = "#C2DD9B", "Fibrobacterota" = "#669C4C", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#C3B4D1", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C")) +
# My_Theme +
# labs(x = "NMDS1", colour = "Phylum", y = "NMDS2", shape = "Type")


input_season_cluster <- subset(data.scores, select = c(Season, cluster))
input_season_cluster <- dplyr::rename(count(input_season_cluster, Season, cluster), Freq = n)

cluster_count <- data.frame(cluster = c("clust1", "clust2", "clust3"),
                            clust_count = c(26, 122, 858))

cluster_count <- merge(input_season_cluster, cluster_count, by = "cluster")

input_season_cluster$norm <- cluster_count$Freq/cluster_count$clust_count

ggplot(input_season_cluster, aes(x = factor(Season, level=c('Autumn','Winter','Spring','Summer')), y = norm, group = cluster, color = cluster)) +
  geom_line(stat = "identity", size = 1.5) +
  scale_color_brewer(palette="Paired") +
  labs(x = "Season", y = "Normalized count") +
  My_Theme

kruskal.test(norm ~ cluster, data = input_season_cluster)
kruskal.test(Freq ~ Season, data = input_season_cluster)

nest <- aov(input_season_cluster$norm ~ input_season_cluster$Season / factor(input_season_cluster$clust))



##add in season with abundance
###data.scores_env clust
abundance <- psmelt(mots.phylo)


data.scores_env_clust <- merge(x = data.scores_env, y = data.scores_clust, by.x='bin_ID', by.y='bin_ID')
data.scores_env_clust_season <- merge (x = data.scores_env_clust, y = with_season_metadata, by.x='Sample.x', by.y='Sample')

data.scores_env_clust_season <- subset(data.scores_env_clust_season, select = c(bin_ID, clust, Season))

abundance_clust <- merge(x = abundance, y = data.scores_env_clust_season, by.x='OTU', by.y='bin_ID')

cluster_count <- data.frame(cluster = c("Cluster 1", "Cluster 2", "Cluster 3"),
                            clust_count = c(26, 122, 858))

abundance_clust <- merge(abundance_clust, cluster_count, by.x = "clust", by.y="cluster")

abundance_clust$norm <- abundance_clust$Abundance/abundance_clust$clust_count

abundance_clust_subset <- subset(abundance_clust, select = c(clust, Watermass, Phylum, Season.y, norm))


abundance_merged <- abundance_clust_subset %>%
  group_by(clust, Season.y) %>%
  summarise(across(c(norm), sum))

ggplot(abundance_merged, aes(x = factor(Season.y, level=c('Spring','Summer','Autumn','Winter')), y = norm, group = clust, color = clust)) +
  geom_line(stat = "identity", size = 1.5) +
  scale_color_brewer(palette="Paired") +
  labs(x = "Season", y = "Normalized change in absolute abundance") +
  theme(axis.title.x = element_text(face="bold",size=24),
        axis.text.x = element_text(angle = 90, colour = "black", size=22), 
        axis.text.y = element_text(colour = "black", size=22),
        axis.title.y = element_text(face="bold",size=24),
        plot.title = element_text(size = 24),
        legend.title =element_text(face="bold",size = 14),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1, "cm"),
        strip.text.x = element_text(size=16, face="bold"),
        strip.text.y = element_text(size=22, face="bold"),
        panel.border = element_rect(fill = NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

abundance_merged <- abundance_clust %>%
  group_by(clust, Phylum,Season.y) %>%
  summarise(across(c(norm), sum))

ggplot(abundance_merged, aes(x = factor(Season.y, level=c('Spring','Summer','Autumn','Winter')), y = norm, group = clust, color = clust)) +
  geom_line(stat = "identity", size = 1.5) +
  facet_wrap(~ Phylum, ncol = 6) +
  scale_color_brewer(palette="Paired") +
  labs(x = "Season", y = "Normalized change in absolute abundance") +
  theme(axis.title.x = element_text(face="bold",size=24),
        axis.text.x = element_text(angle = 90, colour = "black", size=22), 
        axis.text.y = element_text(colour = "black", size=22),
        axis.title.y = element_text(face="bold",size=24),
        plot.title = element_text(size = 24),
        legend.title =element_text(face="bold",size = 14),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1, "cm"),
        strip.text.x = element_text(size=16, face="bold"),
        strip.text.y = element_text(size=22, face="bold"),
        panel.border = element_rect(fill = NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

nest <- aov(abundance_merged$norm ~ abundance_merged$Season.y / factor(abundance_merged$clust))
summary(nest)

###data.scores_env clust

data.scores_env_clust <- merge(x = data.scores_env, y = data.scores_clust, by.x='bin_ID', by.y='bin_ID')

ggplot(data.scores_env_clust, aes(x = Phylum.x, y = Temp)) +
  geom_boxplot() +
  facet_grid(~clust) +
  My_Theme


###NMDS by Tranpsorters function only
transporters_metabolism_summary<-organic_nitrogen_combined

t_transporters_metabolism_summary <- t(transporters_metabolism_summary)

transporters_com = transporters_metabolism_summary[,7:ncol(transporters_metabolism_summary)]

transporters_com <- t(transporters_com)
transporters_m_com <- as.matrix(transporters_com)


transporters_tot <- rowSums(transporters_m_com)
transporters_m_com <- transporters_m_com[tot > 0, ]

transporters_m_com <- na.omit(transporters_m_com)

set.seed(123)
transporters_nmds = metaMDS(transporters_m_com, distance = "bray")
transporters_nmds

transporters_data.scores = as.data.frame(scores(transporters_nmds)$sites)

#add columns to data frame 
NMDS_metadata <- read.csv("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/read_recruitment_to_bins/Yugen_working_files/for_NMDS_metadata.csv")



transporters_data.scores <- cbind(NMDS_metadata, transporters_data.scores)

head(transporters_data.scores)

transporters_data.scores <- subset(transporters_data.scores, select = c(bin_ID, Sample, Watermass, Phylum, cluster, NMDS1, NMDS2))

rownames(data.scores ) <- NULL

transporters_data.scores_clust <- read.csv("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/read_recruitment_to_bins/Yugen_working_files/data.scores_clust.csv")

ggplot(transporters_data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(color = factor(Watermass, level=c('Neritic','STW','Front','SAW','Deep')))) + 
  My_Theme +
  labs(x = "NMDS1", colour = "Time", y = "NMDS2", shape = "Type")  + 
  scale_color_viridis(discrete=TRUE, option = "C")

ggplot(transporters_data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(color = factor(Year, level=c('2014','2015','2016','2017', '2018')))) + 
  My_Theme +
  labs(x = "NMDS1", colour = "Time", y = "NMDS2", shape = "Type")  + 
  scale_color_viridis(discrete=TRUE, option = "C")



transporters_data.scores$Phylum<-factor(data.scores$Phylum,levels =c( "Proteobacteria","Bacteroidota","Thermoplasmatota","Verrucomicrobiota","Actinobacteriota","Cyanobacteria","Chloroflexota","Crenarchaeota","Planctomycetota","SAR324","Gemmatimonadota","Acidobacteriota","Nitrospinota","Myxococcota","Binatota","Dadabacteria","Fibrobacterota","Marinisomatota") )

transporters_data.scores_clust$Phylum<-factor(data.scores_clust$Phylum,levels =c( "Proteobacteria","Bacteroidota","Thermoplasmatota","Verrucomicrobiota","Actinobacteriota","Cyanobacteria","Chloroflexota","Crenarchaeota","Planctomycetota","SAR324","Gemmatimonadota","Acidobacteriota","Nitrospinota","Myxococcota","Binatota","Dadabacteria","Fibrobacterota","Marinisomatota") )

functional_transporters_NMDS_Phylum <- ggplot(transporters_data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(color = Phylum)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#C3B4D1", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "NMDS2", shape = "Type") +
  annotate("text", x=-1.3, y=-1.3, label= "0.13",
           col="red", size=10, parse=TRUE) +
  ylim(-1.5, 1.1)

functional_transporters_NMDS_Phylum

##NMDS for gene metabolism colour
metabolism_summary<-dplyr::bind_rows(MISC_combined,carbon_utilization_combined,Transporters_combined,Energy_combined,organic_nitrogen_combined,Woodcroft_combined)

metabolism_summary <- na.omit(metabolism_summary)

com = metabolism_summary[,7:ncol(metabolism_summary)]

m_com <- as.matrix(com)


tot <- rowSums(m_com)
m_com <- m_com[tot > 0, ]

m_com <- na.omit(m_com)

set.seed(123)
nmds = metaMDS(m_com, distance = "bray")
nmds

data.scores = as.data.frame(scores(nmds)$sites)

data.scores <- cbind(metabolism_summary[1:5], data.scores)


###NMDS1 by Bin Size

for_bin_size <- read.csv("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/read_recruitment_to_bins/Yugen_working_files/for_bin_size_bin_metadata.csv")


#nmds_bin_size <- cbind(for_bin_size, data.scores_clust)
#for_bin_size %>%
#join(data.scores_clust, by = c("Bin_ID", "bin_ID"))

nmds_bin_size <- merge(x=for_bin_size, y=data.scores, by.x='Bin_ID', by.y='bin_ID')



abundance_for_bin_size <- subset(abundance, select = -c(Abundance, Sample_ID_2, MOTS_ID,Year,Month,M.Y,sampDate,Season,Published_station,Species, Species_ANI,mOTU,X16S))

data.scores_for_bin_size <- subset(data.scores, select = c(bin_ID, NMDS1, NMDS2))

for_bin_size <- merge(x = data.scores_for_bin_size, y = abundance_for_bin_size, by.x = "bin_ID", by.y = "OTU", all.y = FALSE, no.dups = TRUE)

for_bin_size <- for_bin_size[!duplicated(for_bin_size$bin_ID), ]

for_bin_size$Bin_Size..Mbp. <- as.numeric(for_bin_size$Bin_Size..Mbp.)
for_bin_size$Completeness.... <- as.numeric(for_bin_size$Completeness....)
for_bin_size$Community.... <- as.numeric(for_bin_size$Community....)
for_bin_size$Lon..oE. <- as.numeric(for_bin_size$Lon..oE.)
for_bin_size$Lat..oN. <- as.numeric(for_bin_size$Lat..oN.)
for_bin_size$Distance..km.from.TH. <- as.numeric(for_bin_size$Distance..km.from.TH.)
for_bin_size$Temp..oC. <- as.numeric(for_bin_size$Temp..oC.)
for_bin_size$Salinity..psu. <- as.numeric(for_bin_size$Salinity..psu.)
for_bin_size$Phosphate..mol.m.3. <- as.numeric(for_bin_size$Phosphate..mol.m.3.)
for_bin_size$Nitrate..mol.m.3. <- as.numeric(for_bin_size$Nitrate..mol.m.3.)
for_bin_size$Silicate..mol.m.3. <- as.numeric(for_bin_size$Silicate..mol.m.3.)
for_bin_size$chlA..mg.m3. <- as.numeric(for_bin_size$chlA..mg.m3.)


#NMDS1
ggplot(for_bin_size, aes(x = NMDS1, y = Bin_Size..Mbp.)) + 
  geom_point(size = 4, aes(color = Phylum)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "Bin size (Mbp)", shape = "Type")

ggplot(for_bin_size, aes(x = NMDS1, y = ORF_count)) + 
  geom_point(size = 4, aes(color = Phylum)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "ORF count", shape = "Type")

ggplot(for_bin_size, aes(x = NMDS1, y = Completeness)) + 
  geom_point(size = 4, aes(color = Phylum)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "Completeness (%)", shape = "Type")

ggplot(for_bin_size, aes(x = NMDS1, y = percent_Community)) + 
  geom_point(size = 4, aes(color = Phylum)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "Percentage of community (%)", shape = "Type")

#NMDS2
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

ggplot(for_bin_size, aes(x = NMDS2, y = Completeness....)) + 
  geom_point(size = 4, aes(color = Phylum)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "Completeness (%)", shape = "Type")

ggplot(for_bin_size, aes(x = NMDS2, y = for_bin_size$Community....)) + 
  geom_point(size = 4, aes(color = Phylum)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "Percentage of community (%)", shape = "Type")



ggplot(for_bin_size, aes(x = NMDS2, y = Lon..oE.)) + 
  geom_point(size = 4, aes(color = Phylum)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "Longitude", shape = "Type")

ggplot(for_bin_size, aes(x = NMDS2, y = Lat..oN.)) + 
  geom_point(size = 4, aes(color = Phylum)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "Latitude", shape = "Type")

ggplot(for_bin_size, aes(x = NMDS2, y = Distance..km.from.TH.)) + 
  geom_point(size = 4, aes(color = Phylum)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "Distance from the shore (km)", shape = "Type")


ggplot(for_bin_size, aes(x = NMDS2, y = Temp..oC.)) + 
  geom_point(size = 4, aes(color = Phylum)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "Temperature (celcius)", shape = "Type")


ggplot(for_bin_size, aes(x = NMDS2, y = Salinity..psu.)) + 
  geom_point(size = 4, aes(color = Phylum)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "Salinity (psu)", shape = "Type")


ggplot(for_bin_size, aes(x = NMDS2, y = Phosphate..mol.m.3.)) + 
  geom_point(size = 4, aes(color = Phylum)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "Phosphate (mol/m3)", shape = "Type")



ggplot(for_bin_size, aes(x = NMDS2, y = Nitrate..mol.m.3.)) + 
  geom_point(size = 4, aes(color = Phylum)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "Nitrate (mol/m3)", shape = "Type")

ggplot(for_bin_size, aes(x = NMDS2, y = Silicate..mol.m.3.)) + 
  geom_point(size = 4, aes(color = Phylum)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "Silicate (mol/m3)", shape = "Type")

ggplot(for_bin_size, aes(x = NMDS2, y = chlA..mg.m3.)) + 
  geom_point(size = 4, aes(color = Phylum)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "chlA (mg/mg3)", shape = "Type")
###T-test for bin size and cluster

pairwise.wilcox.test(nmds_bin_size$Bin_Size.Mbp., nmds_bin_size$clust,
                     p.adjust.method = "BH")



#T-test


wilcox.test(Bin_Size.Mbp. ~ clust, data = nmds_bin_size,
            exact = TRUE)

wilcox.test(Completeness ~ cluster, data = nmds_bin_size,
            exact = FALSE)
wilcox.test(ORF_count ~ cluster, data = nmds_bin_size,
            exact = FALSE)


wilcox.test(Bin_Size.Mbp. ~ cluster, data = nmds_bin_size,
            exact = TRUE)

wilcox.test(Bin_Size.Mbp. ~ cluster, data = all_cluster_bin_size_2_3,
            exact = TRUE)





#boxplot
nmds_bin_size$cluster <- str_replace(nmds_bin_size$cluster, "clust1", "Cluster 1")
nmds_bin_size$cluster <- str_replace(nmds_bin_size$cluster, "clust2", "Cluster 2")
nmds_bin_size$cluster <- str_replace(nmds_bin_size$cluster, "clust3", "Cluster 3")

library('rstatix')
Bin_Size_clust<- ggplot(nmds_bin_size, aes(x = cluster, y = Bin_Size..Mbp.)) +
  geom_boxplot() +
  labs(y = "Bin size (Mbp)", x = "Cluster", tag = "C") +
  stat_pvalue_manual(nmds_bin_size %>% 
                       wilcox_test(Bin_Size..Mbp. ~ cluster, exact = FALSE) %>% 
                       add_xy_position()) +
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
        panel.background = element_blank())

Bin_Size_clust



ORF_count_clust <- ggplot(nmds_bin_size, aes(x = cluster, y = ORF_count)) +
  geom_boxplot() +
  labs(y = "Number of ORFs", x = "Cluster", tag = "D") +
  stat_pvalue_manual(nmds_bin_size %>% 
                       wilcox_test(ORF_count ~ cluster, exact = NULL) %>% 
                       add_xy_position()) +
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
        panel.background = element_blank())

ORF_count_clust

percent_community_clust <- ggplot(nmds_bin_size, aes(x = clust, y = percent_Community)) +
  geom_boxplot() +
  labs(y = "Percent of community") +
  stat_pvalue_manual(nmds_bin_size %>% 
                       wilcox_test(percent_Community ~ clust, exact = NULL) %>% 
                       add_xy_position())+
  My_Theme

percent_community_clust

ggplot(nmds_bin_size, aes(x = cluster, y = Completeness)) +
  geom_boxplot() +
  stat_pvalue_manual(nmds_bin_size %>% 
                       wilcox_test(Completeness ~ cluster, exact = NULL) %>% 
                       add_xy_position())+
  My_Theme

ggplot(nmds_bin_size, aes(x = clust, y = Contamination)) +
  geom_boxplot() +
  stat_pvalue_manual(nmds_bin_size %>% 
                       wilcox_test(Contamination ~ clust, exact = NULL) %>% 
                       add_xy_position())+
  My_Theme


###Bin size ORF regression
library(ggpmisc)

bin_size_ORF_regression <- ggplot(nmds_bin_size, aes(x = ORF_count, y = Bin_Size.Mbp.)) +
  geom_point() +
  labs(x = "ORF count", y = "Bin size (Mbp)", tag = "B") +
  stat_poly_line() +
  stat_poly_eq() +
  My_Theme

bin_size_ORF_regression

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/bin_sie_ORF_regression.pdf",width=12,height=10)

bin_size_ORF_regression

dev.off()

##Boxplot 
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

##ORF_binsize NMDS
names(nmds_bin_size)[names(nmds_bin_size) == 'Bin_ID'] <- 'bin_ID'
data.scores_ORF_binsize <- merge(data.scores, nmds_bin_size, by = "bin_ID")

NMDS1_ORF_binsize <- ggplot(data.scores_ORF_binsize, aes(x = NMDS1, y = ORF_binsize)) + 
  geom_point(size = 4, aes(color = Phylum.y)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  theme(axis.title.x = element_text(face="bold",size=24),
        axis.text.x = element_text(colour = "black", size=22), 
        axis.text.y = element_text(colour = "black", size=22),
        axis.title.y = element_text(face="bold",size=24),
        plot.title = element_text(size = 24),
        strip.text.x = element_text(size=22, face="bold"),
        strip.text.y = element_text(size=22, face="bold"),
        panel.border = element_rect(fill = NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none") + 
  labs(x = "NMDS1", colour = "Phylum.y", y = "Gene density (ORF/Bin size)", shape = "Type", tag = "A")

NMDS2_ORF_binsize <- ggplot(data.scores_ORF_binsize, aes(x = NMDS2, y = ORF_binsize)) + 
  geom_point(size = 4, aes(color = Phylum.y)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  theme(axis.title.x = element_text(face="bold",size=24),
        axis.text.x = element_text(colour = "black", size=22), 
        axis.text.y = element_text(colour = "black", size=22),
        axis.title.y = element_text(face="bold",size=24),
        plot.title = element_text(size = 24),
        strip.text.x = element_text(size=22, face="bold"),
        strip.text.y = element_text(size=22, face="bold"),
        panel.border = element_rect(fill = NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none") +
  labs(x = "NMDS2", colour = "Phylum.y", y = "Gene density (ORF/Bin size)", shape = "Type", tag = "B")

NMDS2_ORF_binsize
#ANOVA
# Compute the analysis of variance
res.aov <- aov(Bin_Size.Mbp. ~ cluster, data = nmds_bin_size)
# Summary of the analysis
summary(res.aov)



pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/Figure_4.pdf",15, height = 13)

grid.arrange(Cluster1_2_Cren_Therm_gene_count_2_panel, metabolism_summary_a.a_n_p,
             widths = c(2, 2))

#grid.arrange(Cluster1_Cren_Therm_gene_plot, Cluster3_gene_plot, metabolism_summary_a.a_n_p,
#layout_matrix = rbind(c(1, 2),
# c(3)))

dev.off()

###to create correlation by header

data.scores_genes <- data.scores
library(data.table)
ADO_CBL_synthesis <- metabolism_summary[metabolism_summary$header %like% "ADO-CBL synthesis", ]


ADO_CBL_synthesis <- colSums(ADO_CBL_synthesis [ , c(7:1014)], na.rm=TRUE)

ADO_CBL_synthesis <- data.frame(ADO_CBL_synthesis)

library(tibble)
ADO_CBL_synthesis <- tibble::rownames_to_column(ADO_CBL_synthesis, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(ADO_CBL_synthesis, by = c("bin_ID"))

ADO_CBL_synthesis_NMDS <- ggplot(data.scores_genes, aes(x = NMDS1, y = ADO_CBL_synthesis)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "ADO-CBL synthesis", shape = "Type") +
  stat_cor(method="spearman", size = 7)

#+ 
#stat_smooth(method = lm, formula = y ~ poly(x, 3, raw = TRUE))

ADO_CBL_synthesis_NMDS

lm(NMDS1~ADO_CBL_synthesis, I(ADO_CBL_synthesis^2), data = data.scores_genes)
lm(NMDS1~poly(ADO_CBL_synthesis, 1, raw = TRUE), data = data.scores_genes)

corr <- cor.test(x=data.scores_genes$NMDS1, y=data.scores_genes$ADO_CBL_synthesis, method = 'spearman')
corr

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/NMDS/ADO_CBL_synthesis_NMDS.pdf",width=12,height=10) # Open a new pdf file

#view
ADO_CBL_synthesis_NMDS

dev.off()
#Save as PDF


Antibiotic_Resistance <- metabolism_summary[metabolism_summary$header %like% "Antibiotic Resistance", ]
Antibiotic_Resistance <- colSums(Antibiotic_Resistance [ , c(7:1014)], na.rm=TRUE)

Antibiotic_Resistance <- data.frame(Antibiotic_Resistance)

library(tibble)
Antibiotic_Resistance <- tibble::rownames_to_column(Antibiotic_Resistance, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(Antibiotic_Resistance, by = c("bin_ID"))

AMR_NMDS <- ggplot(data.scores_genes, aes(x = NMDS1, y = Antibiotic_Resistance)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "Antibiotic Resistance", shape = "Type") +
  stat_cor(method="spearman", size = 7)

AMR_NMDS

corr <- cor.test(x=data.scores_genes$NMDS1, y=data.scores_genes$Antibiotic_Resistance, method = 'spearman')
corr

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/NMDS/AMR_NMDS.pdf",width=12,height=10) # Open a new pdf file

#view
AMR_NMDS

dev.off()
#Save as PDF


CRISPR <- metabolism_summary[metabolism_summary$header %like% "CRISPR", ]
CRISPR <- colSums(CRISPR [ , c(7:1014)], na.rm=TRUE)

CRISPR <- data.frame(CRISPR)

library(tibble)
CRISPR <- tibble::rownames_to_column(CRISPR, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(CRISPR, by = c("bin_ID"))

ggplot(data.scores_genes, aes(x = NMDS1, y = CRISPR)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "CRISPR", shape = "Type")

corr <- cor.test(x=data.scores_genes$NMDS1, y=data.scores_genes$CRISPR, method = 'spearman')
corr


Flagella_Structure <- metabolism_summary[metabolism_summary$header %like% "Flagella Structure", ]
Flagella_Structure <- colSums(Flagella_Structure [ , c(7:1014)], na.rm=TRUE)

Flagella_Structure <- data.frame(Flagella_Structure)

library(tibble)
Flagella_Structure <- tibble::rownames_to_column(Flagella_Structure, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(Flagella_Structure, by = c("bin_ID"))

Flagella_Structure_NMDS <- ggplot(data.scores_genes, aes(x = NMDS1, y = Flagella_Structure)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "Flagella structure gene counts per bin", shape = "Type") +
  stat_cor(method="spearman", size = 7)

corr <- cor.test(x=data.scores_genes$NMDS1, y=data.scores_genes$Flagella_Structure, method = 'spearman')
corr

Flagella_Structure_NMDS

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/NMDS/Flagella_Structure_NMDS.pdf",width=12,height=10) # Open a new pdf file

#view
Flagella_Structure_NMDS

dev.off()
#Save as PDF

#Flagellar_ctyoplasmic_structure <- metabolism_summary[metabolism_summary$header %like% "Flagellar cytoplasmic structure", ]
#Flagellar_ctyoplasmic_structure <- colSums(Flagellar_ctyoplasmic_structure [ , c(7:1014)], na.rm=TRUE)

#data.scores_genes$Flagellar_ctyoplasmic_structure = c(Flagellar_ctyoplasmic_structure)

#ggplot(data.scores_genes, aes(x = NMDS1, y = Flagellar_ctyoplasmic_structure)) + 
geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "Flagellar ctyoplasmic structure", shape = "Type")

#corr <- cor.test(x=data.scores_genes$NMDS1, y=data.scores_genes$Flagellar_ctyoplasmic_structure, method = 'spearman')
#corr

Information_systems <- metabolism_summary[metabolism_summary$header %like% "Information systems", ]
Information_systems <- colSums(Information_systems [ , c(7:1014)], na.rm=TRUE)

Information_systems <- data.frame(Information_systems)

library(tibble)
Information_systems <- tibble::rownames_to_column(Information_systems, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(Information_systems, by = c("bin_ID"))

ggplot(data.scores_genes, aes(x = NMDS1, y = Information_systems)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#B4CCDF", "Cyanobacteria" = "#4A75AA", "Dadabacteria" = "#C2DD9B", "Fibrobacterota" = "#669C4C", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#C3B4D1", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "Information systems", shape = "Type")

corr <- cor.test(x=data.scores_genes$NMDS1, y=data.scores_genes$Information_systems, method = 'spearman')
corr

##Trying to see what is in IS


Information_systems <- metabolism_summary[metabolism_summary$header %like% "Information systems", ]

I_S <- subset(Information_systems, select = -c(gene_id, gene_description,function_description,subheader))



I_S <- gather(I_S, bin, count, '10_bin.1':'9_bin.9', factor_key=TRUE)



I_S_cluster <- merge(I_S, data.scores, by.x = 'bin', by.y = 'bin_ID')

I_S_cluster <- subset(I_S_cluster, select = -c(bin, Sample, Watermass.y,Season,NMDS1, NMDS2))

library(data.table)
setDT(I_S_cluster)
I_S_cluster1  = I_S_cluster [ , .(counted = sum(count)), by = .(module, header, Phylum, cluster)]

ggplot(I_S_cluster1, aes(x = reorder(module, -counted), y = log(counted))) +
  geom_boxplot() +
  facet_wrap(~cluster, ncol = 1)+
  theme(axis.text.x = element_text(angle=90))

ggboxplot(I_S_cluster1 , x = "module", y = "counted", 
          add = "jitter",
          facet.by = "cluster", width = .5) +
  rotate_x_text(90)

MISC <- metabolism_summary[metabolism_summary$header %like% "MISC", ]
MISC <- colSums(MISC [ , c(7:1014)], na.rm=TRUE)

MISC <- data.frame(MISC)

library(tibble)
MISC <- tibble::rownames_to_column(MISC, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(MISC, by = c("bin_ID"))

MISC_NMDS <- ggplot(data.scores_genes, aes(x = NMDS1, y = MISC)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "MISC", shape = "Type")+
  stat_cor(method="spearman", size = 7)

corr <- cor.test(x=data.scores_genes$NMDS1, y=data.scores_genes$MISC, method = 'spearman')
corr

MISC_NMDS

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/NMDS/MISC_NMDS.pdf",width=12,height=10) # Open a new pdf file

#view
MISC_NMDS

dev.off()
#Save as PDF

#SCFA_and_alcohol_conversions <- metabolism_summary[metabolism_summary$header %like% "SCFA and alcohol conversions", ]
#SCFA_and_alcohol_conversions <- colSums(SCFA_and_alcohol_conversions [ , c(7:1014)], na.rm=TRUE)

#data.scores_genes$SCFA_and_alcohol_conversions = c(SCFA_and_alcohol_conversions)

#ggplot(data.scores_genes, aes(x = NMDS1, y = SCFA_and_alcohol_conversions)) + 
geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "SCFA and alcohol conversions", shape = "Type")

#corr <- cor.test(x=data.scores_genes$NMDS1, y=data.scores_genes$SCFA_and_alcohol_conversions, method = 'spearman')
#corr

aerobic_corrin_ring_synthesis <- metabolism_summary[metabolism_summary$header %like% "aerobic corrin ring synthesis", ]
aerobic_corrin_ring_synthesis <- colSums(aerobic_corrin_ring_synthesis [ , c(7:1014)], na.rm=TRUE)

aerobic_corrin_ring_synthesis <- data.frame(aerobic_corrin_ring_synthesis)

library(tibble)
aerobic_corrin_ring_synthesis <- tibble::rownames_to_column(aerobic_corrin_ring_synthesis, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(aerobic_corrin_ring_synthesis, by = c("bin_ID"))

aerobic_corrin_ring_synthesis_NMDS <- ggplot(data.scores_genes, aes(x = NMDS1, y = aerobic_corrin_ring_synthesis)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "aerobic corrin ring synthesis", shape = "Type")+
  stat_cor(method="spearman", size = 7)

corr <- cor.test(x=data.scores_genes$NMDS1, y=data.scores_genes$aerobic_corrin_ring_synthesis, method = 'spearman')
corr

aerobic_corrin_ring_synthesis_NMDS

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/NMDS/aerobic_corrin_ring_synthesis_NMDS.pdf",width=12,height=10) # Open a new pdf file

#view
aerobic_corrin_ring_synthesis_NMDS

dev.off()
#Save as PDF

#anerobic_corrin_ring_synthesis <- metabolism_summary[metabolism_summary$header %like% "anerobic corrin ring synthesis", ]
#anerobic_corrin_ring_synthesis <- colSums(anerobic_corrin_ring_synthesis [ , c(7:1014)], na.rm=TRUE)

#data.scores_genes$anerobic_corrin_ring_synthesis = c(anerobic_corrin_ring_synthesis)

#ggplot(data.scores_genes, aes(x = NMDS1, y = anerobic_corrin_ring_synthesis)) + 
geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "anerobic corrin ring synthesis", shape = "Type")

#corr <- cor.test(x=data.scores_genes$NMDS1, y=data.scores_genes$anerobic_corrin_ring_synthesis, method = 'spearman')
#corr

hydrocarbon_degradation <- metabolism_summary[metabolism_summary$header %like% "hydrocarbon degradation", ]
hydrocarbon_degradation <- colSums(hydrocarbon_degradation [ , c(7:1014)], na.rm=TRUE)

hydrocarbon_degradation <- data.frame(hydrocarbon_degradation)

library(tibble)
hydrocarbon_degradation <- tibble::rownames_to_column(hydrocarbon_degradation, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(hydrocarbon_degradation, by = c("bin_ID"))

hydrocarbon_degradation_NMDS <- ggplot(data.scores_genes, aes(x = NMDS1, y = hydrocarbon_degradation)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "hydrocarbon degradation", shape = "Type")+
  stat_cor(method="spearman", size = 7)

corr <- cor.test(x=data.scores_genes$NMDS1, y=data.scores_genes$hydrocarbon_degradation, method = 'spearman')
corr

hydrocarbon_degradation_NMDS

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/NMDS/hydrocarbon_degradation_NMDS.pdf",width=12,height=10) # Open a new pdf file

#view
hydrocarbon_degradation_NMDS

dev.off()
#Save as PDF

CAZY <- metabolism_summary[metabolism_summary$header %like% "CAZY", ]
CAZY <- colSums(CAZY [ , c(7:1014)], na.rm=TRUE)

CAZY <- data.frame(CAZY)

library(tibble)
CAZY <- tibble::rownames_to_column(CAZY, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(CAZY, by = c("bin_ID"))

CAZY_NMDS <- ggplot(data.scores_genes, aes(x = NMDS1, y = CAZY)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "CAZY)", shape = "Type")+
  stat_cor(method="spearman", size = 7)

corr <- cor.test(x=data.scores_genes$NMDS1, y=data.scores_genes$CAZY, method = 'spearman')
corr

CAZY_NMDS

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/NMDS/CAZY_NMDS.pdf",width=12,height=10) # Open a new pdf file

#view
CAZY_NMDS

dev.off()
#Save as PDF

central_carbon <- metabolism_summary[metabolism_summary$header %like% "central carbon", ]
central_carbon <- colSums(central_carbon [ , c(7:1014)], na.rm=TRUE)

central_carbon <- data.frame(central_carbon)

library(tibble)
central_carbon <- tibble::rownames_to_column(central_carbon, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(central_carbon, by = c("bin_ID"))

central_carbon_NMDS <- ggplot(data.scores_genes, aes(x = NMDS1, y = central_carbon)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "central_carbon", shape = "Type")+
  stat_cor(method="spearman", size = 7)

corr <- cor.test(x=data.scores_genes$NMDS1, y=data.scores_genes$central_carbon, method = 'spearman')
corr

central_carbon_NMDS

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/NMDS/central_carbon_NMDS.pdf",width=12,height=10) # Open a new pdf file

#view
central_carbon_NMDS

dev.off()
#Save as PDF

pyruvate_metabolism <- metabolism_summary[metabolism_summary$header %like% "pyruvate metabolism", ]
pyruvate_metabolism <- colSums(pyruvate_metabolism [ , c(7:1014)], na.rm=TRUE)

pyruvate_metabolism <- data.frame(pyruvate_metabolism)

library(tibble)
pyruvate_metabolism <- tibble::rownames_to_column(pyruvate_metabolism, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(pyruvate_metabolism, by = c("bin_ID"))

ggplot(data.scores_genes, aes(x = NMDS1, y = pyruvate_metabolism)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "pyruvate metabolism", shape = "Type")+
  stat_cor(method="spearman", size = 7)

corr <- cor.test(x=data.scores_genes$NMDS1, y=data.scores_genes$pyruvate_metabolism, method = 'spearman')
corr

#sugar_utilization_woodcroft <- metabolism_summary[metabolism_summary$header %like% "sugar utilization (woodcroft)", ]
#sugar_utilization_woodcroft <- colSums(sugar_utilization_woodcroft [ , c(7:1014)], na.rm=TRUE)

#data.scores_genes$sugar_utilization_woodcroft = c(sugar_utilization_woodcroft)

#ggplot(data.scores_genes, aes(x = NMDS1, y = sugar_utilization_woodcroft) + 
geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "sugar utilization (woodcroft)", shape = "Type")

#corr <- cor.test(x=data.scores_genes$NMDS1, y=data.scores_genes$sugar_utilization_woodcroft, method = 'spearman')
#corr

Vitamin_B12_transport_system <- metabolism_summary[metabolism_summary$header %like% "Vitamin B12 transport system", ]
Vitamin_B12_transport_system <- colSums(Vitamin_B12_transport_system [ , c(7:1014)], na.rm=TRUE)

Vitamin_B12_transport_system <- data.frame(Vitamin_B12_transport_system)

library(tibble)
Vitamin_B12_transport_system <- tibble::rownames_to_column(Vitamin_B12_transport_system, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(Vitamin_B12_transport_system, by = c("bin_ID"))

ggplot(data.scores_genes, aes(x = NMDS1, y = Vitamin_B12_transport_system)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "Vitamin B12 transport system", shape = "Type")+
  stat_cor(method="spearman", size = 7)

corr <- cor.test(x=data.scores_genes$NMDS1, y=data.scores_genes$Vitamin_B12_transport_system, method = 'spearman')
corr


C1 <- metabolism_summary[metabolism_summary$header %like% "C1", ]
C1 <- colSums(C1 [ , c(7:1014)], na.rm=TRUE)

C1 <- data.frame(C1)

library(tibble)
C1 <- tibble::rownames_to_column(C1, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(C1, by = c("bin_ID"))

ggplot(data.scores_genes, aes(x = NMDS1, y = C1)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "C1", shape = "Type")+
  stat_cor(method="spearman", size = 7)

corr <- cor.test(x=data.scores_genes$NMDS1, y=data.scores_genes$C1, method = 'spearman')
corr


#C1_methane <- metabolism_summary[metabolism_summary$header %like% "C1_methane", ]
#C1_methane <- colSums(C1_methane [ , c(7:1014)], na.rm=TRUE)

#data.scores_genes$C1_methane = c(C1_methane)

#ggplot(data.scores_genes, aes(x = NMDS1, y = C1_methane)) + 
geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "C1-methane", shape = "Type")

#corr <- cor.test(x=data.scores_genes$NMDS1, y=data.scores_genes$C1_methane, method = 'spearman')
#corr


Electron_transport_Chain <- metabolism_summary[metabolism_summary$header %like% "Electron transport Chain", ]
Electron_transport_Chain <- colSums(Electron_transport_Chain [ , c(7:1014)], na.rm=TRUE)

Electron_transport_Chain <- data.frame(Electron_transport_Chain)

library(tibble)
Electron_transport_Chain <- tibble::rownames_to_column(Electron_transport_Chain, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(Electron_transport_Chain, by = c("bin_ID"))

Electron_transport_Chain_NMDS <- ggplot(data.scores_genes, aes(x = NMDS1, y = Electron_transport_Chain)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "Electron transport Chain", shape = "Type")+
  stat_cor(method="spearman", size = 7)

corr <- cor.test(x=data.scores_genes$NMDS1, y=data.scores_genes$Electron_transport_Chain, method = 'spearman')
corr

Electron_transport_Chain_NMDS

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/NMDS/Electron_transport_Chain_NMDS.pdf",width=12,height=10) # Open a new pdf file

#view
Electron_transport_Chain_NMDS

dev.off()
#Save as PDF


Hydrogenases <- metabolism_summary[metabolism_summary$header %like% "Hydrogenases", ]
Hydrogenases <- colSums(Hydrogenases [ , c(7:1014)], na.rm=TRUE)

Hydrogenases <- data.frame(Hydrogenases)

library(tibble)
Hydrogenases <- tibble::rownames_to_column(Hydrogenases, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(Hydrogenases, by = c("bin_ID"))

ggplot(data.scores_genes, aes(x = NMDS1, y = Hydrogenases)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "Hydrogenases", shape = "Type")

corr <- cor.test(x=data.scores_genes$NMDS1, y=data.scores_genes$Hydrogenases, method = 'spearman')
corr


Metal_Reduction <- metabolism_summary[metabolism_summary$header %like% "Metal Reduction", ]
Metal_Reduction <- colSums(Metal_Reduction [ , c(7:1014)], na.rm=TRUE)

Metal_Reduction <- data.frame(Metal_Reduction)

library(tibble)
Metal_Reduction <- tibble::rownames_to_column(Metal_Reduction, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(Metal_Reduction, by = c("bin_ID"))

ggplot(data.scores_genes, aes(x = NMDS1, y = Metal_Reduction)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "Metal Reduction", shape = "Type")

corr <- cor.test(x=data.scores_genes$NMDS1, y=data.scores_genes$Metal_Reduction, method = 'spearman')
corr

Nitrogen <- metabolism_summary[metabolism_summary$header %like% "Nitrogen", ]
Nitrogen <- colSums(Nitrogen [ , c(7:1014)], na.rm=TRUE)

Nitrogen <- data.frame(Nitrogen)

library(tibble)
Nitrogen <- tibble::rownames_to_column(Nitrogen, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(Nitrogen, by = c("bin_ID"))

ggplot(data.scores_genes, aes(x = NMDS1, y = Nitrogen)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "Nitrogen", shape = "Type")+
  stat_cor(method="spearman", size = 7)

corr <- cor.test(x=data.scores_genes$NMDS1, y=data.scores_genes$Nitrogen, method = 'spearman')
corr

Oxygen <- metabolism_summary[metabolism_summary$header %like% "Oxygen", ]
Oxygen <- colSums(Oxygen [ , c(7:1014)], na.rm=TRUE)

Oxygen <- data.frame(Oxygen)

library(tibble)
Oxygen <- tibble::rownames_to_column(Oxygen, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(Oxygen, by = c("bin_ID"))

Oxygen_NMDS <- ggplot(data.scores_genes, aes(x = NMDS1, y = Oxygen)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "Oxygen", shape = "Type")+
  stat_cor(method="spearman", size = 7)

corr <- cor.test(x=data.scores_genes$NMDS1, y=data.scores_genes$Oxygen, method = 'spearman')
corr

Oxygen_NMDS

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/NMDS/Oxygen_NMDS.pdf",width=12,height=10) # Open a new pdf file

#view
Oxygen_NMDS

dev.off()
#Save as PDF


Photosynthesis <- metabolism_summary[metabolism_summary$header %like% "Photosynthesis", ]
Photosynthesis <- colSums(Photosynthesis [ , c(7:1014)], na.rm=TRUE)

Photosynthesis <- data.frame(Photosynthesis)

library(tibble)
Photosynthesis <- tibble::rownames_to_column(Photosynthesis, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(Photosynthesis, by = c("bin_ID"))

Photosynthesis_NMDS <- ggplot(data.scores_genes, aes(x = NMDS1, y = Photosynthesis)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "Photosynthesis", shape = "Type")+
  stat_cor(method="spearman", size = 7)

corr <- cor.test(x=data.scores_genes$NMDS1, y=data.scores_genes$Photosynthesis, method = 'spearman')
corr

Photosynthesis_NMDS

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/NMDS/Photosynthesis_NMDS.pdf",width=12,height=10) # Open a new pdf file

#view
Photosynthesis_NMDS

dev.off()
#Save as PDF


Sulfur <- metabolism_summary[metabolism_summary$header %like% "Sulfur", ]
Sulfur <- colSums(Sulfur [ , c(7:1014)], na.rm=TRUE)

Sulfur <- data.frame(Sulfur)

library(tibble)
Sulfur <- tibble::rownames_to_column(Sulfur, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(Sulfur, by = c("bin_ID"))

Sulfur_NMDS <- ggplot(data.scores_genes, aes(x = NMDS1, y = Sulfur)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "Sulfur", shape = "Type")+
  stat_cor(method="spearman", size = 7)

corr <- cor.test(x=data.scores_genes$NMDS1, y=data.scores_genes$Sulfur, method = 'spearman')
corr

Sulfur_NMDS

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/NMDS/Sulfur_NMDS.pdf",width=12,height=10) # Open a new pdf file

#view
Sulfur_NMDS

dev.off()
#Save as PDF

Amino_Acid <- metabolism_summary[metabolism_summary$header %like% "Amino Acid", ]
Amino_Acid <- colSums(Amino_Acid [ , c(7:1014)], na.rm=TRUE)

Amino_Acid <- data.frame(Amino_Acid)

library(tibble)
Amino_Acid <- tibble::rownames_to_column(Amino_Acid, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(Amino_Acid, by = c("bin_ID"))

Amino_Acid_NMDS <- ggplot(data.scores_genes, aes(x = NMDS1, y = Amino_Acid)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "Amino acid gene counts per bin", shape = "Type")+
  stat_cor(method="spearman", size = 7)

corr <- cor.test(x=data.scores_genes$NMDS1, y=data.scores_genes$Amino_Acid, method = 'spearman')
corr

Amino_Acid_NMDS

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/NMDS/Amino_Acid_NMDS.pdf",width=12,height=10) # Open a new pdf file

#view
Amino_Acid_NMDS

dev.off()
#Save as PDF

Peptidase <- metabolism_summary[metabolism_summary$header %like% "Peptidase", ]
Peptidase <- colSums(Peptidase [ , c(7:1014)], na.rm=TRUE)

Peptidase <- data.frame(Peptidase)

library(tibble)
Peptidase <- tibble::rownames_to_column(Peptidase, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(Peptidase, by = c("bin_ID"))

peptidase_NMDS <- ggplot(data.scores_genes, aes(x = NMDS1, y = Peptidase)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "Peptidase", shape = "Type")+
  stat_cor(method="spearman", size = 7)

corr <- cor.test(x=data.scores_genes$NMDS1, y=data.scores_genes$Peptidase, method = 'spearman')
corr

peptidase_NMDS

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/NMDS/peptidase_NMDS.pdf",width=12,height=10) # Open a new pdf file

#view
peptidase_NMDS

dev.off()
#Save as PDF

###FOR NMDS2 AND GENES
ADO_CBL_synthesis_NMDS2 <- ggplot(data.scores_genes, aes(x = NMDS2, y = ADO_CBL_synthesis)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "ADO-CBL synthesis", shape = "Type") +
  stat_cor(method="spearman", size = 7)

#+ 
#stat_smooth(method = lm, formula = y ~ poly(x, 3, raw = TRUE))

ADO_CBL_synthesis_NMDS2

lm(NMDS1~ADO_CBL_synthesis, I(ADO_CBL_synthesis^2), data = data.scores_genes)
lm(NMDS1~poly(ADO_CBL_synthesis, 1, raw = TRUE), data = data.scores_genes)

corr <- cor.test(x=data.scores_genes$NMDS2, y=data.scores_genes$ADO_CBL_synthesis, method = 'spearman')
corr

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/NMDS/ADO_CBL_synthesis_NMDS.pdf",width=12,height=10) # Open a new pdf file

#view
ADO_CBL_synthesis_NMDS2

dev.off()
#Save as PDF


Antibiotic_Resistance <- metabolism_summary[metabolism_summary$header %like% "Antibiotic Resistance", ]
Antibiotic_Resistance <- colSums(Antibiotic_Resistance [ , c(7:1014)], na.rm=TRUE)

Antibiotic_Resistance <- data.frame(Antibiotic_Resistance)

library(tibble)
Antibiotic_Resistance <- tibble::rownames_to_column(Antibiotic_Resistance, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(Antibiotic_Resistance, by = c("bin_ID"))

AMR_NMDS2 <- ggplot(data.scores_genes, aes(x = NMDS2, y = Antibiotic_Resistance)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "Antibiotic Resistance", shape = "Type") +
  stat_cor(method="spearman", size = 7)

AMR_NMDS2

corr <- cor.test(x=data.scores_genes$NMDS2, y=data.scores_genes$Antibiotic_Resistance, method = 'spearman')
corr

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/NMDS/AMR_NMDS.pdf",width=12,height=10) # Open a new pdf file

#view
AMR_NMDS2

dev.off()
#Save as PDF


CRISPR <- metabolism_summary[metabolism_summary$header %like% "CRISPR", ]
CRISPR <- colSums(CRISPR [ , c(7:1014)], na.rm=TRUE)

CRISPR <- data.frame(CRISPR)

library(tibble)
CRISPR <- tibble::rownames_to_column(CRISPR, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(CRISPR, by = c("bin_ID"))

ggplot(data.scores_genes, aes(x = NMDS2, y = CRISPR)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "CRISPR", shape = "Type")

corr <- cor.test(x=data.scores_genes$NMDS2, y=data.scores_genes$CRISPR, method = 'spearman')
corr


Flagella_Structure <- metabolism_summary[metabolism_summary$header %like% "Flagella Structure", ]
Flagella_Structure <- colSums(Flagella_Structure [ , c(7:1014)], na.rm=TRUE)

Flagella_Structure <- data.frame(Flagella_Structure)

library(tibble)
Flagella_Structure <- tibble::rownames_to_column(Flagella_Structure, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(Flagella_Structure, by = c("bin_ID"))

Flagella_Structure_NMDS2 <- ggplot(data.scores_genes, aes(x = NMDS2, y = Flagella_Structure)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "Flagella structure gene counts per bin", shape = "Type") +
  stat_cor(method="spearman", size = 7)

corr <- cor.test(x=data.scores_genes$NMDS2, y=data.scores_genes$Flagella_Structure, method = 'spearman')
corr

Flagella_Structure_NMDS2

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/NMDS/Flagella_Structure_NMDS.pdf",width=12,height=10) # Open a new pdf file

#view
Flagella_Structure_NMDS2

dev.off()
#Save as PDF

#Flagellar_ctyoplasmic_structure <- metabolism_summary[metabolism_summary$header %like% "Flagellar cytoplasmic structure", ]
#Flagellar_ctyoplasmic_structure <- colSums(Flagellar_ctyoplasmic_structure [ , c(7:1014)], na.rm=TRUE)

#data.scores_genes$Flagellar_ctyoplasmic_structure = c(Flagellar_ctyoplasmic_structure)

ggplot(data.scores_genes, aes(x = NMDS2, y = Flagellar_ctyoplasmic_structure)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "Flagellar ctyoplasmic structure", shape = "Type")

#corr <- cor.test(x=data.scores_genes$NMDS1, y=data.scores_genes$Flagellar_ctyoplasmic_structure, method = 'spearman')
#corr

Information_systems <- metabolism_summary[metabolism_summary$header %like% "Information systems", ]
Information_systems <- colSums(Information_systems [ , c(7:1014)], na.rm=TRUE)

Information_systems <- data.frame(Information_systems)

library(tibble)
Information_systems <- tibble::rownames_to_column(Information_systems, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(Information_systems, by = c("bin_ID"))

data.scores_genes <- merge(data.scores_genes, Information_systems, by = "bin_ID")

Information_systems_NMDS2 <- ggplot(data.scores_genes, aes(x = NMDS2, y = Information_systems)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#B4CCDF", "Cyanobacteria" = "#4A75AA", "Dadabacteria" = "#C2DD9B", "Fibrobacterota" = "#669C4C", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#C3B4D1", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "Information systems", shape = "Type") +
  stat_cor(method="spearman", size = 7)

corr <- cor.test(x=data.scores_genes$NMDS2, y=data.scores_genes$Information_systems, method = 'spearman')
corr

MISC <- metabolism_summary[metabolism_summary$header %like% "MISC", ]
MISC <- colSums(MISC [ , c(7:1014)], na.rm=TRUE)

MISC <- data.frame(MISC)

library(tibble)
MISC <- tibble::rownames_to_column(MISC, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(MISC, by = c("bin_ID"))

MISC_NMDS2 <- ggplot(data.scores_genes, aes(x = NMDS2, y = MISC)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "MISC", shape = "Type")+
  stat_cor(method="spearman", size = 7)

corr <- cor.test(x=data.scores_genes$NMDS2, y=data.scores_genes$MISC, method = 'spearman')
corr

MISC_NMDS2

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/NMDS/MISC_NMDS.pdf",width=12,height=10) # Open a new pdf file

#view
MISC_NMDS2

dev.off()
#Save as PDF

#SCFA_and_alcohol_conversions <- metabolism_summary[metabolism_summary$header %like% "SCFA and alcohol conversions", ]
#SCFA_and_alcohol_conversions <- colSums(SCFA_and_alcohol_conversions [ , c(7:1014)], na.rm=TRUE)

#data.scores_genes$SCFA_and_alcohol_conversions = c(SCFA_and_alcohol_conversions)

ggplot(data.scores_genes, aes(x = NMDS2, y = SCFA_and_alcohol_conversions)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "SCFA and alcohol conversions", shape = "Type")

#corr <- cor.test(x=data.scores_genes$NMDS1, y=data.scores_genes$SCFA_and_alcohol_conversions, method = 'spearman')
#corr

aerobic_corrin_ring_synthesis <- metabolism_summary[metabolism_summary$header %like% "aerobic corrin ring synthesis", ]
aerobic_corrin_ring_synthesis <- colSums(aerobic_corrin_ring_synthesis [ , c(7:1014)], na.rm=TRUE)

aerobic_corrin_ring_synthesis <- data.frame(aerobic_corrin_ring_synthesis)

library(tibble)
aerobic_corrin_ring_synthesis <- tibble::rownames_to_column(aerobic_corrin_ring_synthesis, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(aerobic_corrin_ring_synthesis, by = c("bin_ID"))

aerobic_corrin_ring_synthesis_NMDS2 <- ggplot(data.scores_genes, aes(x = NMDS2, y = aerobic_corrin_ring_synthesis)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "aerobic corrin ring synthesis", shape = "Type")+
  stat_cor(method="spearman", size = 7)

corr <- cor.test(x=data.scores_genes$NMDS1, y=data.scores_genes$aerobic_corrin_ring_synthesis, method = 'spearman')
corr

aerobic_corrin_ring_synthesis_NMDS2

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/NMDS/aerobic_corrin_ring_synthesis_NMDS.pdf",width=12,height=10) # Open a new pdf file

#view
aerobic_corrin_ring_synthesis_NMDS2

dev.off()
#Save as PDF

anerobic_corrin_ring_synthesis <- metabolism_summary[metabolism_summary$header %like% "anerobic corrin ring synthesis", ]
anerobic_corrin_ring_synthesis <- colSums(anerobic_corrin_ring_synthesis [ , c(7:1014)], na.rm=TRUE)

data.scores_genes$anerobic_corrin_ring_synthesis = c(anerobic_corrin_ring_synthesis)

ggplot(data.scores_genes, aes(x = NMDS2, y = anerobic_corrin_ring_synthesis)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "anerobic corrin ring synthesis", shape = "Type")

#corr <- cor.test(x=data.scores_genes$NMDS1, y=data.scores_genes$anerobic_corrin_ring_synthesis, method = 'spearman')
#corr

hydrocarbon_degradation <- metabolism_summary[metabolism_summary$header %like% "hydrocarbon degradation", ]
hydrocarbon_degradation <- colSums(hydrocarbon_degradation [ , c(7:1014)], na.rm=TRUE)

hydrocarbon_degradation <- data.frame(hydrocarbon_degradation)

library(tibble)
hydrocarbon_degradation <- tibble::rownames_to_column(hydrocarbon_degradation, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(hydrocarbon_degradation, by = c("bin_ID"))

hydrocarbon_degradation_NMDS2 <- ggplot(data.scores_genes, aes(x = NMDS2, y = hydrocarbon_degradation)) +
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "hydrocarbon degradation", shape = "Type")+
  stat_cor(method="spearman", size = 7)

corr <- cor.test(x=data.scores_genes$NMDS2, y=data.scores_genes$hydrocarbon_degradation, method = 'spearman')
corr

hydrocarbon_degradation_NMDS2

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/NMDS/hydrocarbon_degradation_NMDS.pdf",width=12,height=10) # Open a new pdf file

#view
hydrocarbon_degradation_NMDS2

dev.off()
#Save as PDF

CAZY <- metabolism_summary[metabolism_summary$header %like% "CAZY", ]
CAZY <- colSums(CAZY [ , c(7:1014)], na.rm=TRUE)

CAZY <- data.frame(CAZY)

library(tibble)
CAZY <- tibble::rownames_to_column(CAZY, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(CAZY, by = c("bin_ID"))

CAZY_NMDS2 <- ggplot(data.scores_genes, aes(x = NMDS2, y = CAZY)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "CAZY)", shape = "Type")+
  stat_cor(method="spearman", size = 7)

corr <- cor.test(x=data.scores_genes$NMDS1, y=data.scores_genes$CAZY, method = 'spearman')
corr

CAZY_NMDS2

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/NMDS/CAZY_NMDS.pdf",width=12,height=10) # Open a new pdf file

#view
CAZY_NMDS2

dev.off()
#Save as PDF

central_carbon <- metabolism_summary[metabolism_summary$header %like% "central carbon", ]
central_carbon <- colSums(central_carbon [ , c(7:1014)], na.rm=TRUE)

central_carbon <- data.frame(central_carbon)

library(tibble)
central_carbon <- tibble::rownames_to_column(central_carbon, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(central_carbon, by = c("bin_ID"))

central_carbon_NMDS2 <- ggplot(data.scores_genes, aes(x = NMDS2, y = central_carbon)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "central_carbon", shape = "Type")+
  stat_cor(method="spearman", size = 7)

corr <- cor.test(x=data.scores_genes$NMDS2, y=data.scores_genes$central_carbon, method = 'spearman')
corr

central_carbon_NMDS2

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/NMDS/central_carbon_NMDS.pdf",width=12,height=10) # Open a new pdf file

#view
central_carbon_NMDS2

dev.off()
#Save as PDF

pyruvate_metabolism <- metabolism_summary[metabolism_summary$header %like% "pyruvate metabolism", ]
pyruvate_metabolism <- colSums(pyruvate_metabolism [ , c(7:1014)], na.rm=TRUE)

pyruvate_metabolism <- data.frame(pyruvate_metabolism)

library(tibble)
pyruvate_metabolism <- tibble::rownames_to_column(pyruvate_metabolism, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(pyruvate_metabolism, by = c("bin_ID"))

pyruvate_metabolism_NMDS2 <- ggplot(data.scores_genes, aes(x = NMDS2, y = pyruvate_metabolism)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "pyruvate metabolism", shape = "Type")+
  stat_cor(method="spearman", size = 7)

corr <- cor.test(x=data.scores_genes$NMDS2, y=data.scores_genes$pyruvate_metabolism, method = 'spearman')
corr

sugar_utilization_woodcroft <- metabolism_summary[metabolism_summary$header %like% "sugar utilization (woodcroft)", ]
sugar_utilization_woodcroft <- colSums(sugar_utilization_woodcroft [ , c(7:1014)], na.rm=TRUE)

data.scores_genes$sugar_utilization_woodcroft = c(sugar_utilization_woodcroft)

ggplot(data.scores_genes, aes(x = NMDS2, y = sugar_utilization_woodcroft)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "sugar utilization (woodcroft)", shape = "Type")

#corr <- cor.test(x=data.scores_genes$NMDS1, y=data.scores_genes$sugar_utilization_woodcroft, method = 'spearman')
#corr

Vitamin_B12_transport_system <- metabolism_summary[metabolism_summary$header %like% "Vitamin B12 transport system", ]
Vitamin_B12_transport_system <- colSums(Vitamin_B12_transport_system [ , c(7:1014)], na.rm=TRUE)

Vitamin_B12_transport_system <- data.frame(Vitamin_B12_transport_system)

library(tibble)
Vitamin_B12_transport_system <- tibble::rownames_to_column(Vitamin_B12_transport_system, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(Vitamin_B12_transport_system, by = c("bin_ID"))

Vitamin_B12_transport_system_NMDS2 <- ggplot(data.scores_genes, aes(x = NMDS2, y = Vitamin_B12_transport_system)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "Vitamin B12 transport system", shape = "Type")+
  stat_cor(method="spearman", size = 7)

corr <- cor.test(x=data.scores_genes$NMDS2, y=data.scores_genes$Vitamin_B12_transport_system, method = 'spearman')
corr


C1 <- metabolism_summary[metabolism_summary$header %like% "C1", ]
C1 <- colSums(C1 [ , c(7:1014)], na.rm=TRUE)

C1 <- data.frame(C1)

library(tibble)
C1 <- tibble::rownames_to_column(C1, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(C1, by = c("bin_ID"))

C1_NMDS2 <- ggplot(data.scores_genes, aes(x = NMDS2, y = C1)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "C1", shape = "Type")+
  stat_cor(method="spearman", size = 7)

corr <- cor.test(x=data.scores_genes$NMDS2, y=data.scores_genes$C1, method = 'spearman')
corr


#C1_methane <- metabolism_summary[metabolism_summary$header %like% "C1_methane", ]
#C1_methane <- colSums(C1_methane [ , c(7:1014)], na.rm=TRUE)

#data.scores_genes$C1_methane = c(C1_methane)

ggplot(data.scores_genes, aes(x = NMDS2, y = C1_methane)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "C1-methane", shape = "Type")

#corr <- cor.test(x=data.scores_genes$NMDS1, y=data.scores_genes$C1_methane, method = 'spearman')
#corr


Electron_transport_Chain <- metabolism_summary[metabolism_summary$header %like% "Electron transport Chain", ]
Electron_transport_Chain <- colSums(Electron_transport_Chain [ , c(7:1014)], na.rm=TRUE)

Electron_transport_Chain <- data.frame(Electron_transport_Chain)

library(tibble)
Electron_transport_Chain <- tibble::rownames_to_column(Electron_transport_Chain, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(Electron_transport_Chain, by = c("bin_ID"))

Electron_transport_Chain_NMDS2 <- ggplot(data.scores_genes, aes(x = NMDS2, y = Electron_transport_Chain)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "Electron transport Chain", shape = "Type")+
  stat_cor(method="spearman", size = 7)

corr <- cor.test(x=data.scores_genes$NMDS2, y=data.scores_genes$Electron_transport_Chain, method = 'spearman')
corr

Electron_transport_Chain_NMDS2

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/NMDS/Electron_transport_Chain_NMDS.pdf",width=12,height=10) # Open a new pdf file

#view
Electron_transport_Chain_NMD2

dev.off()
#Save as PDF


Hydrogenases <- metabolism_summary[metabolism_summary$header %like% "Hydrogenases", ]
Hydrogenases <- colSums(Hydrogenases [ , c(7:1014)], na.rm=TRUE)

Hydrogenases <- data.frame(Hydrogenases)

library(tibble)
Hydrogenases <- tibble::rownames_to_column(Hydrogenases, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(Hydrogenases, by = c("bin_ID"))

ggplot(data.scores_genes, aes(x = NMDS2, y = Hydrogenases)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "Hydrogenases", shape = "Type")+
  stat_cor(method="spearman", size = 7)

corr <- cor.test(x=data.scores_genes$NMDS2, y=data.scores_genes$Hydrogenases, method = 'spearman')
corr


Metal_Reduction <- metabolism_summary[metabolism_summary$header %like% "Metal Reduction", ]
Metal_Reduction <- colSums(Metal_Reduction [ , c(7:1014)], na.rm=TRUE)

Metal_Reduction <- data.frame(Metal_Reduction)

library(tibble)
Metal_Reduction <- tibble::rownames_to_column(Metal_Reduction, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(Metal_Reduction, by = c("bin_ID"))

ggplot(data.scores_genes, aes(x = NMDS2, y = Metal_Reduction)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "Metal Reduction", shape = "Type")

corr <- cor.test(x=data.scores_genes$NMDS2, y=data.scores_genes$Metal_Reduction, method = 'spearman')
corr

Nitrogen <- metabolism_summary[metabolism_summary$header %like% "Nitrogen", ]
Nitrogen <- colSums(Nitrogen [ , c(7:1014)], na.rm=TRUE)

Nitrogen <- data.frame(Nitrogen)

library(tibble)
Nitrogen <- tibble::rownames_to_column(Nitrogen, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(Nitrogen, by = c("bin_ID"))

ggplot(data.scores_genes, aes(x = NMDS2, y = Nitrogen)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "Nitrogen", shape = "Type")+
  stat_cor(method="spearman", size = 7)

corr <- cor.test(x=data.scores_genes$NMDS1, y=data.scores_genes$Nitrogen, method = 'spearman')
corr

Oxygen <- metabolism_summary[metabolism_summary$header %like% "Oxygen", ]
Oxygen <- colSums(Oxygen [ , c(7:1014)], na.rm=TRUE)

Oxygen <- data.frame(Oxygen)

library(tibble)
Oxygen <- tibble::rownames_to_column(Oxygen, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(Oxygen, by = c("bin_ID"))

Oxygen_NMDS2 <- ggplot(data.scores_genes, aes(x = NMDS2, y = Oxygen)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "Oxygen", shape = "Type")+
  stat_cor(method="spearman", size = 7)

corr <- cor.test(x=data.scores_genes$NMDS2, y=data.scores_genes$Oxygen, method = 'spearman')
corr

Oxygen_NMDS2

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/NMDS/Oxygen_NMDS.pdf",width=12,height=10) # Open a new pdf file

#view
Oxygen_NMDS2

dev.off()
#Save as PDF


Photosynthesis <- metabolism_summary[metabolism_summary$header %like% "Photosynthesis", ]
Photosynthesis <- colSums(Photosynthesis [ , c(7:1014)], na.rm=TRUE)

Photosynthesis <- data.frame(Photosynthesis)

library(tibble)
Photosynthesis <- tibble::rownames_to_column(Photosynthesis, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(Photosynthesis, by = c("bin_ID"))

Photosynthesis_NMDS2 <- ggplot(data.scores_genes, aes(x = NMDS2, y = Photosynthesis)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "Photosynthesis", shape = "Type")+
  stat_cor(method="spearman", size = 7)

corr <- cor.test(x=data.scores_genes$NMDS2, y=data.scores_genes$Photosynthesis, method = 'spearman')
corr

Photosynthesis_NMDS2

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/NMDS/Photosynthesis_NMDS.pdf",width=12,height=10) # Open a new pdf file

#view
Photosynthesis_NMDS2

dev.off()
#Save as PDF


Sulfur <- metabolism_summary[metabolism_summary$header %like% "Sulfur", ]
Sulfur <- colSums(Sulfur [ , c(7:1014)], na.rm=TRUE)

Sulfur <- data.frame(Sulfur)

library(tibble)
Sulfur <- tibble::rownames_to_column(Sulfur, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(Sulfur, by = c("bin_ID"))

Sulfur_NMDS2 <- ggplot(data.scores_genes, aes(x = NMDS2, y = Sulfur)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "Sulfur", shape = "Type")+
  stat_cor(method="spearman", size = 7)

corr <- cor.test(x=data.scores_genes$NMDS2, y=data.scores_genes$Sulfur, method = 'spearman')
corr

Sulfur_NMDS2

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/NMDS/Sulfur_NMDS.pdf",width=12,height=10) # Open a new pdf file

#view
Sulfur_NMDS2

dev.off()
#Save as PDF

Amino_Acid <- metabolism_summary[metabolism_summary$header %like% "Amino Acid", ]
Amino_Acid <- colSums(Amino_Acid [ , c(7:1014)], na.rm=TRUE)

Amino_Acid <- data.frame(Amino_Acid)

library(tibble)
Amino_Acid <- tibble::rownames_to_column(Amino_Acid, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(Amino_Acid, by = c("bin_ID"))

Amino_Acid_NMDS2 <- ggplot(data.scores_genes, aes(x = NMDS2, y = Amino_Acid)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "Amino acid gene counts per bin", shape = "Type")+
  stat_cor(method="spearman", size = 7)

corr <- cor.test(x=data.scores_genes$NMDS2, y=data.scores_genes$Amino_Acid, method = 'spearman')
corr

Amino_Acid_NMDS2

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/NMDS/Amino_Acid_NMDS.pdf",width=12,height=10) # Open a new pdf file

#view
Amino_Acid_NMDS2

dev.off()
#Save as PDF

Peptidase <- metabolism_summary[metabolism_summary$header %like% "Peptidase", ]
Peptidase <- colSums(Peptidase [ , c(7:1014)], na.rm=TRUE)

Peptidase <- data.frame(Peptidase)

library(tibble)
Peptidase <- tibble::rownames_to_column(Peptidase, "bin_ID")

data.scores_genes <- data.scores_genes %>%
  join(Peptidase, by = c("bin_ID"))

peptidase_NMDS2 <- ggplot(data.scores_genes, aes(x = NMDS2, y = Peptidase)) + 
  geom_point(size = 4, aes(color = Phylum, shape = cluster)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS2", colour = "Phylum", y = "Peptidase", shape = "Type")+
  stat_cor(method="spearman", size = 7)

corr <- cor.test(x=data.scores_genes$NMDS2, y=data.scores_genes$Peptidase, method = 'spearman')
corr

peptidase_NMDS2

pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/NMDS/peptidase_NMDS.pdf",width=12,height=10) # Open a new pdf file

#view
peptidase_NMDS2

dev.off()
#Save as PDF

###For NMDS1 to different environmental parameters
environmental_params <- read_excel("/Volumes/micro-shared$/MoralesLab/Projects/MOTS/MOTS_MAGs/MOTS_BINS.xlsx", sheet = "MOTS_metadata")

colnames(environmental_params)[2] <- "Sample"

data.scores_env <- data.scores

data.scores_env <- merge(data.scores_env, environmental_params, by=c("Sample","Sample"))

colnames(data.scores_env)[19] <- "Temp"
colnames(data.scores_env)[20] <- "Salinity"
colnames(data.scores_env)[21] <- "Phosphate"
colnames(data.scores_env)[22] <- "Nitrate"
colnames(data.scores_env)[23] <- "Silicate"
colnames(data.scores_env)[24] <- "chlA"

ggplot(data.scores_env, aes(x = NMDS1, y = Temp)) + 
  geom_point(size = 4, aes(color = Phylum)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "Temp", shape = "Type")

corr <- cor.test(x=data.scores_env$NMDS1, y=data.scores_env$Temp, method = 'spearman')
corr

ggplot(data.scores_env, aes(x = NMDS1, y = Salinity)) + 
  geom_point(size = 4, aes(color = Phylum)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "Salinity", shape = "Type")

corr <- cor.test(x=data.scores_env$NMDS1, y=data.scores_env$Salinity, method = 'spearman')
corr

ggplot(data.scores_env, aes(x = NMDS1, y = Phosphate)) + 
  geom_point(size = 4, aes(color = Phylum)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "Phosphate", shape = "Type")

corr <- cor.test(x=data.scores_env$NMDS1, y=data.scores_env$Phosphate, method = 'spearman')
corr

ggplot(data.scores_env, aes(x = NMDS1, y = Nitrate)) + 
  geom_point(size = 4, aes(color = Phylum)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "Nitrate", shape = "Type")

corr <- cor.test(x=data.scores_env$NMDS1, y=data.scores_env$Nitrate, method = 'spearman')
corr

ggplot(data.scores_env, aes(x = NMDS1, y = Silicate)) + 
  geom_point(size = 4, aes(color = Phylum)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "Silicate", shape = "Type")

corr <- cor.test(x=data.scores_env$NMDS1, y=data.scores_env$Silicate, method = 'spearman')
corr

ggplot(data.scores_env, aes(x = NMDS1, y = chlA)) + 
  geom_point(size = 4, aes(color = Phylum)) +
  scale_color_manual(values = c("Chloroflexota" = "#C2DD9B", "Cyanobacteria" = "#669C4C", "Dadabacteria" = "#B4CCDF", "Fibrobacterota" = "#4A75AA", "Gemmatimonadota" = "#E3A29D", "Marinisomatota" = "#C04335", "Myxococcota" = "#ECC485", "Nitrospinota" = "#E18B46", "Planctomycetota" = "#4B9B7A", "Proteobacteria" = "#5F4190", "SAR324" = "#FFFFB0", "Verrucomicrobiota" = "#9C623C", "Acidobacteriota" = "#666666", "Actinobacteriota" = "#D43F88", "Bacteroidota" = "#7470AF", "Binatota" = "#CA6728", "Crenarchaeota" = "#C3B4D1", "Thermoplasmatota" = "lightgrey")) +
  My_Theme +
  labs(x = "NMDS1", colour = "Phylum", y = "chlA", shape = "Type")

corr <- cor.test(x=data.scores_env$NMDS1, y=data.scores_env$chlA, method = 'spearman')
corr

legend <- cowplot::get_legend(Amino_Acid_NMDS)

##NMDS1
pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/all_NMDS.pdf", width = 30, height = 40) # Open a new pdf file



grid.arrange(Amino_Acid_NMDS + theme(legend.position="none"), Sulfur_NMDS + theme(legend.position="none"), CAZY_NMDS+ theme(legend.position="none"), central_carbon_NMDS+ theme(legend.position="none"), Electron_transport_Chain_NMDS + theme(legend.position="none"), MISC_NMDS+ theme(legend.position="none"), peptidase_NMDS + theme(legend.position="none"), AMR_NMDS+ theme(legend.position="none"), Oxygen_NMDS + theme(legend.position="none"), aerobic_corrin_ring_synthesis_NMDS+ theme(legend.position="none"), Flagella_Structure_NMDS+ theme(legend.position="none"), ADO_CBL_synthesis_NMDS+ theme(legend.position="none"), Photosynthesis_NMDS + theme(legend.position="none"), hydrocarbon_degradation_NMDS+ theme(legend.position="none"),
             ncol=3, nrow=5)

dev.off()

###NMDS2
pdf("/Volumes/micro-shared$/MoralesLab/Manuscripts/MOTS_global_paper/temp_figs/all_NMDS2.pdf", width = 40, height = 40) # Open a new pdf file



grid.arrange(Amino_Acid_NMDS2 + theme(legend.position="none"), 
             Sulfur_NMDS2 + theme(legend.position="none"), 
             CAZY_NMDS2+ theme(legend.position="none"), 
             central_carbon_NMDS2+ theme(legend.position="none"), 
             Electron_transport_Chain_NMDS2 + theme(legend.position="none"), 
             MISC_NMDS2+ theme(legend.position="none"), 
             peptidase_NMDS2 + theme(legend.position="none"), 
             AMR_NMDS2+ theme(legend.position="none"), 
             Information_systems_NMDS2+ theme(legend.position="none"), 
             Oxygen_NMDS2 + theme(legend.position="none"), 
             aerobic_corrin_ring_synthesis_NMDS2+ theme(legend.position="none"), 
             Flagella_Structure_NMDS2+ theme(legend.position="none"), 
             ADO_CBL_synthesis_NMDS2+ theme(legend.position="none"), 
             Photosynthesis_NMDS2 + theme(legend.position="none"), 
             hydrocarbon_degradation_NMDS2+ theme(legend.position="none"), 
             pyruvate_metabolism_NMDS2+ theme(legend.position="none"), 
             C1_NMDS2+ theme(legend.position="none"), ncol=4, nrow=5)

dev.off()
