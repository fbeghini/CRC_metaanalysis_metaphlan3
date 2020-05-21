### --- Directory ------------------------------------------------------------------------------
setwd("/news/users/athomas/athomas_crc_reanalysis/CRC_metaanalysis_metaphlan3/")

### --- Libraries ------------------------------------------------------------------------------
library("textshape")
library(tidyr)
library(tibble)
library(dplyr)
library(readr)
library(stringr)
library("metafor")
library("meta")
library("plyr")
library("phyloseq") 
library("ExperimentHub")
library("Biobase")
library("curatedMetagenomicData")
library("reshape2")
library("ggplot2")
library("vegan")

### --- Functions ------------------------------------------------------------------------------
effect_size_calc <- function(vector1, vector2) {
  group1_mean <- mean(as.numeric(vector1))
  group2_mean <- mean(as.numeric(vector2))
  group1_sd <- sd(as.numeric(vector1))
  group2_sd <- sd(as.numeric(vector2))
  effect_size <- (group2_mean - group1_mean) / sqrt((group1_sd^2 + group2_sd^2) / 2)
  return(effect_size)
  
}

### --- Species abundance meta analysis mpa3 stat_q 0.1 ------------------------------------------------------------------------------
metadata <- read_tsv("tables/CRC_analysis_metadata_final_version.tsv", col_names = T) %>% 
  as.data.frame(row.names = FALSE, stringsAsFactors=FALSE)
metadata$study_condition <- gsub("CRC", "carcinoma", metadata$study_condition)
rownames(metadata)[metadata$dataset_name == "FengQ_2015"] <- as.character(metadata$sampleID[metadata$dataset_name == "FengQ_2015"])
rownames(metadata)[metadata$dataset_name == "ThomasAM_2019_c"] <- as.character(metadata$sampleID[metadata$dataset_name == "ThomasAM_2019_c"])
rownames(metadata)[metadata$dataset_name == "ZellerG_2014"] <- as.character(metadata$sampleID[metadata$dataset_name == "ZellerG_2014"])
rownames(metadata)[metadata$dataset_name == "YachidaS_2019"] <- as.character(metadata$sampleID[metadata$dataset_name == "YachidaS_2019"])
rownames(metadata)[metadata$dataset_name == "YuJ_2015"] <- as.character(metadata$sampleID[metadata$dataset_name == "YuJ_2015"])
rownames(metadata)[metadata$dataset_name == "ThomasAM_2019_a"] <- as.character(metadata$subjectID[metadata$dataset_name == "ThomasAM_2019_a"])
rownames(metadata)[metadata$dataset_name == "ThomasAM_2019_b"] <- as.character(metadata$subjectID[metadata$dataset_name == "ThomasAM_2019_b"])
rownames(metadata)[metadata$dataset_name == "VogtmannE_2016"] <- as.character(metadata$subjectID[metadata$dataset_name == "VogtmannE_2016"])
rownames(metadata)[metadata$dataset_name == "WirbelJ_2018"] <- as.character(metadata$subjectID[metadata$dataset_name == "WirbelJ_2018"])
rownames(metadata)[metadata$dataset_name == "GuptaA_2019"] <- as.character(metadata$sampleID[metadata$dataset_name == "GuptaA_2019"])

metaphlan3 <- read.table("tables/CRC_mpa3_statq01.tsv", stringsAsFactors = F, header = T, check.names = F, row.names = 1, as.is = T)
colnames(metaphlan3) <- gsub("_profile", "", colnames(metaphlan3))
colnames(metaphlan3)[grep("LILT", colnames(metaphlan3))] <- unlist(lapply( colnames(metaphlan3)[grep("LILT", colnames(metaphlan3))], function(x) unlist(strsplit(x, split = "_"))[[2]]))
colnames(metaphlan3)[grep("CRC_MR", colnames(metaphlan3))] <- unlist(lapply( colnames(metaphlan3)[grep("CRC_MR", colnames(metaphlan3))], function(x) unlist(strsplit(x, split = "_"))[[3]]))
metaphlan3 <- as.data.frame(metaphlan3)
metaphlan3_2 <- apply(metaphlan3, 2, as.numeric)
rownames(metaphlan3_2) <- rownames(metaphlan3)
metaphlan3 <- metaphlan3_2
metaphlan3 <- apply(metaphlan3, 2, function(x) asin(sqrt(x/100)))
metaphlan3 <- as.matrix(metaphlan3)
common.ids <- intersect(rownames(metadata), colnames(metaphlan3))
metadata <- metadata[common.ids,]
metaphlan3 <- metaphlan3[,common.ids]

metadata$sampleID <- rownames(metadata)
labelDescriptors <- data.frame(labelDescription=c(colnames(metadata)))
metadata_annotated <- AnnotatedDataFrame(data = metadata, varMetadata = labelDescriptors, dimLabels= c("sampleNames", "varLabels"))
metaphlan3 <- ExpressionSet(phenoData = metadata_annotated, assayData = metaphlan3)
metaphlan3 <- ExpressionSet2phyloseq(metaphlan3)

#Prune taxa
metaphlan_species <- subset_taxa(metaphlan3, !is.na(Species) & is.na(Strain))

# Remove 2 samples that had unclassified reads
metaphlan_species <- subset_samples(physeq = metaphlan_species, !(sample_names(metaphlan_species) %in% names(which(colSums(otu_table(metaphlan_species)) == 0))))
sample_data(metaphlan_species)$study_condition <- gsub("CRC", "carcinoma", sample_data(metaphlan_species)$study_condition)
sample_data(metaphlan_species)$dataset_name <- gsub("ThomasAM_2019_c", "YachidaS_2019", sample_data(metaphlan_species)$dataset_name)
meta_analysis_results <- matrix(NA, ncol = 18, nrow=nrow(otu_table(metaphlan_species)))
colnames(meta_analysis_results) <- c("Species", "Number_Samples", "I^2", "P-value Fit", "Pvalue", "Standard Error", "Coefficient", "Confidence Interval Lower Limit", "Confidence Interval Upper Limit", "ZellerG_2014", "FengQ_2015", "YuJ_2015", "ThomasAM_2019_a", "ThomasAM_2019_b", "VogtmannE_2016","WirbelJ_2019", "YachidaS_2019", 'GuptaA_2019')
meta_analysis_results <- as.data.frame(meta_analysis_results)
for (k in 1:nrow(otu_table(metaphlan_species))) {
  species <- otu_table(metaphlan_species)[k,]
  species_table <- matrix(NA, ncol = 8, nrow = 10)
  colnames(species_table) <- c("Study", "Dataset", "SampleSize_Cancer", "Mean_Cancer", "SD_Cancer", "SampleSize_Control", "Mean_Control", "SD_Control")
  species_table <- as.data.frame(species_table)
  datasets <- as.character(unique(sample_data(metaphlan_species)$dataset_name))
  meta_analysis_results$VogtmannE_2016[k] <- wilcox.test(as.numeric(species[,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "carcinoma" & sample_data(metaphlan_species)$dataset_name == "VogtmannE_2016"]]),
                                                         as.numeric(species[,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "control" & sample_data(metaphlan_species)$dataset_name == "VogtmannE_2016"]]))$p.value
  meta_analysis_results$ThomasAM_2019_b[k] <- wilcox.test(as.numeric(species[,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "carcinoma" & sample_data(metaphlan_species)$dataset_name == "ThomasAM_2019_b"]]),
                                                          as.numeric(species[,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "control" & sample_data(metaphlan_species)$dataset_name == "ThomasAM_2019_b"]]))$p.value
  meta_analysis_results$ThomasAM_2019_a[k] <- wilcox.test(as.numeric(species[,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "carcinoma" & sample_data(metaphlan_species)$dataset_name == "ThomasAM_2019_a"]]),
                                                          as.numeric(species[,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "control" & sample_data(metaphlan_species)$dataset_name == "ThomasAM_2019_a"]]))$p.value
  meta_analysis_results$YuJ_2015[k] <- wilcox.test(as.numeric(species[,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "carcinoma" & sample_data(metaphlan_species)$dataset_name == "YuJ_2015"]]),
                                                   as.numeric(species[,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "control" & sample_data(metaphlan_species)$dataset_name == "YuJ_2015"]]))$p.value
  meta_analysis_results$FengQ_2015[k] <- wilcox.test(as.numeric(species[,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "carcinoma" & sample_data(metaphlan_species)$dataset_name == "FengQ_2015"]]),
                                                     as.numeric(species[,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "control" & sample_data(metaphlan_species)$dataset_name == "FengQ_2015"]]))$p.value
  meta_analysis_results$ZellerG_2014[k] <- wilcox.test(as.numeric(species[,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "carcinoma" & sample_data(metaphlan_species)$dataset_name == "ZellerG_2014"]]),
                                                       as.numeric(species[,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "control" & sample_data(metaphlan_species)$dataset_name == "ZellerG_2014"]]))$p.value
  meta_analysis_results$WirbelJ_2019[k] <- wilcox.test(as.numeric(species[,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "carcinoma" & sample_data(metaphlan_species)$dataset_name == "WirbelJ_2018"]]),
                                                       as.numeric(species[,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "control" & sample_data(metaphlan_species)$dataset_name == "WirbelJ_2018"]]))$p.value
  meta_analysis_results$YachidaS_2019[k] <- wilcox.test(as.numeric(species[,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "carcinoma" & sample_data(metaphlan_species)$dataset_name == "YachidaS_2019"]]),
                                                        as.numeric(species[,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "control" & sample_data(metaphlan_species)$dataset_name == "YachidaS_2019"]]))$p.value
  meta_analysis_results$GuptaA_2019[k] <- wilcox.test(as.numeric(species[,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "carcinoma" & sample_data(metaphlan_species)$dataset_name == "GuptaA_2019"]]),
                                                      as.numeric(species[,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "control" & sample_data(metaphlan_species)$dataset_name == "GuptaA_2019"]]))$p.value
  for (i in 1:length(datasets)) {
    dataset <- datasets[i]
    species_table$Study[i] <- i
    species_table$Dataset[i] <- dataset
    species_table$SampleSize_Cancer[i] <- length(which(sample_data(metaphlan_species)$study_condition == "carcinoma" & sample_data(metaphlan_species)$dataset_name == dataset))
    species_table$Mean_Cancer[i] <- mean(species[,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "carcinoma" & sample_data(metaphlan_species)$dataset_name == dataset]])
    species_table$SD_Cancer[i] <- sd(species[,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "carcinoma" & sample_data(metaphlan_species)$dataset_name == dataset]])
    species_table$SampleSize_Control[i] <- length(which(sample_data(metaphlan_species)$study_condition == "control" & sample_data(metaphlan_species)$dataset_name == dataset))
    species_table$Mean_Control[i] <- mean(species[,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "control" & sample_data(metaphlan_species)$dataset_name == dataset]])
    species_table$SD_Control[i] <- sd(species[,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "control" & sample_data(metaphlan_species)$dataset_name == dataset]])
  }
  dat1 <- escalc(measure="SMD", m1i=Mean_Cancer, sd1i=SD_Cancer, n1i=SampleSize_Cancer,
                 m2i=Mean_Control, sd2i=SD_Control, n2i=SampleSize_Control, data=species_table, vtype="UB")
  if (length(which((dat1$Mean_Cancer == dat1$Mean_Control) == TRUE)) > 3) {
    meta_analysis_results$Species[k] <- rownames(otu_table(metaphlan_species))[k]
    meta_analysis_results$Number_Samples[k] <- 0  
    next()
  }
  else{
    res1 <- rma(yi, vi, data=dat1, control=list(stepadj=0.5, maxiter = 10000), verbose=TRUE, digits=5)
    meta_analysis_results$Species[k] <- rownames(otu_table(metaphlan_species))[k]
    meta_analysis_results$Number_Samples[k] <- sum(res1$ni)
    meta_analysis_results$`I^2`[k] <- res1$I2
    meta_analysis_results$`P-value Fit`[k] <- res1$QEp
    meta_analysis_results$Pvalue[k] <- res1$pval
    meta_analysis_results$`Standard Error`[k] <- res1$b
    meta_analysis_results$Coefficient[k] <- res1$b
    meta_analysis_results$`Confidence Interval Lower Limit`[k] <- res1$ci.lb
    meta_analysis_results$`Confidence Interval Upper Limit`[k] <- res1$ci.ub
  }
}
meta_analysis_results$PvalueAdjusted <- p.adjust(meta_analysis_results$Pvalue, method = "fdr")
meta_analysis_results <- meta_analysis_results[order(abs(meta_analysis_results$Coefficient), decreasing = T), ]
write.table(meta_analysis_results, "results/meta_analysis_results_metaphlan3_arcsine_all_datasets.txt", col.names = T, row.names = F, quote = F, sep = "\t")

rownames(meta_analysis_results) <- meta_analysis_results$Species
meta_analysis_results <- meta_analysis_results[!(is.na(meta_analysis_results$`I^2`)),]

meta_analysis_results_filt <- meta_analysis_results[meta_analysis_results$PvalueAdjusted < 0.05,]
meta_analysis_results_filt <- meta_analysis_results_filt[meta_analysis_results_filt$`P-value Fit` > 0.05,]
write.table(meta_analysis_results, "results/meta_analysis_results_filt_metaphlan3_arcsine_all_datasets.txt", col.names = T, row.names = F, quote = F, sep = "\t")

meta_analysis_results_filt <- meta_analysis_results_filt[1:20,]
species <- meta_analysis_results_filt$Species
write.table(meta_analysis_results_filt, "results/meta_analysis_results_top20_metaphlan3_arcsine_all_datasets.txt", col.names = T, row.names = F, quote = F, sep = "\t")

species_meta <- otu_table(metaphlan_species)[species,]
effect_sizes <- matrix(NA, ncol = 10, nrow=nrow(species_meta))
colnames(effect_sizes) <- c("Species", "ZellerG_2014", "FengQ_2015", "YuJ_2015", "ThomasAM_2019_b", "ThomasAM_2019_b", "VogtmannE_2016", "WirbelJ_2018", "YachidaS_2019", 'GuptaA_2019')
effect_sizes <- as.data.frame(effect_sizes)
for (i in 1:nrow(species_meta)) {
  effect_sizes$Species[i] <- rownames(species_meta)[i]
  effect_sizes$VogtmannE_2016[i] <- effect_size_calc(as.numeric(species_meta[i,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "carcinoma" & sample_data(metaphlan_species)$dataset_name == "VogtmannE_2016"]]), as.numeric(species_meta[i,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "control" & sample_data(metaphlan_species)$dataset_name == "VogtmannE_2016"]]))
  effect_sizes$ThomasAM_2019_a[i] <- effect_size_calc(as.numeric(species_meta[i,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "carcinoma" & sample_data(metaphlan_species)$dataset_name == "ThomasAM_2019_a"]]), as.numeric(species_meta[i,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "control" & sample_data(metaphlan_species)$dataset_name == "ThomasAM_2019_a"]]))
  effect_sizes$ThomasAM_2019_b[i] <- effect_size_calc(as.numeric(species_meta[i,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "carcinoma" & sample_data(metaphlan_species)$dataset_name == "ThomasAM_2019_b"]]), as.numeric(species_meta[i,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "control" & sample_data(metaphlan_species)$dataset_name == "ThomasAM_2019_b"]]))
  effect_sizes$YuJ_2015[i] <- effect_size_calc(as.numeric(species_meta[i,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "carcinoma" & sample_data(metaphlan_species)$dataset_name == "YuJ_2015"]]), as.numeric(species_meta[i,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "control" & sample_data(metaphlan_species)$dataset_name == "YuJ_2015"]]))
  effect_sizes$FengQ_2015[i] <- effect_size_calc(as.numeric(species_meta[i,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "carcinoma" & sample_data(metaphlan_species)$dataset_name == "FengQ_2015"]]), as.numeric(species_meta[i,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "control" & sample_data(metaphlan_species)$dataset_name == "FengQ_2015"]]))
  effect_sizes$ZellerG_2014[i] <- effect_size_calc(as.numeric(species_meta[i,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "carcinoma" & sample_data(metaphlan_species)$dataset_name == "ZellerG_2014"]]),as.numeric(species_meta[i,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "control" & sample_data(metaphlan_species)$dataset_name == "ZellerG_2014"]]))
  effect_sizes$YachidaS_2019[i] <- effect_size_calc(as.numeric(species_meta[i,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "carcinoma" & sample_data(metaphlan_species)$dataset_name == "YachidaS_2019"]]),as.numeric(species_meta[i,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "control" & sample_data(metaphlan_species)$dataset_name == "YachidaS_2019"]]))
  effect_sizes$GuptaA_2019[i] <- effect_size_calc(as.numeric(species_meta[i,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "carcinoma" & sample_data(metaphlan_species)$dataset_name == "GuptaA_2019"]]),as.numeric(species_meta[i,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "control" & sample_data(metaphlan_species)$dataset_name == "GuptaA_2019"]]))
  effect_sizes$WirbelJ_2018[i] <- effect_size_calc(as.numeric(species_meta[i,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "carcinoma" & sample_data(metaphlan_species)$dataset_name == "WirbelJ_2018"]]),as.numeric(species_meta[i,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$study_condition == "control" & sample_data(metaphlan_species)$dataset_name == "WirbelJ_2018"]]))
}

effect_sizes$Species <- gsub("s__", "", effect_sizes$Species)
effect_sizes$Species <- gsub("_", " ", effect_sizes$Species)
rownames(meta_analysis_results_filt) <- gsub("s__", "", rownames(meta_analysis_results_filt))
rownames(meta_analysis_results_filt) <- gsub("_", " ", rownames(meta_analysis_results_filt))
effect_sizes$Species <- factor(effect_sizes$Species, levels = rev(c(rownames(meta_analysis_results_filt)[order(meta_analysis_results_filt$Coefficient, decreasing = T)])))

effect_sizes[,-1] <- effect_sizes[,-1] * -1

effect_sizes$RandomEffects <- meta_analysis_results_filt$Coefficient
effect_sizes <- melt(effect_sizes)
effect_sizes$CI_Low <- NA 
effect_sizes$CI_UP <- NA
effect_sizes[effect_sizes$variable == "RandomEffects","CI_Low"] <- meta_analysis_results_filt$`Confidence Interval Lower Limit`
effect_sizes[effect_sizes$variable == "RandomEffects","CI_UP"] <- meta_analysis_results_filt$`Confidence Interval Upper Limit`
effect_sizes$variable <- revalue(effect_sizes$variable, replace = c("WirbelJ_2018" = "WirbelJ_2019"))
effect_sizes$variable <- factor(effect_sizes$variable, levels = c("RandomEffects", "ThomasAM_2019_a", "ThomasAM_2019_b", "YachidaS_2019", "GuptaA_2019", "WirbelJ_2019", "HanniganGD_2017", "VogtmannE_2016", "FengQ_2015", "YuJ_2015", "ZellerG_2014"))
effect_sizes$variable <- revalue(effect_sizes$variable, replace = c("RandomEffects" = "Random Effects"))

p1 <- ggplot(effect_sizes, aes(Species, value, shape = variable, ymin=CI_Low, ymax=CI_UP)) + 
  geom_point(size = 2.5) +
  geom_linerange() + 
  scale_shape_manual(values= rev(c(3,4,8,0, 13, 11, 12, 2,5,19))) +
  coord_flip() +
  theme_minimal() +
  geom_hline(aes(yintercept = 0), color = "gray", linetype = "dashed") +
  theme(axis.text.y = element_text(face = "italic", size = 8), axis.title.y = element_blank(), axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 12), legend.title = element_blank(),  
        legend.key.size = unit(0.7, "cm"), legend.position=c(0.1,0.75), legend.text = element_text(size = 7.5)) + 
  ylab("Effect size") 
ggsave("figures/meta_analysis_mpa3_statq01_top20.svg", p1, width = 10, height = 8, device = "svg")

## ----------------- Species richness mpa3 rarefied stat_q 0.2 ----------------------------------------------------------------------------------------------------------------------------------------
metadata <- read_tsv("tables/CRC_analysis_metadata_final_version.tsv", col_names = T) %>% 
  as.data.frame(row.names = FALSE, stringsAsFactors=FALSE)
metadata$study_condition <- gsub("CRC", "carcinoma", metadata$study_condition)
rownames(metadata)[metadata$dataset_name == "FengQ_2015"] <- as.character(metadata$sampleID[metadata$dataset_name == "FengQ_2015"])
rownames(metadata)[metadata$dataset_name == "ThomasAM_2019_c"] <- as.character(metadata$sampleID[metadata$dataset_name == "ThomasAM_2019_c"])
rownames(metadata)[metadata$dataset_name == "ZellerG_2014"] <- as.character(metadata$sampleID[metadata$dataset_name == "ZellerG_2014"])
rownames(metadata)[metadata$dataset_name == "YachidaS_2019"] <- as.character(metadata$sampleID[metadata$dataset_name == "YachidaS_2019"])
rownames(metadata)[metadata$dataset_name == "YuJ_2015"] <- as.character(metadata$sampleID[metadata$dataset_name == "YuJ_2015"])
rownames(metadata)[metadata$dataset_name == "ThomasAM_2019_a"] <- as.character(metadata$subjectID[metadata$dataset_name == "ThomasAM_2019_a"])
rownames(metadata)[metadata$dataset_name == "ThomasAM_2019_b"] <- as.character(metadata$subjectID[metadata$dataset_name == "ThomasAM_2019_b"])
rownames(metadata)[metadata$dataset_name == "VogtmannE_2016"] <- as.character(metadata$subjectID[metadata$dataset_name == "VogtmannE_2016"])
rownames(metadata)[metadata$dataset_name == "WirbelJ_2018"] <- as.character(metadata$subjectID[metadata$dataset_name == "WirbelJ_2018"])
rownames(metadata)[metadata$dataset_name == "GuptaA_2019"] <- as.character(metadata$sampleID[metadata$dataset_name == "GuptaA_2019"])

metaphlan3 <- read.table("tables/CRC_mpa3_rarefied_0.1_percentile_statq_02.tsv", stringsAsFactors = F, header = T, row.names = 1, as.is = T, check.names = F)
colnames(metaphlan3) <- gsub("_rarefied_0.1_percentile_profile", "", colnames(metaphlan3))
colnames(metaphlan3)[grep("LILT", colnames(metaphlan3))] <- unlist(lapply( colnames(metaphlan3)[grep("LILT", colnames(metaphlan3))], function(x) unlist(strsplit(x, split = "_"))[[2]]))
colnames(metaphlan3)[grep("CRC_MR", colnames(metaphlan3))] <- unlist(lapply( colnames(metaphlan3)[grep("CRC_MR", colnames(metaphlan3))], function(x) unlist(strsplit(x, split = "_"))[[3]]))
common.ids <- intersect(rownames(metadata), colnames(metaphlan3))
metadata <- metadata[common.ids,]
metaphlan3 <- metaphlan3[,common.ids]
metaphlan3 <- as.matrix(metaphlan3)

labelDescriptors <- data.frame(labelDescription=c(colnames(metadata)))
metadata_annotated <- AnnotatedDataFrame(data = metadata, varMetadata = labelDescriptors, dimLabels= c("sampleNames", "varLabels"))
metaphlan3 <- ExpressionSet(phenoData = metadata_annotated, assayData = metaphlan3)
metaphlan3 <- ExpressionSet2phyloseq(metaphlan3)

#Prune taxa
metaphlan_species <- subset_taxa(metaphlan3, !is.na(Species) & is.na(Strain))

table <- otu_table(metaphlan_species)
species_diversity <- data.frame(colnames(table))
samples <- colnames(table)
rownames(species_diversity) <- colnames(table)
species_diversity$Observed <- NA

for (sample in samples) {
  number_species <- length(which(table[,sample] > 0))
  species_diversity[sample,"Observed"] <- number_species
}

species_diversity$Group <- sample_data(metaphlan_species)$study_condition
species_diversity$Dataset <- sample_data(metaphlan_species)$dataset_name
species_diversity$Group <- revalue(species_diversity$Group, replace = c("control" = "Control","carcinoma" = "Carcinoma"))
species_diversity$Group <- factor(species_diversity$Group, levels = c("Control", "Carcinoma"))
species_diversity$Dataset <- revalue(species_diversity$Dataset, replace =c("ThomasAM_2019_c" = "YachidaS_2019"))
species_diversity$Dataset <- factor(species_diversity$Dataset, levels=c("YuJ_2015", "ThomasAM_2019_b", "FengQ_2015", "VogtmannE_2016", "ZellerG_2014", "YachidaS_2019", "WirbelJ_2018", "ThomasAM_2019_a", "GuptaA_2019"))

df_richness <- data.frame(Dataset=c(unique(as.character(species_diversity$Dataset))), pvalue = rep(NA, 9), y = c(rep(3, 9)))
rownames(df_richness) <- unique(as.character(species_diversity$Dataset))
for (i in unique(as.character(species_diversity$Dataset)) ) {
  df_richness[i,"pvalue"] <- round(wilcox.test(species_diversity$Observed[species_diversity$Group == "Carcinoma" & species_diversity$Dataset == i],
                                               species_diversity$Observed[species_diversity$Group == "Control" & species_diversity$Dataset == i])$p.value, digits = 8)
}
df_richness$Dataset <- factor(df_richness$Dataset, levels=c("YuJ_2015", "ThomasAM_2019_b", "FengQ_2015", "VogtmannE_2016", "ZellerG_2014", "YachidaS_2019", "WirbelJ_2018", "ThomasAM_2019_a", "GuptaA_2019"))
df_richness <- df_richness[c("YuJ_2015", "ThomasAM_2019_b", "FengQ_2015", "VogtmannE_2016", "ZellerG_2014", "YachidaS_2019", "WirbelJ_2018", "ThomasAM_2019_a", "GuptaA_2019"),]
df_richness$x <- c(1:9)
df_richness$pvalue <-formatC(df_richness$pvalue, format = "e", digits = 1)
df_richness$pvalue <- paste("P = ", df_richness$pvalue, sep = "")

p1 <- ggplot(species_diversity, aes(Dataset, Observed, fill = Group, colour = Group)) +
  geom_boxplot(aes(colour = Group), outlier.colour = NA, position = "dodge", fill="white") +
  geom_point(position=position_jitterdodge(dodge.width=0.8), alpha = 0.4, size = 0.5) +     
  theme_bw() +
  ylab("Number of species") +
  theme(axis.title.x= element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_text(size = 8), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 7, angle = 45, hjust = 1), legend.key.size = unit(0.3, "cm"), legend.position=c(0.25,0.2), legend.direction = "horizontal", legend.title = element_blank(), legend.text = element_text(size = 7), legend.background = element_rect(fill="transparent")) +
  scale_y_continuous(breaks = c(0, 50, 100, 150, 200), limits = c(0, 200), labels = c(0, 50, 100, 150, 200)) +
  scale_color_manual(labels = c("Control","Carcinoma"), values = c("#6F835A", "#aa483b")) +
  geom_text(data=df_richness, aes(x, y, label=pvalue), color = "black", size = 4, inherit.aes= F)
ggsave("figures/species_richness_mpa3_rarefied_statq02.svg", p1, width = 10, height = 8, device = "svg")

## ----------------- Species richness mpa3 rarefied stat_q 0.1 ----------------------------------------------------------------------------------------------------------------------------------------
metadata <- read_tsv("tables/CRC_analysis_metadata_final_version.tsv", col_names = T) %>% 
  as.data.frame(row.names = FALSE, stringsAsFactors=FALSE)
metadata$study_condition <- gsub("CRC", "carcinoma", metadata$study_condition)
rownames(metadata)[metadata$dataset_name == "FengQ_2015"] <- as.character(metadata$sampleID[metadata$dataset_name == "FengQ_2015"])
rownames(metadata)[metadata$dataset_name == "ThomasAM_2019_c"] <- as.character(metadata$sampleID[metadata$dataset_name == "ThomasAM_2019_c"])
rownames(metadata)[metadata$dataset_name == "ZellerG_2014"] <- as.character(metadata$sampleID[metadata$dataset_name == "ZellerG_2014"])
rownames(metadata)[metadata$dataset_name == "YachidaS_2019"] <- as.character(metadata$sampleID[metadata$dataset_name == "YachidaS_2019"])
rownames(metadata)[metadata$dataset_name == "YuJ_2015"] <- as.character(metadata$sampleID[metadata$dataset_name == "YuJ_2015"])
rownames(metadata)[metadata$dataset_name == "ThomasAM_2019_a"] <- as.character(metadata$subjectID[metadata$dataset_name == "ThomasAM_2019_a"])
rownames(metadata)[metadata$dataset_name == "ThomasAM_2019_b"] <- as.character(metadata$subjectID[metadata$dataset_name == "ThomasAM_2019_b"])
rownames(metadata)[metadata$dataset_name == "VogtmannE_2016"] <- as.character(metadata$subjectID[metadata$dataset_name == "VogtmannE_2016"])
rownames(metadata)[metadata$dataset_name == "WirbelJ_2018"] <- as.character(metadata$subjectID[metadata$dataset_name == "WirbelJ_2018"])
rownames(metadata)[metadata$dataset_name == "GuptaA_2019"] <- as.character(metadata$sampleID[metadata$dataset_name == "GuptaA_2019"])

metaphlan3 <- read.table("tables/CRC_mpa3_rarefied_0.1_percentile_statq_01.tsv", stringsAsFactors = F, header = T, row.names = 1, as.is = T, check.names = F)
colnames(metaphlan3) <- gsub("_rarefied_0.1_percentile_profile_statq_01", "", colnames(metaphlan3))
colnames(metaphlan3)[grep("LILT", colnames(metaphlan3))] <- unlist(lapply( colnames(metaphlan3)[grep("LILT", colnames(metaphlan3))], function(x) unlist(strsplit(x, split = "_"))[[2]]))
colnames(metaphlan3)[grep("CRC_MR", colnames(metaphlan3))] <- unlist(lapply( colnames(metaphlan3)[grep("CRC_MR", colnames(metaphlan3))], function(x) unlist(strsplit(x, split = "_"))[[3]]))
common.ids <- intersect(rownames(metadata), colnames(metaphlan3))
metadata <- metadata[common.ids,]
metaphlan3 <- metaphlan3[,common.ids]
metaphlan3 <- as.matrix(metaphlan3)

labelDescriptors <- data.frame(labelDescription=c(colnames(metadata)))
metadata_annotated <- AnnotatedDataFrame(data = metadata, varMetadata = labelDescriptors, dimLabels= c("sampleNames", "varLabels"))
metaphlan3 <- ExpressionSet(phenoData = metadata_annotated, assayData = metaphlan3)
metaphlan3 <- ExpressionSet2phyloseq(metaphlan3)

#Prune taxa
metaphlan_species <- subset_taxa(metaphlan3, !is.na(Species) & is.na(Strain))

table <- otu_table(metaphlan_species)
species_diversity <- data.frame(colnames(table))
samples <- colnames(table)
rownames(species_diversity) <- colnames(table)
species_diversity$Observed <- NA

for (sample in samples) {
  number_species <- length(which(table[,sample] > 0))
  species_diversity[sample,"Observed"] <- number_species
}

species_diversity$Group <- sample_data(metaphlan_species)$study_condition
species_diversity$Dataset <- sample_data(metaphlan_species)$dataset_name
species_diversity$Group <- revalue(species_diversity$Group, replace = c("control" = "Control","carcinoma" = "Carcinoma"))
species_diversity$Group <- factor(species_diversity$Group, levels = c("Control", "Carcinoma"))
species_diversity$Dataset <- revalue(species_diversity$Dataset, replace =c("ThomasAM_2019_c" = "YachidaS_2019", "CM_lilt" = "ThomasAM_2019_a", "CM_rescignocrc" = "ThomasAM_2019_b"))
species_diversity$Dataset <- factor(species_diversity$Dataset, levels=c("YuJ_2015", "ThomasAM_2019_b", "FengQ_2015", "VogtmannE_2016", "ZellerG_2014", "YachidaS_2019", "WirbelJ_2018", "ThomasAM_2019_a", "GuptaA_2019"))

df_richness <- data.frame(Dataset=c(unique(as.character(species_diversity$Dataset))), pvalue = rep(NA, 9), y = c(rep(3, 9)))
rownames(df_richness) <- unique(as.character(species_diversity$Dataset))
for (i in unique(as.character(species_diversity$Dataset)) ) {
  df_richness[i,"pvalue"] <- round(wilcox.test(species_diversity$Observed[species_diversity$Group == "Carcinoma" & species_diversity$Dataset == i],
                                               species_diversity$Observed[species_diversity$Group == "Control" & species_diversity$Dataset == i])$p.value, digits = 8)
}
df_richness$Dataset <- factor(df_richness$Dataset, levels=c("YuJ_2015", "ThomasAM_2019_b", "FengQ_2015", "VogtmannE_2016", "ZellerG_2014", "YachidaS_2019", "WirbelJ_2018", "ThomasAM_2019_a", "GuptaA_2019"))
df_richness <- df_richness[c("YuJ_2015", "ThomasAM_2019_b", "FengQ_2015", "VogtmannE_2016", "ZellerG_2014", "YachidaS_2019", "WirbelJ_2018", "ThomasAM_2019_a", "GuptaA_2019"),]
df_richness$x <- c(1:9)
df_richness$pvalue <-formatC(df_richness$pvalue, format = "e", digits = 1)
df_richness$pvalue <- paste("P = ", df_richness$pvalue, sep = "")

p1 <- ggplot(species_diversity, aes(Dataset, Observed, fill = Group, colour = Group)) +
  geom_boxplot(aes(colour = Group), outlier.colour = NA, position = "dodge", fill="white") +
  geom_point(position=position_jitterdodge(dodge.width=0.8), alpha = 0.4, size = 0.5) +     
  theme_bw() +
  ylab("Number of species") +
  theme(axis.title.x= element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_text(size = 8), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 7, angle = 45, hjust = 1), legend.key.size = unit(0.3, "cm"), legend.position=c(0.25,0.2), legend.direction = "horizontal", legend.title = element_blank(), legend.text = element_text(size = 7), legend.background = element_rect(fill="transparent")) +
  scale_y_continuous(breaks = c(0, 50, 100, 150, 200, 250), limits = c(0, 250), labels = c(0, 50, 100, 150, 200, 250)) +
  scale_color_manual(labels = c("Control","Carcinoma"), values = c("#6F835A", "#aa483b")) +
  geom_text(data=df_richness, aes(x, y, label=pvalue), color = "black", size = 4, inherit.aes= F)
ggsave("figures/species_richness_mpa3_rarefied_statq01.svg", p1, width = 10, height = 8, device = "svg")

## ----------------- species richness mpa2 rarefied ----------------------------------------------------------------------------------------------------------------------------------------
metadata <- read_tsv("tables/CRC_analysis_metadata_final_version.tsv", col_names = T) %>% 
  as.data.frame(row.names = FALSE, stringsAsFactors=FALSE)
metadata$study_condition <- gsub("CRC", "carcinoma", metadata$study_condition)
rownames(metadata)[metadata$dataset_name == "FengQ_2015"] <- as.character(metadata$sampleID[metadata$dataset_name == "FengQ_2015"])
rownames(metadata)[metadata$dataset_name == "ThomasAM_2019_c"] <- as.character(metadata$sampleID[metadata$dataset_name == "ThomasAM_2019_c"])
rownames(metadata)[metadata$dataset_name == "ZellerG_2014"] <- as.character(metadata$sampleID[metadata$dataset_name == "ZellerG_2014"])
rownames(metadata)[metadata$dataset_name == "YachidaS_2019"] <- as.character(metadata$sampleID[metadata$dataset_name == "YachidaS_2019"])
rownames(metadata)[metadata$dataset_name == "YuJ_2015"] <- as.character(metadata$sampleID[metadata$dataset_name == "YuJ_2015"])
rownames(metadata)[metadata$dataset_name == "ThomasAM_2019_a"] <- as.character(metadata$subjectID[metadata$dataset_name == "ThomasAM_2019_a"])
rownames(metadata)[metadata$dataset_name == "ThomasAM_2019_b"] <- as.character(metadata$subjectID[metadata$dataset_name == "ThomasAM_2019_b"])
rownames(metadata)[metadata$dataset_name == "VogtmannE_2016"] <- as.character(metadata$subjectID[metadata$dataset_name == "VogtmannE_2016"])
rownames(metadata)[metadata$dataset_name == "WirbelJ_2018"] <- as.character(metadata$subjectID[metadata$dataset_name == "WirbelJ_2018"])
rownames(metadata)[metadata$dataset_name == "GuptaA_2019"] <- as.character(metadata$sampleID[metadata$dataset_name == "GuptaA_2019"])

metaphlan2 <- read.table("tables/CRC_mpa2_rarefied_0.1_percentile.tsv", stringsAsFactors = F, header = T, check.names = F, row.names = 1)
colnames(metaphlan2) <- gsub("_metaphlan2.7_rarefied_0.1_percentile_profile", "", colnames(metaphlan2))
colnames(metaphlan2) <- gsub("_profile", "", colnames(metaphlan2))
colnames(metaphlan2)[grep("LILT", colnames(metaphlan2))] <- unlist(lapply( colnames(metaphlan2)[grep("LILT", colnames(metaphlan2))], function(x) unlist(strsplit(x, split = "_"))[[2]]))
colnames(metaphlan2)[grep("CRC_MR", colnames(metaphlan2))] <- unlist(lapply( colnames(metaphlan2)[grep("CRC_MR", colnames(metaphlan2))], function(x) unlist(strsplit(x, split = "_"))[[3]]))
common.ids <- intersect(rownames(metadata), colnames(metaphlan2))
metadata <- metadata[common.ids,]
metaphlan2 <- metaphlan2[,common.ids]

metaphlan2 <- as.data.frame(metaphlan2)
metaphlan2_2 <- apply(metaphlan2, 2, as.numeric)
rownames(metaphlan2_2) <- rownames(metaphlan2)
metaphlan2 <- metaphlan2_2
labelDescriptors <- data.frame(labelDescription=c(colnames(metadata)))
metadata_annotated <- AnnotatedDataFrame(data = metadata, varMetadata = labelDescriptors, dimLabels= c("sampleNames", "varLabels"))
metaphlan2 <- ExpressionSet(phenoData = metadata_annotated, assayData = metaphlan2)
metaphlan2 <- ExpressionSet2phyloseq(metaphlan2)

#Prune taxa
metaphlan_species <- subset_taxa(metaphlan2, !is.na(Species) & is.na(Strain))

table <- otu_table(metaphlan_species)
species_diversity <- data.frame(colnames(table))
samples <- colnames(table)
rownames(species_diversity) <- colnames(table)
species_diversity$Observed <- NA

for (sample in samples) {
  number_species <- length(which(table[,sample] > 0))
  species_diversity[sample,"Observed"] <- number_species
}

species_diversity$Group <- sample_data(metaphlan_species)$study_condition
species_diversity$Dataset <- sample_data(metaphlan_species)$dataset_name
species_diversity$Group <- revalue(species_diversity$Group, replace = c("control" = "Control","carcinoma" = "Carcinoma"))
species_diversity$Group <- factor(species_diversity$Group, levels = c("Control", "Carcinoma"))
species_diversity$Dataset <- revalue(species_diversity$Dataset, replace =c("ThomasAM_2019_c" = "YachidaS_2019"))
species_diversity$Dataset <- factor(species_diversity$Dataset, levels=c("YuJ_2015", "ThomasAM_2019_b", "FengQ_2015", "VogtmannE_2016", "ZellerG_2014", "YachidaS_2019", "WirbelJ_2018", "ThomasAM_2019_a", "GuptaA_2019"))

df_richness <- data.frame(Dataset=c(unique(as.character(species_diversity$Dataset))), pvalue = rep(NA, 9), y = c(rep(199, 9)))
rownames(df_richness) <- unique(as.character(species_diversity$Dataset))
for (i in unique(as.character(species_diversity$Dataset)) ) {
  df_richness[i,"pvalue"] <- round(wilcox.test(species_diversity$Observed[species_diversity$Group == "Carcinoma" & species_diversity$Dataset == i],
                                               species_diversity$Observed[species_diversity$Group == "Control" & species_diversity$Dataset == i])$p.value, digits = 8)
}
df_richness$Dataset <- factor(df_richness$Dataset, levels=c("YuJ_2015", "ThomasAM_2019_b", "FengQ_2015", "VogtmannE_2016", "ZellerG_2014", "YachidaS_2019", "WirbelJ_2018", "ThomasAM_2019_a", "GuptaA_2019"))
df_richness <- df_richness[c("YuJ_2015", "ThomasAM_2019_b", "FengQ_2015", "VogtmannE_2016", "ZellerG_2014", "YachidaS_2019", "WirbelJ_2018", "ThomasAM_2019_a", "GuptaA_2019"),]
df_richness$x <- c(1:9)
df_richness$pvalue <-formatC(df_richness$pvalue, format = "e", digits = 1)
df_richness$pvalue <- paste("P = ", df_richness$pvalue, sep = "")

p1 <- ggplot(species_diversity, aes(Dataset, Observed, fill = Group, colour = Group)) +
  geom_boxplot(aes(colour = Group), outlier.colour = NA, position = "dodge", fill="white") +
  geom_point(position=position_jitterdodge(dodge.width=0.8), alpha = 0.4, size = 0.5) +     
  theme_bw() +
  ylab("Number of species") +
  theme(axis.title.x= element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_text(size = 8), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 7, angle = 45, hjust = 1), legend.key.size = unit(0.3, "cm"), legend.position=c(0.25,0.2), legend.direction = "horizontal", legend.title = element_blank(), legend.text = element_text(size = 7), legend.background = element_rect(fill="transparent")) +
  scale_y_continuous(breaks = c(0, 50, 100, 150, 200), limits = c(0, 200), labels = c(0, 50, 100, 150, 200)) +
  scale_color_manual(labels = c("Control","Carcinoma"), values = c("#6F835A", "#aa483b")) +
  geom_text(data=df_richness, aes(x, y, label=pvalue), color = "black", size = 4, inherit.aes= F)
ggsave("figures/species_richness_mpa2_rarefied.svg", p1, width = 10, height = 8, device = "svg")

## ----------------- HumaNn3 pathways meta analysis ----------------------------------------------------------------------------------------------------------------------------------------
metadata <- read_tsv("tables/CRC_analysis_metadata_final_version.tsv", col_names = T) %>% 
  as.data.frame(row.names = FALSE, stringsAsFactors=FALSE)
metadata$study_condition <- gsub("CRC", "carcinoma", metadata$study_condition)
rownames(metadata)[metadata$dataset_name == "FengQ_2015"] <- as.character(metadata$sampleID[metadata$dataset_name == "FengQ_2015"])
rownames(metadata)[metadata$dataset_name == "ThomasAM_2019_c"] <- as.character(metadata$sampleID[metadata$dataset_name == "ThomasAM_2019_c"])
rownames(metadata)[metadata$dataset_name == "ZellerG_2014"] <- as.character(metadata$sampleID[metadata$dataset_name == "ZellerG_2014"])
rownames(metadata)[metadata$dataset_name == "YachidaS_2019"] <- as.character(metadata$sampleID[metadata$dataset_name == "YachidaS_2019"])
rownames(metadata)[metadata$dataset_name == "YuJ_2015"] <- as.character(metadata$sampleID[metadata$dataset_name == "YuJ_2015"])
rownames(metadata)[metadata$dataset_name == "ThomasAM_2019_a"] <- as.character(metadata$subjectID[metadata$dataset_name == "ThomasAM_2019_a"])
rownames(metadata)[metadata$dataset_name == "ThomasAM_2019_b"] <- as.character(metadata$subjectID[metadata$dataset_name == "ThomasAM_2019_b"])
rownames(metadata)[metadata$dataset_name == "VogtmannE_2016"] <- as.character(metadata$subjectID[metadata$dataset_name == "VogtmannE_2016"])
rownames(metadata)[metadata$dataset_name == "WirbelJ_2018"] <- as.character(metadata$subjectID[metadata$dataset_name == "WirbelJ_2018"])
rownames(metadata)[metadata$dataset_name == "GuptaA_2019"] <- as.character(metadata$sampleID[metadata$dataset_name == "GuptaA_2019"])

pathways <- read.table("tables/humann3_pathabundance_CRC_cohorts_cpm_unstratified.tsv", stringsAsFactors = F, header = F)
colnames(pathways) <- pathways[1,]
rownames(pathways) <- pathways[,1]
pathways <- pathways[-1,]
pathways <- pathways[,-1]
colnames(pathways) <- gsub("_Abundance", "", colnames(pathways))
colnames(pathways)[grep("LILT", colnames(pathways))] <- unlist(lapply( colnames(pathways)[grep("LILT", colnames(pathways))], function(x) unlist(strsplit(x, split = "_"))[[2]]))
colnames(pathways)[grep("CRC_MR", colnames(pathways))] <- unlist(lapply( colnames(pathways)[grep("CRC_MR", colnames(pathways))], function(x) unlist(strsplit(x, split = "_"))[[3]]))
common.ids <- intersect(colnames(pathways), rownames(metadata))

pathways <- as.data.frame(pathways)
pathways_2 <- apply(pathways, 2, as.numeric)
rownames(pathways_2) <- rownames(pathways)
pathways <- pathways_2
pathways <- apply(pathways, 2, function(x) asin(sqrt(x/sum(x))))
pathways <- pathways[,common.ids]
metadata <- metadata[common.ids,]

meta_analysis_results <- matrix(NA, ncol = 18, nrow=nrow(pathways))
colnames(meta_analysis_results) <- c("Species", "Number_Samples", "I^2", "P-value Fit", "Pvalue", "Standard Error", "Coefficient", "Confidence Interval Lower Limit", "Confidence Interval Upper Limit", "ZellerG_2014", "FengQ_2015", "YuJ_2015", "ThomasAM_2019_a", "ThomasAM_2019_b", "VogtmannE_2016","WirbelJ_2019", "YachidaS_2019", 'GuptaA_2019')
meta_analysis_results <- as.data.frame(meta_analysis_results)
species <- pathways
for (k in 1:nrow(pathways)) {
  species_table <- matrix(NA, ncol = 8, nrow = 10)
  colnames(species_table) <- c("Study", "Dataset", "SampleSize_Cancer", "Mean_Cancer", "SD_Cancer", "SampleSize_Control", "Mean_Control", "SD_Control")
  species_table <- as.data.frame(species_table)
  datasets <- as.character(unique(metadata$dataset_name))
  meta_analysis_results$VogtmannE_2016[k] <- wilcox.test(as.numeric(species[k,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "VogtmannE_2016"]]),
                                                         as.numeric(species[k,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "VogtmannE_2016"]]))$p.value
  meta_analysis_results$ThomasAM_2019_b[k] <- wilcox.test(as.numeric(species[k,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "ThomasAM_2019_b"]]),
                                                          as.numeric(species[k,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "ThomasAM_2019_b"]]))$p.value
  meta_analysis_results$ThomasAM_2019_a[k] <- wilcox.test(as.numeric(species[k,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "ThomasAM_2019_a"]]),
                                                          as.numeric(species[k,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "ThomasAM_2019_a"]]))$p.value
  meta_analysis_results$YuJ_2015[k] <- wilcox.test(as.numeric(species[k,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "YuJ_2015"]]),
                                                   as.numeric(species[k,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "YuJ_2015"]]))$p.value
  meta_analysis_results$FengQ_2015[k] <- wilcox.test(as.numeric(species[k,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "FengQ_2015"]]),
                                                     as.numeric(species[k,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "FengQ_2015"]]))$p.value
  meta_analysis_results$ZellerG_2014[k] <- wilcox.test(as.numeric(species[k,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "ZellerG_2014"]]),
                                                       as.numeric(species[k,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "ZellerG_2014"]]))$p.value
  meta_analysis_results$WirbelJ_2019[k] <- wilcox.test(as.numeric(species[k,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "WirbelJ_2018"]]),
                                                       as.numeric(species[k,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "WirbelJ_2018"]]))$p.value
  meta_analysis_results$YachidaS_2019[k] <- wilcox.test(as.numeric(species[k,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "YachidaS_2019"]]),
                                                        as.numeric(species[k,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "YachidaS_2019"]]))$p.value
  meta_analysis_results$GuptaA_2019[k] <- wilcox.test(as.numeric(species[k,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "GuptaA_2019"]]),
                                                      as.numeric(species[k,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "GuptaA_2019"]]))$p.value
  for (i in 1:length(datasets)) {
    dataset <- datasets[i]
    species_table$Study[i] <- i
    species_table$Dataset[i] <- dataset
    species_table$SampleSize_Cancer[i] <- length(which(metadata$study_condition == "carcinoma" & metadata$dataset_name == dataset))
    species_table$Mean_Cancer[i] <- mean(species[k,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == dataset]])
    species_table$SD_Cancer[i] <- sd(species[k,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == dataset]])
    species_table$SampleSize_Control[i] <- length(which(metadata$study_condition == "control" & metadata$dataset_name == dataset))
    species_table$Mean_Control[i] <- mean(species[k,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == dataset]])
    species_table$SD_Control[i] <- sd(species[k,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == dataset]])
  }
  dat1 <- escalc(measure="SMD", m1i=Mean_Cancer, sd1i=SD_Cancer, n1i=SampleSize_Cancer,
                 m2i=Mean_Control, sd2i=SD_Control, n2i=SampleSize_Control, data=species_table, vtype="UB")
  if (length(which((dat1$Mean_Cancer == dat1$Mean_Control) == TRUE)) > 3) {
    meta_analysis_results$Species[k] <- rownames(pathways)[k]
    meta_analysis_results$Number_Samples[k] <- 0  
    next()
  }
  else{
    res1 <- rma(yi, vi, data=dat1, control=list(stepadj=0.5, maxiter = 10000), verbose=TRUE, digits=5)
    meta_analysis_results$Species[k] <- rownames(pathways)[k]
    meta_analysis_results$Number_Samples[k] <- sum(res1$ni)
    meta_analysis_results$`I^2`[k] <- res1$I2
    meta_analysis_results$`P-value Fit`[k] <- res1$QEp
    meta_analysis_results$Pvalue[k] <- res1$pval
    meta_analysis_results$`Standard Error`[k] <- res1$b
    meta_analysis_results$Coefficient[k] <- res1$b
    meta_analysis_results$`Confidence Interval Lower Limit`[k] <- res1$ci.lb
    meta_analysis_results$`Confidence Interval Upper Limit`[k] <- res1$ci.ub
  }
}
meta_analysis_results$PvalueAdjusted <- p.adjust(meta_analysis_results$Pvalue, method = "fdr")
meta_analysis_results <- meta_analysis_results[order(abs(meta_analysis_results$Coefficient), decreasing = T), ]
meta_analysis_results <- meta_analysis_results[meta_analysis_results$PvalueAdjusted < 0.05,]
write.table(meta_analysis_results, "results/meta_analysis_results_humann3_arcsine_all_datasets.txt", col.names = T, row.names = F, quote = F, sep = "\t")

meta_analysis_results <- meta_analysis_results[!(is.na(meta_analysis_results$`I^2`)),]
meta_analysis_results <- meta_analysis_results[order(abs(meta_analysis_results$Coefficient), decreasing = T), ]
meta_analysis_results <- meta_analysis_results[meta_analysis_results$`P-value Fit` > 0.05,]
write.table(meta_analysis_results, "results/meta_analysis_results_humann3_arcsine_all_datasets_filt.txt", col.names = T, row.names = F, quote = F, sep = "\t")

meta_analysis_results_2 <- meta_analysis_results[1:20,]
write.table(meta_analysis_results_2, "results/meta_analysis_results_humann3_arcsine_all_datasets_top20.txt", col.names = T, row.names = F, quote = F, sep = "\t")

species_meta <- pathways[meta_analysis_results_2$Species,]
effect_sizes <- matrix(NA, ncol = 10, nrow=nrow(species_meta))
colnames(effect_sizes) <- c("Species", "ZellerG_2014", "FengQ_2015", "YuJ_2015", "ThomasAM_2019_b", "ThomasAM_2019_b", "VogtmannE_2016", "WirbelJ_2018", "YachidaS_2019", 'GuptaA_2019')
effect_sizes <- as.data.frame(effect_sizes)
for (i in 1:nrow(species_meta)) {
  effect_sizes$Species[i] <- rownames(species_meta)[i]
  effect_sizes$VogtmannE_2016[i] <- effect_size_calc(as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "VogtmannE_2016"]]), as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "VogtmannE_2016"]]))
  effect_sizes$ThomasAM_2019_a[i] <- effect_size_calc(as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "ThomasAM_2019_a"]]), as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "ThomasAM_2019_a"]]))
  effect_sizes$ThomasAM_2019_b[i] <- effect_size_calc(as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "ThomasAM_2019_b"]]), as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "ThomasAM_2019_b"]]))
  effect_sizes$YuJ_2015[i] <- effect_size_calc(as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "YuJ_2015"]]), as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "YuJ_2015"]]))
  effect_sizes$FengQ_2015[i] <- effect_size_calc(as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "FengQ_2015"]]), as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "FengQ_2015"]]))
  effect_sizes$ZellerG_2014[i] <- effect_size_calc(as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "ZellerG_2014"]]),as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "ZellerG_2014"]]))
  effect_sizes$YachidaS_2019[i] <- effect_size_calc(as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "YachidaS_2019"]]),as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "YachidaS_2019"]]))
  effect_sizes$GuptaA_2019[i] <- effect_size_calc(as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "GuptaA_2019"]]),as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "GuptaA_2019"]]))
  effect_sizes$WirbelJ_2018[i] <- effect_size_calc(as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "WirbelJ_2018"]]),as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "WirbelJ_2018"]]))
}
meta_analysis_results_2$Species <- gsub("_", " ", meta_analysis_results_2$Species)
effect_sizes$Species <- gsub("_", " ", effect_sizes$Species)
effect_sizes$Species <- factor(effect_sizes$Species, levels = rev(c(meta_analysis_results_2$Species[order(meta_analysis_results_2$Coefficient, decreasing = T)])))


effect_sizes[,-1] <- effect_sizes[,-1] * -1

effect_sizes$RandomEffects <- meta_analysis_results_2$Coefficient
effect_sizes <- melt(effect_sizes)
effect_sizes$CI_Low <- NA 
effect_sizes$CI_UP <- NA
effect_sizes[effect_sizes$variable == "RandomEffects","CI_Low"] <- meta_analysis_results_2$`Confidence Interval Lower Limit`
effect_sizes[effect_sizes$variable == "RandomEffects","CI_UP"] <- meta_analysis_results_2$`Confidence Interval Upper Limit`
effect_sizes$variable <- revalue(effect_sizes$variable, replace = c("WirbelJ_2018" = "WirbelJ_2019"))
effect_sizes$variable <- factor(effect_sizes$variable, levels = c("RandomEffects", "ThomasAM_2019_a", "ThomasAM_2019_b", "YachidaS_2019", "GuptaA_2019", "WirbelJ_2019", "HanniganGD_2017", "VogtmannE_2016", "FengQ_2015", "YuJ_2015", "ZellerG_2014"))
effect_sizes$variable <- revalue(effect_sizes$variable, replace = c("RandomEffects" = "Random Effects"))

p1 <- ggplot(effect_sizes, aes(Species, value, shape = variable, ymin=CI_Low, ymax=CI_UP)) + 
  geom_point(size = 2.5) +
  geom_linerange() + 
  scale_shape_manual(values= rev(c(3,4,8,0, 13, 11, 12, 2,5,19))) +
  coord_flip() +
  theme_minimal() +
  geom_hline(aes(yintercept = 0), color = "gray", linetype = "dashed") +
  theme(axis.text.y = element_text(face = "italic", size = 8), axis.title.y = element_blank(), axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 12), legend.title = element_blank(),  
        legend.key.size = unit(0.7, "cm"), legend.position=c(0.1,0.75), legend.text = element_text(size = 7.5)) + 
  ylab("Effect size") 
ggsave("figures/meta_analysis_pathways_top20.svg", p1, width = 10, height = 8, device = "svg")

## ----------------- HumaNn3 ECs meta analysis ----------------------------------------------------------------------------------------------------------------------------------------
metadata <- read_tsv("tables/CRC_analysis_metadata_final_version.tsv", col_names = T) %>% 
  as.data.frame(row.names = FALSE, stringsAsFactors=FALSE)
metadata$study_condition <- gsub("CRC", "carcinoma", metadata$study_condition)
rownames(metadata)[metadata$dataset_name == "FengQ_2015"] <- as.character(metadata$sampleID[metadata$dataset_name == "FengQ_2015"])
rownames(metadata)[metadata$dataset_name == "ThomasAM_2019_c"] <- as.character(metadata$sampleID[metadata$dataset_name == "ThomasAM_2019_c"])
rownames(metadata)[metadata$dataset_name == "ZellerG_2014"] <- as.character(metadata$sampleID[metadata$dataset_name == "ZellerG_2014"])
rownames(metadata)[metadata$dataset_name == "YachidaS_2019"] <- as.character(metadata$sampleID[metadata$dataset_name == "YachidaS_2019"])
rownames(metadata)[metadata$dataset_name == "YuJ_2015"] <- as.character(metadata$sampleID[metadata$dataset_name == "YuJ_2015"])
rownames(metadata)[metadata$dataset_name == "ThomasAM_2019_a"] <- as.character(metadata$subjectID[metadata$dataset_name == "ThomasAM_2019_a"])
rownames(metadata)[metadata$dataset_name == "ThomasAM_2019_b"] <- as.character(metadata$subjectID[metadata$dataset_name == "ThomasAM_2019_b"])
rownames(metadata)[metadata$dataset_name == "VogtmannE_2016"] <- as.character(metadata$subjectID[metadata$dataset_name == "VogtmannE_2016"])
rownames(metadata)[metadata$dataset_name == "WirbelJ_2018"] <- as.character(metadata$subjectID[metadata$dataset_name == "WirbelJ_2018"])
rownames(metadata)[metadata$dataset_name == "GuptaA_2019"] <- as.character(metadata$sampleID[metadata$dataset_name == "GuptaA_2019"])

pathways <- read.table("tables/CRC_humann3_ecs_cpm_ungrouped.tsv", stringsAsFactors = F, header = T, row.names = 1, check.names = F)
colnames(pathways) <- gsub("_Abundance-CPM", "", colnames(pathways))
colnames(pathways)[grep("LILT", colnames(pathways))] <- unlist(lapply( colnames(pathways)[grep("LILT", colnames(pathways))], function(x) unlist(strsplit(x, split = "_"))[[2]]))
colnames(pathways)[grep("CRC_MR", colnames(pathways))] <- unlist(lapply( colnames(pathways)[grep("CRC_MR", colnames(pathways))], function(x) unlist(strsplit(x, split = "_"))[[3]]))
common.ids <- intersect(colnames(pathways), rownames(metadata))

pathways <- as.data.frame(pathways)
pathways_2 <- apply(pathways, 2, as.numeric)
rownames(pathways_2) <- rownames(pathways)
pathways <- pathways_2
pathways <- apply(pathways, 2, function(x) asin(sqrt(x/sum(x))))
pathways <- pathways[,common.ids]
metadata <- metadata[common.ids,]

meta_analysis_results <- matrix(NA, ncol = 18, nrow=nrow(pathways))
colnames(meta_analysis_results) <- c("Species", "Number_Samples", "I^2", "P-value Fit", "Pvalue", "Standard Error", "Coefficient", "Confidence Interval Lower Limit", "Confidence Interval Upper Limit", "ZellerG_2014", "FengQ_2015", "YuJ_2015", "ThomasAM_2019_a", "ThomasAM_2019_b", "VogtmannE_2016","WirbelJ_2019", "YachidaS_2019", 'GuptaA_2019')
meta_analysis_results <- as.data.frame(meta_analysis_results)
species <- pathways
for (k in 1:nrow(pathways)) {
  species_table <- matrix(NA, ncol = 8, nrow = 10)
  colnames(species_table) <- c("Study", "Dataset", "SampleSize_Cancer", "Mean_Cancer", "SD_Cancer", "SampleSize_Control", "Mean_Control", "SD_Control")
  species_table <- as.data.frame(species_table)
  datasets <- as.character(unique(metadata$dataset_name))
  meta_analysis_results$VogtmannE_2016[k] <- wilcox.test(as.numeric(species[k,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "VogtmannE_2016"]]),
                                                         as.numeric(species[k,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "VogtmannE_2016"]]))$p.value
  meta_analysis_results$ThomasAM_2019_b[k] <- wilcox.test(as.numeric(species[k,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "ThomasAM_2019_b"]]),
                                                          as.numeric(species[k,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "ThomasAM_2019_b"]]))$p.value
  meta_analysis_results$ThomasAM_2019_a[k] <- wilcox.test(as.numeric(species[k,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "ThomasAM_2019_a"]]),
                                                          as.numeric(species[k,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "ThomasAM_2019_a"]]))$p.value
  meta_analysis_results$YuJ_2015[k] <- wilcox.test(as.numeric(species[k,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "YuJ_2015"]]),
                                                   as.numeric(species[k,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "YuJ_2015"]]))$p.value
  meta_analysis_results$FengQ_2015[k] <- wilcox.test(as.numeric(species[k,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "FengQ_2015"]]),
                                                     as.numeric(species[k,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "FengQ_2015"]]))$p.value
  meta_analysis_results$ZellerG_2014[k] <- wilcox.test(as.numeric(species[k,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "ZellerG_2014"]]),
                                                       as.numeric(species[k,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "ZellerG_2014"]]))$p.value
  meta_analysis_results$WirbelJ_2019[k] <- wilcox.test(as.numeric(species[k,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "WirbelJ_2018"]]),
                                                       as.numeric(species[k,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "WirbelJ_2018"]]))$p.value
  meta_analysis_results$YachidaS_2019[k] <- wilcox.test(as.numeric(species[k,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "YachidaS_2019"]]),
                                                        as.numeric(species[k,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "YachidaS_2019"]]))$p.value
  meta_analysis_results$GuptaA_2019[k] <- wilcox.test(as.numeric(species[k,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "GuptaA_2019"]]),
                                                      as.numeric(species[k,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "GuptaA_2019"]]))$p.value
  for (i in 1:length(datasets)) {
    dataset <- datasets[i]
    species_table$Study[i] <- i
    species_table$Dataset[i] <- dataset
    species_table$SampleSize_Cancer[i] <- length(which(metadata$study_condition == "carcinoma" & metadata$dataset_name == dataset))
    species_table$Mean_Cancer[i] <- mean(species[k,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == dataset]])
    species_table$SD_Cancer[i] <- sd(species[k,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == dataset]])
    species_table$SampleSize_Control[i] <- length(which(metadata$study_condition == "control" & metadata$dataset_name == dataset))
    species_table$Mean_Control[i] <- mean(species[k,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == dataset]])
    species_table$SD_Control[i] <- sd(species[k,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == dataset]])
  }
  dat1 <- escalc(measure="SMD", m1i=Mean_Cancer, sd1i=SD_Cancer, n1i=SampleSize_Cancer,
                 m2i=Mean_Control, sd2i=SD_Control, n2i=SampleSize_Control, data=species_table, vtype="UB")
  if (length(which((dat1$Mean_Cancer == dat1$Mean_Control) == TRUE)) > 3) {
    meta_analysis_results$Species[k] <- rownames(pathways)[k]
    meta_analysis_results$Number_Samples[k] <- 0  
    next()
  }
  else{
    res1 <- rma(yi, vi, data=dat1, control=list(stepadj=0.5, maxiter = 10000), verbose=TRUE, digits=5)
    meta_analysis_results$Species[k] <- rownames(pathways)[k]
    meta_analysis_results$Number_Samples[k] <- sum(res1$ni)
    meta_analysis_results$`I^2`[k] <- res1$I2
    meta_analysis_results$`P-value Fit`[k] <- res1$QEp
    meta_analysis_results$Pvalue[k] <- res1$pval
    meta_analysis_results$`Standard Error`[k] <- res1$b
    meta_analysis_results$Coefficient[k] <- res1$b
    meta_analysis_results$`Confidence Interval Lower Limit`[k] <- res1$ci.lb
    meta_analysis_results$`Confidence Interval Upper Limit`[k] <- res1$ci.ub
  }
}
meta_analysis_results$PvalueAdjusted <- p.adjust(meta_analysis_results$Pvalue, method = "fdr")
meta_analysis_results <- meta_analysis_results[order(abs(meta_analysis_results$Coefficient), decreasing = T), ]
meta_analysis_results <- meta_analysis_results[meta_analysis_results$PvalueAdjusted < 0.05,]
write.table(meta_analysis_results, "results/meta_analysis_results_humann3_ecs_arcsine_all_datasets.txt", col.names = T, row.names = F, quote = F, sep = "\t")

meta_analysis_results <- meta_analysis_results[!(is.na(meta_analysis_results$`I^2`)),]
meta_analysis_results <- meta_analysis_results[order(abs(meta_analysis_results$Coefficient), decreasing = T), ]
meta_analysis_results <- meta_analysis_results[meta_analysis_results$`P-value Fit` > 0.05,]
write.table(meta_analysis_results, "results/meta_analysis_results_humann3_ecs_arcsine_all_datasets_filt.txt", col.names = T, row.names = F, quote = F, sep = "\t")

pwy_names <- read.delim("tables/map_ec_name.txt", header = F, sep = "\t", as.is = T, check.names = F, stringsAsFactors = F)
meta_analysis_results$Description <- NA
for (i in 1:nrow(meta_analysis_results)) {
  new_name <- pwy_names$V2[grep(paste("^", as.character(meta_analysis_results$Species[i]), "$", sep = ""), pwy_names$V1)]
  meta_analysis_results$Description[i] <- new_name 
}

meta_analysis_results_2 <- meta_analysis_results[1:20,]
write.table(meta_analysis_results_2, "results/meta_analysis_results_humann3_ecs_arcsine_all_datasets_top20.txt", col.names = T, row.names = F, quote = F, sep = "\t")

species_meta <- pathways[meta_analysis_results_2$Species,]
effect_sizes <- matrix(NA, ncol = 10, nrow=nrow(species_meta))
colnames(effect_sizes) <- c("Species", "ZellerG_2014", "FengQ_2015", "YuJ_2015", "ThomasAM_2019_a", "ThomasAM_2019_b", "VogtmannE_2016", "WirbelJ_2018", "YachidaS_2019", 'GuptaA_2019')
effect_sizes <- as.data.frame(effect_sizes)
for (i in 1:nrow(species_meta)) {
  effect_sizes$Species[i] <- rownames(species_meta)[i]
  effect_sizes$VogtmannE_2016[i] <- effect_size_calc(as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "VogtmannE_2016"]]), as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "VogtmannE_2016"]]))
  effect_sizes$ThomasAM_2019_a[i] <- effect_size_calc(as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "ThomasAM_2019_a"]]), as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "ThomasAM_2019_a"]]))
  effect_sizes$ThomasAM_2019_b[i] <- effect_size_calc(as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "ThomasAM_2019_b"]]), as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "ThomasAM_2019_b"]]))
  effect_sizes$YuJ_2015[i] <- effect_size_calc(as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "YuJ_2015"]]), as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "YuJ_2015"]]))
  effect_sizes$FengQ_2015[i] <- effect_size_calc(as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "FengQ_2015"]]), as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "FengQ_2015"]]))
  effect_sizes$ZellerG_2014[i] <- effect_size_calc(as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "ZellerG_2014"]]),as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "ZellerG_2014"]]))
  effect_sizes$YachidaS_2019[i] <- effect_size_calc(as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "YachidaS_2019"]]),as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "YachidaS_2019"]]))
  effect_sizes$GuptaA_2019[i] <- effect_size_calc(as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "GuptaA_2019"]]),as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "GuptaA_2019"]]))
  effect_sizes$WirbelJ_2018[i] <- effect_size_calc(as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "carcinoma" & metadata$dataset_name == "WirbelJ_2018"]]),as.numeric(species_meta[i,rownames(metadata)[metadata$study_condition == "control" & metadata$dataset_name == "WirbelJ_2018"]]))
}

effect_sizes$Description <- meta_analysis_results_2$Description
effect_sizes$Description <- factor(effect_sizes$Description, levels = effect_sizes$Description[order(meta_analysis_results_2$Coefficient)]) 

effect_sizes[,c(-1,-11)] <- effect_sizes[,c(-1,-11)] * -1

effect_sizes$RandomEffects <- meta_analysis_results_2$Coefficient
effect_sizes <- melt(effect_sizes)
effect_sizes$CI_Low <- NA 
effect_sizes$CI_UP <- NA
effect_sizes[effect_sizes$variable == "RandomEffects","CI_Low"] <- meta_analysis_results_2$`Confidence Interval Lower Limit`
effect_sizes[effect_sizes$variable == "RandomEffects","CI_UP"] <- meta_analysis_results_2$`Confidence Interval Upper Limit`
effect_sizes$variable <- revalue(effect_sizes$variable, replace = c("WirbelJ_2018" = "WirbelJ_2019"))
effect_sizes$variable <- factor(effect_sizes$variable, levels = c("RandomEffects", "ThomasAM_2019_a", "ThomasAM_2019_b", "YachidaS_2019", "GuptaA_2019", "WirbelJ_2019", "HanniganGD_2017", "VogtmannE_2016", "FengQ_2015", "YuJ_2015", "ZellerG_2014"))
effect_sizes$variable <- revalue(effect_sizes$variable, replace = c("RandomEffects" = "Random Effects"))

p1 <- ggplot(effect_sizes, aes(Description, value, shape = variable, ymin=CI_Low, ymax=CI_UP)) + 
  geom_point(size = 2.5) +
  geom_linerange() + 
  scale_shape_manual(values= rev(c(3,4,8,0, 13, 11, 12, 2,5,19))) +
  coord_flip() +
  theme_minimal() +
  geom_hline(aes(yintercept = 0), color = "gray", linetype = "dashed") +
  theme(axis.text.y = element_text(face = "italic", size = 8), axis.title.y = element_blank(), axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 12), legend.title = element_blank(),  
        legend.key.size = unit(0.7, "cm"), legend.position=c(0.1,0.75), legend.text = element_text(size = 7.5)) + 
  ylab("Effect size") 
ggsave("figures/meta_analysis_ecs_top20.svg", p1, width = 10, height = 8, device = "svg")

## -----------------humann3 cutC and yeaW ----------------------------------------------------------------------------------------------------------------------------------------
metadata <- read_tsv("tables/CRC_analysis_metadata_final_version.tsv", col_names = T) %>% 
  as.data.frame(row.names = FALSE, stringsAsFactors=FALSE)
metadata$study_condition <- gsub("CRC", "carcinoma", metadata$study_condition)
rownames(metadata)[metadata$dataset_name == "FengQ_2015"] <- as.character(metadata$sampleID[metadata$dataset_name == "FengQ_2015"])
rownames(metadata)[metadata$dataset_name == "ThomasAM_2019_c"] <- as.character(metadata$sampleID[metadata$dataset_name == "ThomasAM_2019_c"])
rownames(metadata)[metadata$dataset_name == "ZellerG_2014"] <- as.character(metadata$sampleID[metadata$dataset_name == "ZellerG_2014"])
rownames(metadata)[metadata$dataset_name == "YachidaS_2019"] <- as.character(metadata$sampleID[metadata$dataset_name == "YachidaS_2019"])
rownames(metadata)[metadata$dataset_name == "YuJ_2015"] <- as.character(metadata$sampleID[metadata$dataset_name == "YuJ_2015"])
rownames(metadata)[metadata$dataset_name == "ThomasAM_2019_a"] <- as.character(metadata$subjectID[metadata$dataset_name == "ThomasAM_2019_a"])
rownames(metadata)[metadata$dataset_name == "ThomasAM_2019_b"] <- as.character(metadata$subjectID[metadata$dataset_name == "ThomasAM_2019_b"])
rownames(metadata)[metadata$dataset_name == "VogtmannE_2016"] <- as.character(metadata$subjectID[metadata$dataset_name == "VogtmannE_2016"])
rownames(metadata)[metadata$dataset_name == "WirbelJ_2018"] <- as.character(metadata$subjectID[metadata$dataset_name == "WirbelJ_2018"])
rownames(metadata)[metadata$dataset_name == "GuptaA_2019"] <- as.character(metadata$sampleID[metadata$dataset_name == "GuptaA_2019"])

yeaW_NR90_ids <- strsplit(read_file("yeaW_NR90"), '\n')[[1]]
NR90_yeaW <- read_tsv("merged_yeaW_RPK_CPM.tsv", col_names = c('NR90', 'RPK', 'sample')) %>% 
  filter(NR90 %in% yeaW_NR90_ids) %>% 
  mutate(sample = str_remove(sample, "_genefamilies.*")) %>% 
  filter(sample %in% metadata$sampleID) %>% 
  pivot_wider(names_from = 'sample', values_from = 'RPK', values_fill = list(RPK = 0)) %>% 
  column_to_rownames('NR90')
missing_samples_lst <- metadata[!(metadata$sampleID %in% colnames(NR90_yeaW)), 2]
missing_samples <- as.data.frame(matrix(data = 0, nrow = nrow(NR90_yeaW), ncol = length(missing_samples_lst)))
colnames(missing_samples) <- missing_samples_lst
rownames(missing_samples) <- rownames(NR90_yeaW)
NR90_yeaW <- cbind(missing_samples, NR90_yeaW)

richness_meta_analysis_results <- matrix(NA, ncol = 9, nrow=1)
colnames(richness_meta_analysis_results) <- c("Metric", "Number_Samples", "I^2", "P-value Fit", "Pvalue", "Standard Error", "Coefficient", "Confidence Interval Lower Limit", "Confidence Interval Upper Limit")
richness_meta_analysis_results <- as.data.frame(richness_meta_analysis_results)
richness_table <- matrix(NA, ncol = 8, nrow = 9)
colnames(richness_table) <- c("Study", "Dataset", "SampleSize_Cancer", "Mean_Cancer", "SD_Cancer", "SampleSize_Control", "Mean_Control", "SD_Control")
richness_table <- as.data.frame(richness_table)
richness_table$Study <- 1:9
metadata$dataset_name <- gsub("ThomasAM_2019_c", "YachidaS_2019", metadata$dataset_name)
datasets <- unique(metadata$dataset_name)

for (i in 1:length(datasets)) {
  dataset <- datasets[i]
  richness_table$Dataset[i] <- dataset
  richness_table$SampleSize_Cancer[i] <- length(colSums(NR90_yeaW[,metadata$sampleID[metadata$dataset_name == dataset & metadata$study_condition == "carcinoma"]]))
  richness_table$Mean_Cancer[i] <- mean(colSums(NR90_yeaW[,metadata$sampleID[metadata$dataset_name == dataset & metadata$study_condition == "carcinoma"]]))
  richness_table$SD_Cancer[i] <- sd(colSums(NR90_yeaW[,metadata$sampleID[metadata$dataset_name == dataset & metadata$study_condition == "carcinoma"]]))
  richness_table$SampleSize_Control[i] <- length(colSums(NR90_yeaW[,metadata$sampleID[metadata$dataset_name == dataset & metadata$study_condition == "control"]]))
  richness_table$Mean_Control[i] <- mean(colSums(NR90_yeaW[,metadata$sampleID[metadata$dataset_name == dataset & metadata$study_condition == "control"]]))
  richness_table$SD_Control[i] <- sd(colSums(NR90_yeaW[,metadata$sampleID[metadata$dataset_name == dataset & metadata$study_condition == "control"]]))
}
dat1 <- escalc(measure="SMD", m1i=Mean_Cancer, sd1i=SD_Cancer, n1i=SampleSize_Cancer,
               m2i=Mean_Control, sd2i=SD_Control, n2i=SampleSize_Control, data=richness_table)
res1 <- rma(yi, vi, data=dat1)
res1
m.hksj.raw <- metacont(n.e = SampleSize_Cancer, mean.e = Mean_Cancer, 
                       sd.e = SD_Cancer, n.c = SampleSize_Control, mean.c = Mean_Control, 
                       sd.c = SD_Control, data = dat1,
                       studlab = paste(Dataset),
                       comb.fixed = FALSE,
                       comb.random = TRUE,
                       prediction = TRUE, sm = "SMD")
forest(m.hksj.raw, sortvar = m.hksj.raw$TE, lab.e = "CRC", comb.random = T, print.tau2 = F)
forest(res1, slab=dat1$Dataset, col = "dark blue", annotate = T, order = order(res1$yi))
text(-3, -1.6, pos=4, "P-value = 0.07")
text(-0.5, 10.5, pos=4, "HUMAnN v3 yeaW UniRef90")

cutC_NR90_ids <- strsplit(read_file("cutC_NR90"), '\n')[[1]]
NR90_cutC <- read_tsv("merged_cutC_RPK_CPM.tsv", col_names = c('NR90', 'RPK', 'sample')) %>% 
  filter(NR90 %in% cutC_NR90_ids) %>% 
  mutate(sample = str_remove(sample, "_genefamilies.*")) %>% 
  filter(sample %in% metadata$sampleID) %>% 
  pivot_wider(names_from = 'sample', values_from = 'RPK', values_fill = list(RPK = 0)) %>% 
  column_to_rownames('NR90')

missing_samples_lst <- metadata[!(metadata$sampleID %in% colnames(NR90_cutC)),2]
missing_samples <- as.data.frame(matrix(data = 0, nrow = nrow(NR90_cutC), ncol = length(missing_samples_lst)))
colnames(missing_samples) <- missing_samples_lst
rownames(missing_samples) <- rownames(NR90_cutC)

NR90_cutC <- cbind(missing_samples, NR90_cutC)

richness_meta_analysis_results <- matrix(NA, ncol = 9, nrow=1)
colnames(richness_meta_analysis_results) <- c("Metric", "Number_Samples", "I^2", "P-value Fit", "Pvalue", "Standard Error", "Coefficient", "Confidence Interval Lower Limit", "Confidence Interval Upper Limit")
richness_meta_analysis_results <- as.data.frame(richness_meta_analysis_results)
richness_table <- matrix(NA, ncol = 8, nrow = 9)
colnames(richness_table) <- c("Study", "Dataset", "SampleSize_Cancer", "Mean_Cancer", "SD_Cancer", "SampleSize_Control", "Mean_Control", "SD_Control")
richness_table <- as.data.frame(richness_table)
richness_table$Study <- 1:9
datasets <- unique(metadata$dataset_name)

for (i in 1:length(datasets)) {
  dataset <- datasets[i]
  richness_table$Dataset[i] <- dataset
  richness_table$SampleSize_Cancer[i] <- length(colSums(NR90_cutC[,metadata$sampleID[metadata$dataset_name == dataset & metadata$study_condition == "carcinoma"]]))
  richness_table$Mean_Cancer[i] <- mean(colSums(NR90_cutC[,metadata$sampleID[metadata$dataset_name == dataset & metadata$study_condition == "carcinoma"]]))
  richness_table$SD_Cancer[i] <- sd(colSums(NR90_cutC[,metadata$sampleID[metadata$dataset_name == dataset & metadata$study_condition == "carcinoma"]]))
  richness_table$SampleSize_Control[i] <- length(colSums(NR90_cutC[,metadata$sampleID[metadata$dataset_name == dataset & metadata$study_condition == "control"]]))
  richness_table$Mean_Control[i] <- mean(colSums(NR90_cutC[,metadata$sampleID[metadata$dataset_name == dataset & metadata$study_condition == "control"]]))
  richness_table$SD_Control[i] <- sd(colSums(NR90_cutC[,metadata$sampleID[metadata$dataset_name == dataset & metadata$study_condition == "control"]]))
}
dat1 <- escalc(measure="SMD", m1i=Mean_Cancer, sd1i=SD_Cancer, n1i=SampleSize_Cancer,
               m2i=Mean_Control, sd2i=SD_Control, n2i=SampleSize_Control, data=richness_table)
res1 <- rma(yi, vi, data=dat1)
res1
m.hksj.raw <- metacont(n.e = SampleSize_Cancer, mean.e = Mean_Cancer, 
                       sd.e = SD_Cancer, n.c = SampleSize_Control, mean.c = Mean_Control, 
                       sd.c = SD_Control, data = dat1,
                       studlab = paste(Dataset),
                       comb.fixed = FALSE,
                       comb.random = TRUE,
                       prediction = TRUE, sm = "SMD")
forest(m.hksj.raw, sortvar = m.hksj.raw$TE, lab.e = "CRC", comb.random = T, text.fixed = "p = 0.02", print.tau2 = F)
forest(res1, slab=dat1$Dataset, col = "dark blue", annotate = T, order = order(res1$yi))
text(-2.45, -1.6, pos=4, "P-value <.0001")
text(-0.5, 10.5, pos=4, "HUMAnN v3 cutC UniRef90")