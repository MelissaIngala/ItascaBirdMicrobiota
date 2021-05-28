#install.packages('phyloseq')
library(phyloseq)
library(ggplot2)
library(biomformat)
library(wesanderson)
library(ape)
library(lme4)
library(glmmTMB)
library(cowplot)

########################## DATA IMPORT & PREPROCESSING ##########################

#Set working directory
setwd("~/Dropbox/Itasca_Feces2018/R Analyses/ItascaFeces")

#Import the BIOM file, tree, mapping file
dat <- read_biom(biom_file = "table-with-taxonomy.biom")
otu_table <- as.data.frame(as.matrix(biom_data(dat)))
taxonomy <- observation_metadata(dat)
drop <- "confidence"
taxonomy<-taxonomy[ , !(names(taxonomy) %in% drop)]
metadata <- read.delim(file = "~/Dropbox/Itasca_Feces2018/QIIME2 Files/Itasca_Metadata.txt",
                       header = T, row.names = 1)

phy<-read.tree(file = "Itasca_unrooted_16Stree.nwk")

#Import to phyloseq obj
SAM<-sample_data(metadata)
TAX<-tax_table(as.matrix(taxonomy))
OTU<- otu_table(otu_table, taxa_are_rows=TRUE)
PHY<-phy_tree(phy)
physeq<-merge_phyloseq(OTU, TAX, SAM, PHY)
colnames(tax_table(physeq))<-c("kingdom", "phylum", "class", "order", "family", "genus", "species")

#Remove samples with fewer than 1000 total reads & unassinged phyla
physeq.sub = prune_samples(sample_sums(physeq)>=1000, physeq) #lost 3 samples total
physeq.sub = subset_taxa(physeq.sub, phylum!="p__")
physeq.sub = subset_taxa(physeq.sub, phylum!="")
physeq.sub = prune_samples(physeq@sam_data$Sample_or_control!="Control", physeq)

#Calculate Relative Abundance
relative  = transform_sample_counts(physeq.sub, function(OTU) OTU / sum(OTU))

#Visually inspect barplot of top 250 OTUs
Top50OTUs = names(sort(taxa_sums(relative), TRUE)[1:550])
comparetop50 = prune_taxa(Top50OTUs, relative)
plot_bar(comparetop50, fill = "phylum", title = "Bacterial Phylum by Sample")

#################### DECONTAM: FILTER CONTAMS USING NEG CONTROL #######################
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("decontam")
library(decontam)

#Inspect library sizes
df <- as.data.frame(sample_data(physeq)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(physeq)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_control)) + geom_point()

# Use prevalence method to filter out suspected contaminants
sample_data(physeq)$is.neg <- sample_data(physeq)$Sample_or_control == "Control"
contamdf.prev <- isContaminant(physeq, method="prevalence", neg="is.neg", threshold = 0.5)
table(contamdf.prev$contaminant)

head(which(contamdf.prev$contaminant))

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(physeq, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_control == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_control == "Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

#Remove contams from phyloseq obj
ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, physeq)
ps.noncontam

##################### ALPHA DIVERSITY ####################################
#For this, do NOT use transformed dataset
#Subset out only nestlings (n = 107)
require(cowplot)
require(wesanderson)

nestlings<-subset_samples(physeq, stage=="Nestling")
richness<-estimate_richness(physeq = nestlings, split = T, measures = c("Chao1","Shannon"))
richness$host_species<-nestlings@sam_data$host_species
Foxy<- wes_palette("FantasticFox1", 3, type = "continuous")
shannon.nest.plot<-ggplot(richness, aes(host_species, Shannon, fill= host_species))+
  stat_boxplot(geom='errorbar', linetype=1, width=0.5) +  #whiskers
  geom_boxplot(outlier.shape=1, outlier.colour = "red") + scale_fill_manual(values = Foxy[2:3]) +
  theme_classic() + ggtitle("Nestlings")

#Subset out only adults (n = 29)
adults<-subset_samples(physeq, stage=="Adult")
richness2<-estimate_richness(physeq = adults, split = T, measures = c("Chao1","Shannon"))
richness2$host_species<-adults@sam_data$host_species
shannon.adults.plot<-ggplot(richness2, aes(host_species, Shannon, fill= host_species))+
  stat_boxplot(geom='errorbar', linetype=1, width=0.5) +  #whiskers
  geom_boxplot(outlier.shape=1, outlier.colour = "red") + scale_fill_manual(values = Foxy[2:3]) +
  scale_y_continuous(limits=c(0,6)) +
  theme_classic() + ggtitle("Adults")

###########TAXONOMIC BARPLOT FOR ALL BIRDS#####################

install.packages("remotes")
remotes::install_github("gmteunisse/Fantaxtic")
library(fantaxtic)

phyglom<-tax_glom(physeq = ps.noncontam, taxrank = "phylum")
phyglom<-subset_samples(physeq = phyglom, Sample_or_control!="Control")
tmp<-get_top_taxa(phyglom, 7, relative = TRUE, discard_other = FALSE,
             other_label = "Other")
ps_tmp <- name_taxa(tmp, label = "Unkown", species = T, other_label = "Other")

Wes<-wes_palette("Zissou1", 7, type = "continuous")
#pdf(file = "Tax_barplot.pdf", width = 7, height = 5)
fantaxtic_bar(ps_tmp, color_by = "phylum", label_by = "phylum", facet_by = "host_species",other_label = "Other", palette = Wes)
#dev.off()

###################### SCALE DATA TO ACCOUNT FOR DIFF LIBRARY SIZES #############
#Apply Hellinger transformation to scale library
#The Hellinger transform is square root of the relative abundance but instead given at the scale [0,1].
require(microbiome)
physeq.trans<-transform(x = ps.noncontam, transform = "hellinger", target = "OTU", shift = 0, scale = 1)

#Subset EABL apart from TRES
bluebirds<-subset_samples(nestlings, host_species=="EABL") #n=55
treeswallows<-subset_samples(nestlings, host_species=="TRES") # n=52

bluebirds.trans<-transform(x = bluebirds, transform = "hellinger", target = "OTU", shift = 0, scale = 1)
treeswallows.trans<-transform(x = treeswallows, transform = "hellinger", target = "OTU", shift = 0, scale = 1)

###################### BETA DIVERSITY: ORDINATION & PERMANOVA #########################
#Set seed for reproducibility
set.seed(42)

#Ordination using weighted and unweighted Unifrace dist (Lozupone et al. 2009)

#Allbirds


############################ BLUEBIRDS ##################################

# make a data frame from the sample_data
EABLsampledf <- data.frame(sample_data(bluebirds.trans))

#PERMANOVA of Species
set.seed(1)
# Calculate distance matrix
EABL_unifrac <- phyloseq::distance(bluebirds.trans, method = "unifrac")
EABL_wunifrac<- phyloseq::distance(bluebirds.trans, method = "wunifrac")
EABL_richness<-estimate_richness(bluebirds)
require(vegan)
# PERMANOVA test- composition
UNIFRAC.EABL<-adonis(EABL_unifrac ~ ParasiteTrtmt + HeatTreatmt + ParasiteTrtmt*HeatTreatmt, data = EABLsampledf)
UNIFRAC.EABL

WUNIFRAC.EABL<-adonis(EABL_wunifrac ~ ParasiteTrtmt + HeatTreatmt + ParasiteTrtmt*HeatTreatmt, data = EABLsampledf)
WUNIFRAC.EABL

p.adjust(p = c(UNIFRAC.EABL$aov.tab$`Pr(>F)`[1:3],WUNIFRAC.EABL$aov.tab$`Pr(>F)`[1:3] ))

########################### TREE SWALLOWS ##############################
# make a data frame from the sample_data
TRESsampledf <- data.frame(sample_data(treeswallows.trans))
#PERMANOVA of Species
set.seed(2)
# Calculate distance matrix
TRES_unifrac <- phyloseq::distance(treeswallows.trans, method = "unifrac")
TRES_wunifrac<- phyloseq::distance(treeswallows.trans, method = "wunifrac")
TRES_richness<-estimate_richness(treeswallows)

UNIFRAC.TRES<-adonis(TRES_unifrac ~ ParasiteTrtmt + HeatTreatmt + ParasiteTrtmt*HeatTreatmt, data = TRESsampledf)
UNIFRAC.TRES

WUNIFRAC.TRES<-adonis(TRES_wunifrac ~ ParasiteTrtmt + HeatTreatmt + ParasiteTrtmt*HeatTreatmt, data = TRESsampledf)
WUNIFRAC.TRES

p.adjust(p = c(UNIFRAC.TRES$aov.tab$`Pr(>F)`[1:3],WUNIFRAC.TRES$aov.tab$`Pr(>F)`[1:3] ))

#ALL BIRDS (supplemental info)
Itasca.ord.unifrac <- ordinate(
  physeq = physeq.trans, 
  method = "PCoA", 
  distance = "unifrac"
)

Itasca.ord.wunifrac <- ordinate(
  physeq = physeq.trans, 
  method = "PCoA", 
  distance = "wunifrac"
)

#PLOT
require(ggplot2)
require(wesanderson)

# Calculate distance matrix
spp_unifrac <- phyloseq::distance(physeq.trans, method = "unifrac")
spp_wunifrac<- phyloseq::distance(physeq.trans, method = "wunifrac")
study.bray<- phyloseq::distance(physeq.trans, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(physeq.trans))

Foxy<- wes_palette("FantasticFox1", 3, type = "continuous")
#pdf(file = "Itasca_unifrac_PCoA.pdf", width = 6, height = 5)
membership<-plot_ordination(
  physeq = physeq.trans,
  ordination = Itasca.ord.unifrac,
  axes = c(1,2), 
  color = "host_species",
  shape = "stage", title = "Membership") + scale_color_manual(values = Foxy) +
  theme_classic() + geom_point(aes(color = host_species, shape = stage), size = 2.0) +
  geom_jitter()
membership
#dev.off()

#pdf(file = "Itasca_wunifrac_PCoA.pdf", width = 6, height = 5)
composition<-plot_ordination(
  physeq = physeq.trans,
  ordination = Itasca.ord.wunifrac,
  axes = c(1,2), 
  color = "host_species",
  shape = "stage", title = "Composition") + scale_color_manual(values = Foxy) +
  theme_classic() + geom_point(aes(color = host_species, shape = stage), size = 2.0) +
  geom_jitter()
composition
#dev.off()

###### Plot alpha and beta diversity together
#pdf(file = "Alpha_beta_figs.pdf", width = 9, height = 6)
cowplot::plot_grid(shannon.nest.plot, shannon.adults.plot, membership, composition, labels = "AUTO")
#dev.off()

#PERMANOVA of Species
set.seed(24)

require(vegan)
# PERMANOVA test- composition
compspp<-adonis(spp_unifrac ~ host_species + ParasiteTrtmt + HeatTreatmt + stage, data = sampledf, strata = sampledf$nestID)
compspp

#PERMANOVA test- structure
structspp<-adonis(spp_wunifrac ~ host_species + ParasiteTrtmt + HeatTreatmt + stage, data = sampledf, strata = sampledf$nestID)
structspp

brayspp<-adonis(study.bray ~ host_species + ParasiteTrtmt + HeatTreatmt + stage, data = sampledf, strata = sampledf$nestID)
brayspp

#Adjust p value for multiple comparisons (Benjamini Hochberg)
spp.ps<-c(compspp$aov.tab$`Pr(>F)`[1:4], structspp$aov.tab$`Pr(>F)`[1:4], brayspp$aov.tab$`Pr(>F)`[1:4])
new.spp.ps<-p.adjust(spp.ps, method = "BH", n = 12)
new.spp.ps
#All P values corr to 0.004

## Homogeneity of dispersion test
beta <- betadisper(spp_unifrac, sampledf$host_species)
permutest(beta)
beta2 <- betadisper(spp_wunifrac, sampledf$host_species)
permutest(beta2)
#Differences may be due to differences in dispersion

##################### GLMMs: NEST PARASITISM & HEAT #########################

#FIXED EFFECTS: Parasite treatment: Parasitized (W) or Non-parasitized (P)
#Heat treatment: Heat (H) or No heat (NH)
#RANDOM EFFECT: Nest ID

#install.packages("glmmTMB")
library("glmmTMB")
#install.packages("bbmle")
library("bbmle") ## for AICtab
library("ggplot2")
## cosmetic
theme_set(theme_classic()+
            theme(panel.spacing=grid::unit(0,"lines")))

# Run model on nestlings (n = 107) only first
# Pre processing, prune OTUs not detected in at least 10 individual birds, and that do not have at least 50 total reads
nestlings.pruned = prune_taxa(taxa_sums(nestlings) > 10, nestlings) #reduced from 8108 to 3487 ASVs
nestlings.pruned = prune_taxa(taxa_sums(nestlings.pruned) > 100, nestlings.pruned) # further reduced to 1484 ASVs

#agglomerate to phylum and family level for GLMM models
nestlings.pruned.phy<- tax_glom(physeq = nestlings.pruned, taxrank = "phylum") 
length(unique(nestlings.pruned.phy@tax_table[,2])) # 16 phyla

nestlings.pruned.fam<- tax_glom(physeq = nestlings.pruned, taxrank = "family") 
length(unique(nestlings.pruned.fam@tax_table[,5])) #112 families

############################ GENERAL LINEAR MIXED MODELS #######################

######### "total is a df of all the sample data and phylum level OTU counts
#Split into EABL and TRES dfs
TRES<- subset(total, host_species=="TRES")
TRES$ParasiteTrtmt<-as.factor(TRES$ParasiteTrtmt)
TRES$HeatTreatmt<-as.factor(TRES$HeatTreatmt)


#LOOP THROUGH ALL ASVS  
models <- list()
dvnames <- noquote(paste0('`', colnames(TRES[5:20]), '`'))
ivnames <- paste("ParasiteTrtmt + HeatTreatmt + ParasiteTrtmt*HeatTreatmt + (1|nestID)", sep = ",")
library(glmmTMB)
for (y in dvnames){
  form <- formula(paste(y, "~", ivnames))
  models[[y]] <- glmmTMB(form, data= TRES, family = nbinom2(link = "log"))
}

TRES.summarylist<-lapply(models, summary) ## summarize each model

##################NOW EABL

EABL<- subset(total, host_species=="EABL")
EABL$ParasiteTrtmt<-as.factor(EABL$ParasiteTrtmt)
EABL$HeatTreatmt<-as.factor(EABL$HeatTreatmt)

models2 <- list()
dvnames2 <- noquote(paste0('`', colnames(EABL[5:20]), '`'))
ivnames2 <- paste("ParasiteTrtmt + HeatTreatmt + ParasiteTrtmt*HeatTreatmt + (1|nestID)", sep = ",")
library(glmmTMB)
for (y in dvnames2){
  form <- formula(paste(y, "~", ivnames2))
  models2[[y]] <- glmmTMB(form, data= EABL, family = nbinom2(link = "log"))
}

lapply(models2, summary) ## summarize each model

######################### WITHIN BROOD VARIATION BY SPECIES ########################
require(phyloseq)
require(ggplot2)
require(cowplot)

bluebirds.richnessbynest<-plot_richness(physeq = bluebirds, x = "nestID", measures = "Shannon") + geom_boxplot(aes(x=nestID, y=value, color=NULL), alpha=0.1)
treeswallows.richnessbynest<-plot_richness(physeq = treeswallows, x = "nestID", measures = "Shannon") + geom_boxplot(aes(x=nestID, y=value, color=NULL), alpha=0.1)

pdf(file = "Nestlings_withinnest_alpha.pdf", width = 6, height = 7)
plot_grid(x = bluebirds.richnessbynest, treeswallows.richnessbynest,nrow = 2)
dev.off()

# Calculate distance matrix
EABL_bray <- phyloseq::distance(bluebirds.trans, method = "bray")
TRES_bray<- phyloseq::distance(treeswallows.trans, method = "bray")

# make a data frame from the sample_data
EABLsampledf <- data.frame(sample_data(bluebirds.trans))
TRESsampledf<- data.frame(sample_data(treeswallows.trans))

EABL.ord<-ordinate(physeq = bluebirds.trans, method = "NMDS", distance = "bray")
TRES.ord<-ordinate(physeq = treeswallows.trans, method = "NMDS", distance = "bray")

EABL.NMDS<-plot_ordination(
  physeq = bluebirds.trans,
  ordination = EABL.ord,
  axes = c(1,2), 
  color = "nestID",
  title = "EABL by Nest") +
  theme_classic() +
  geom_jitter()

TRES.NMDS<-plot_ordination(
  physeq = treeswallows.trans,
  ordination = TRES.ord,
  axes = c(1,2), 
  color = "nestID",
  title = "TRES by Nest") +
  theme_classic() +
  geom_jitter()
dev.off()

pdf(file = "NMDS_Nests.pdf", width = 12, height = 6)
plot_grid(x = EABL.NMDS, TRES.NMDS,nrow = 1)
dev.off()
