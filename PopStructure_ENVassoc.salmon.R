#load libraries
#data manip and plot
library(data.table)  
library(tidyverse)
library(ggman)
library(marmap)
library(parallel)

#pop genomics
library(pcadapt)
library(LEA)
library(lfmm)
library(qvalue)
library(vegan)
library(corrplot)


#set a working directory
setwd("~/Desktop/Working/SalmonENV_tests/")

#### data required ####
##plink format files - .bed 
                      ##.ped in 12 format with --recode 12 option in plink, 
                      ###.raw format --recode A option in plink
                      ##.fam with individual info
                      ### .bim for marker information

##metadata with environmental coordinates, 

##for polysel - ssalar_genes_result.txt gene info and ssalar_kegg.csv kegg pathways, biosystems_gene file

###Data import, cleaning, and mapping ####

#read in metadata
salmon.env.meta <- fread("Salmo_220K_Merged2022_ENVmatchyear_Metadata.txt")


#match metadata to individual IDs in order of genotype file

salmon.env.meta.fam <- fread("Salmo_220K_Merged2022_COMPLETE_ENVmatchyear.fam") %>%  mutate(SiteCode = V1, ID = V2) %>%  
  select(SiteCode, ID) %>% 
  inner_join(., salmon.env.meta)

#do any filtering required - here removing a batch of individuals non-randomly sampled and without DU, and making uniq for mapping
salmon.env.meta.map <- salmon.env.meta.fam %>%  
  filter (!Batch %in% "April 2022 Mixed Samples Array", !DU.COSEWIC.2010 %in% "NA") %>%  
  select(-ID) %>%  distinct()


#Map
bathydata <- getNOAA.bathy(-71, -43, 40, 57, res=2,keep=T)
plot(bathydata)
map=autoplot(bathydata, geom=c("r", "c"), colour="grey", size=0.1) +
  scale_fill_etopo(guide = FALSE) +
  geom_contour(aes(z=z), breaks=c(-100, -200, -500, -1000, -2000, -4000), colour="grey90", size=0.01) +
  xlab("Degrees Longitude") +
  ylab("Degrees Latitude") +
  theme(legend.position = "none")

map + geom_point(data = salmon.env.meta.map, aes(x = as.numeric(Lon), y =as.numeric(Lat), colour = DU.COSEWIC.2010),  size = 3, inherit.aes = F) + theme(legend.position="right") 
  

#get SNP map
chrom_snp_map <- fread("Salmo_220K_Merged2022_ENVmatchyear.bim") %>%  
  mutate(CHROM =V1, SNP = V2, BP = V4) %>%  
  select(CHROM, BP, SNP)

#get alt chromosome IDs
chrom_snp_map <- fread("ssa_chromconvert.txt") %>% 
  inner_join(., chrom_snp_map)

#Population structure inference #####

##PCA ####
#PCA functions
PCADAPT_import_PCA <- function(filename, Knum) {
  to_pc <- read.pcadapt(paste0(filename, ".bed"), type = "bed")
  PCA <- pcadapt(to_pc, K = Knum, min.maf = 0.01)
  return(PCA)
} # a function to import data and run a PCA - currently set to maf 0.01

PCADAPT_scores_meta <- function(PCAobj, filename, Knum, meta) {
  FAM <- fread(paste0(filename, ".fam")) %>% 
    select(V1,V2) %>% 
    mutate(SiteCode = V1, ID = V2) %>%  
    select(-V1, -V2)
  meta_ordered <- inner_join(FAM, meta)
  PCAscores <- data.frame(PCAobj$scores)
  colnames(PCAscores) <- paste0("PC", rep(1:Knum))
  PCA_meta <- bind_cols(meta_ordered, PCAscores)
  return(PCA_meta)
} ## a function to match PCA scores to metadata

PCADAPT_var_per_axis <- function(PCAobj) {
  var_per_axis <- (PCAobj$singular.values^2)
  return(var_per_axis)} #a function to get PCA per axis variation

#run PCA  - different K values to check scree plot, explore different associations with axes 
PC_K5 <- PCADAPT_import_PCA(filename = "Salmo_220K_Merged2022_COMPLETE_ENVmatchyear", Knum = 5)
PC_K20 <- PCADAPT_import_PCA(filename = "Salmo_220K_Merged2022_ENVmatchyear", Knum = 20)

#screeplot
plot(PC_K20, option = "screeplot") + theme_classic()  #Check plateau, here K = 5 

##Population structure analysis with PCA
PC_K5_meta <-  PCADAPT_scores_meta(PCAobj = PC_K5, filename =  "Salmo_220K_Merged2022_COMPLETE_ENVmatchyear", Knum = 5, meta = salmon.env.meta.fam)


#plot spatial patterns in PCA
ggplot() +
  geom_point(data = PC_K5_meta, aes(x = PC1, y = PC2, colour = as.numeric(Lat))) +
  theme_classic() + scale_colour_gradient(low = "blue", high = "red")

ggplot() +
  geom_point(data = PC_K5_meta, aes(x = PC1, y = PC2, colour = as.numeric(Lon))) +
  theme_classic() + scale_colour_gradient(low = "blue", high = "red")

#plot by DU
ggplot() +
  geom_point(data = PC_K5_meta, aes(x = PC1, y = PC2, colour = DU.COSEWIC.2010)) +
  theme_classic()


#snmf admixture proportions ####
ped2lfmm(input.file = "Salmo_220K_Merged2022_COMPLETE_ENVmatchyear.ped") #conversion only needs to be done once


envmatch.lea = NULL
envmatch.lea= snmf(input.file = "Salmo_220K_Merged2022_COMPLETE_ENVmatchyear.lfmm", K = 1:10, entropy = TRUE, project = "new", CPU = 16 )

#plot CV scores
plot(envmatch.lea, col = "blue", pch = 19, cex = 1.2)  #big declines until ~6? will use 5 for consistency, could use 10.. art not a science, hierarchical structure is likely

snmf_Qscores_meta <- function(leaproj, Knum, meta) {
  leaQ <-Q(object = leaproj, K = Knum, run = 1)
  colnames(leaQ) <- paste0("Q", rep(1:Knum))
  leaQ_meta <- bind_cols(meta, leaQ)
  return(leaQ_meta) }  # a function to import Q scores from snmf and match to metadata

envmatch.snmf.K5.meta<- snmf_Qscores_meta(leaproj = envmatch.lea, Knum = 5, meta = salmon.env.meta.fam)

Make_admix_table<-  function(Admix_meta_object, Knum){Admix_table <- Admix_meta_object
rownames(Admix_table) <-Admix_table$ID
plot_data <-  Admix_table  %>% 
  gather('pop', 'prob', Q1:paste0("Q", Knum))
return(plot_data)} #a function to make a plottable table for snmf 

envmatch.snmf.K5.admixtable <- Make_admix_table(envmatch.snmf.K5.meta, Knum = 5)

ggplot(envmatch.snmf.K5.admixtable, aes(ID, prob, fill = pop)) +
  geom_col() +
  facet_grid(~DU.COSEWIC.2010, scales = 'free', space = 'free')

#RDA - grilse proportion
#get geno matrix for RDA
genos.dose <- fread("Salmo_220K_Merged2022_COMPLETE_ENVmatchyear.raw", sep = " ") %>%  select (-FID, -IID, -PAT, -MAT, -SEX, -PHENOTYPE)
genos.dose.imp<- apply(genos.dose, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))

#Check predictor correlation
salmon.env.bysite <- salmon.env.meta.fam %>%  select(Minimum.Air.Temperature, 
                                Air.Temperature,
                                Maximum.Air.Temperature,
                                Total.Precipitation,
                                Dew.Point.Temperature, 
                                Relative.Humidity,
                                Wind.Direction,
                                Solar.Radiation,
                                Snow.Precipitation,
                                Snow.Depth.Accumulation,
                                Snow.Water.Equivalent,
                                Wind.Speed.at.2.meters) %>%  distinct()



ENVcor <- cor(salmon.env.bysite)
corrplot(ENVcor, method = 'number') # visualize variable correlation

#remove correlated variables - some room here for picking favorites rather than doing things algorithmically
salmon.env_uncorrelated.meta.fam <- salmon.env.meta.fam %>%  select(Maximum.Air.Temperature, 
                                                                Total.Precipitation, 
                                                                Wind.Direction, 
                                                                Solar.Radiation, 
                                                                Snow.Water.Equivalent, 
                                                                Wind.Speed.at.2.meters )
##RDA
env.rda <- rda(genos.dose.imp ~ salmon.env_uncorrelated.meta.fam$Maximum.Air.Temperature +
              salmon.env_uncorrelated.meta.fam$Total.Precipitation +
              salmon.env_uncorrelated.meta.fam$Wind.Direction + 
              salmon.env_uncorrelated.meta.fam$Solar.Radiation + 
              salmon.env_uncorrelated.meta.fam$Snow.Water.Equivalent + 
              salmon.env_uncorrelated.meta.fam$Wind.Speed.at.2.meters)

#adjusted R2 and significance check 
RsquareAdj(env.rda)
anova.cca(env.rda, parallel= 16)
