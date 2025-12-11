# Purpose: Statistical analysis changes in truffle-like fungal community composition between seasons and with changes in climate 

# Load packages 


library(phyloseq)
library(ggplot2)
library(vegan)
library(patchwork)
library(pairwiseAdonis)
library(tibble)
library(indicspecies)
library(ggrepel)
library(pscl)
library(emmeans)
library(MASS)
library(dendextend)
library(dplyr)
library(geosphere)
library(tidyr)

# 4 datasets:
## 1: Samples from 5 sites between 1993-95
## 2: Samples from Watchmaker and Bellbird between 1998-2005
## 3: Samples from Watchmaker to test LFP traits
## 4: Samples from Bellbird with topographic data 


# 1. Samples at all sites between 1993-1995 ----


## 1.1. Ectomycorrhizal fungi ----

data <- prune_samples(sample_data(ITSrel.count2)$Year < 1996, ITSrel.count2)
data <- prune_samples(!(sample_data(data)$Location %in% c("Ellery Two","Puggaree Rd")), data) # Remove sites with < 5 samples

ecm.list=row.names(guild_table)[which(guild_table$primary_lifestyle=="ectomycorrhizal")] # Filter for ECM taxa
data_ecm=prune_taxa(ecm.list, data)


data_ecm = prune_taxa(taxa_sums(data_ecm) > 0, data_ecm) # Remove OTUs with zero count
ecm <- prune_samples(sample_sums(data_ecm) > 0, data_ecm) # Remove samples with zero OTUs 

### 1.1.1 Ordinations  ----


iDist <- distance((ecm), method="bray")
iMDS  <- ordinate(ecm, "NMDS", distance=iDist)

data.scores <- as.data.frame(scores(iMDS))
sort(data.scores$NMDS1) 

ecm <- subset_samples(ecm, !agrf. %in% c("1", "52","29")) # Remove outliers that obscure patterns in ordination

iDist <- distance((ecm), method="bray")
iMDS  <- ordinate(ecm, "NMDS", distance=iDist)

# Ordinate samples by Location
ordination_df <- as.data.frame(scores(iMDS, display = "sites"))
ord_ecm <- cbind(ordination_df, sample_data(ecm))

location_colours <- c("#D1B08F", "#95A58E","#EAA48F","#92c9c4","#5EA69D")


ordplot <- ggplot(ord_ecm, aes(x = NMDS1, y = NMDS2, color = Location)) +
  geom_point(size = 3) +
  scale_colour_manual(values = location_colours)+
  theme_light() +
  theme(aspect.ratio = 1) +
  ggtitle("Spatial Differences in Community Composition of ECM in LFP scats collected in 1993-1995") +
  theme(plot.title = element_text(size = 16))
ordplot

find_hull <- function(df) df[chull(df$NMDS1, df$NMDS2), ]

# Calculate the polygons for each site
hulls <- ord_ecm %>%
  group_by(Location) %>%
  do(find_hull(.))

# Add polygons to ordination plot
ordplot + 
  geom_polygon(data = hulls, aes(x = NMDS1, y = NMDS2, fill = Location, group = Location), alpha = 0.05)  +
  theme(plot.title = element_text(size = 13)) +scale_fill_manual(values = location_colours) + theme_light()




# Ordinate samples by Region
ordination_df <- as.data.frame(scores(iMDS, display = "sites"))
ord_ecm <- cbind(ordination_df, sample_data(ecm))

location_colours <- c("#D1B08F", "#95A58E","#EAA48F","#92c9c4","#5EA69D")


ordplot <- ggplot(ord_ecm, aes(x = NMDS1, y = NMDS2, color = Region)) +
  geom_point(size = 3) +
  scale_colour_manual(values = location_colours)+
  theme_light() +
  theme(aspect.ratio = 1) +
  ggtitle("Regional Differences in Community Composition of \nECM in LFP scats collected in 1993-1995") +
  theme(plot.title = element_text(size = 16))
ordplot

find_hull <- function(df) df[chull(df$NMDS1, df$NMDS2), ]

# Calculate the polygons for each site
hulls <- ord_ecm %>%
  group_by(Region) %>%
  do(find_hull(.))

# Add polygons to ordination plot
ordplot + 
  geom_polygon(data = hulls, aes(x = NMDS1, y = NMDS2, fill = Region, group = Region), alpha = 0.05)  +
  theme(plot.title = element_text(size = 13)) +scale_fill_manual(values = location_colours) + theme_light()



### 1.1.2 PERMANOVAs  ----

sampledf <- data.frame(sample_data(ecm))
otu_table <- as(t(otu_table(ecm)), "matrix")
bray_dist <- vegdist(otu_table, method = "bray")

# PERMANOVA to test for differences in community composition between locations and seasons

permanova_result <- adonis2(bray_dist ~ Location, data = sampledf)
print(permanova_result)  

permanova_result <- adonis2(bray_dist ~ Location * Season, data = sampledf)
print(permanova_result) # Significant interaction between location and season 


# Pairwise differences in community composition between locations 
pairwise.adonis(bray_dist, phyloseq::sample_data(ecm)$Location)

# Test for difference in group dispersion between location
beta_location <- betadisper(bray_dist, sampledf$Location)
permutest(beta_location) # Not significant
plot(beta_location, hull=FALSE, ellipse=TRUE) 
boxplot(beta_location)

# Test for difference in group dispersion between season
beta_season <- betadisper(bray_dist, sampledf$Season)
permutest(beta_season) # Dispersion differs significantly between seasons
plot(beta_season, hull=FALSE, ellipse=TRUE) 
boxplot(beta_season)



### 1.1.3 Indicspecies ----

# Test for genera that are indicative of particular sample groupings

sampledf <- data.frame(sample_data(ecm))
data_genus <- tax_glom(ecm, taxrank = "Genus") # Group taxa at the level of genus

otu <- as.data.frame(t(otu_table(data_genus)))
tax <- as.data.frame(t(tax_table(data_genus))) %>% rownames_to_column(var = "Taxa")
genus <- tax %>% filter(Taxa == "Genus")
genus <- genus[-1]
genus_otu <- as.data.frame(rbind(genus, otu))
colnames(genus_otu) <- genus_otu[1,]
genus_otu <- genus_otu[-1,]
genus_otu <- genus_otu %>%
  mutate(across(everything(), as.numeric))

genus_otu_bind <- cbind(sampledf, genus_otu)
Location <- (sample_data(genus_otu_bind)$Location) # Replace with $Season or $Region to test for genera that are indicative of differences between these groupings 


indicator_taxa <- multipatt(genus_otu, Location, func = "r", 
                 control = how(nperm=999)) 
summary(indicator_taxa)


### 1.1.4 DB-RDA ----

data_genus <- tax_glom(ecm, taxrank = "Genus") # Group taxa at the level of genus
otu_table <- as(t(otu_table(data_genus)), "matrix") 
sample_data <- data.frame(sample_data(data_genus))

col_names_indices <- data.frame(Index = seq_along(colnames(sample_data)), Column = colnames(sample_data))
print(col_names_indices)

sample_data[9:29] <- lapply(sample_data[9:29], function(x) as.numeric(as.character(x)))
sample_data[9:29]  <- scale(sample_data[9:29]) # Scale climate variables


bray_distance <- vegdist(otu_table, method = "bray")

#use ordistep to chose model by permutation tests.

null<- capscale(bray_distance ~ 1, data = sample_data, distance = "bray")

# Prepare a model with all non correlated climate variables (selected in script 2.alpha_diversity.R , section 1.1.4)
 

full_model <- capscale(bray_distance ~ Location + Region + Season + week_precip + three_mon_mint_av + three_mon_AHMI + twelve_mon_precip + twelve_mon_mint_av + wk_av_temp + mon_precip + mon_AHMI, data = sample_data, distance = "bray")

# Select model
model_select <- ordiR2step(null, scope = formula(full_model), R2scope = TRUE, pstep = 0.05)

formula(model_select)
print(model_select)
summary(model_select)

# Check for VIF
vif_values <- vif.cca(model_select)
print(vif_values)
                            
marginal_effects <- anova(model_select, by = "margin")
marginal_effects

anova_overall <- anova.cca(model_select, permutations = 999)
print(anova_overall)
db_rda_scores <- scores(model_select, display = "sites")
db_rda_env <- scores(model_select, display = "bp")
scores <- as.data.frame(db_rda_env)

# Combine the scores with the sample data
sites_data <- data.frame(db_rda_scores, sample_data) 
env_data <- data.frame(db_rda_env)

location_colours <- c("#D1B08F", "#95A58E","#EAA48F","#92c9c4","#5EA69D")

rownames(env_data)
env_data$name <- c("Spring","Summer","Winter", "Ellery One", "Riley", "West Buffalo", "Quarterly MINT")

p <- ggplot(sites_data, aes(x = CAP1, y = CAP2, color = Location)) +
  geom_point(size = 3) 

find_hull <- function(df) df[chull(df$CAP1, df$CAP2), ]

# Calculate polygon for each site
hulls <- sites_data %>%
  group_by(Location) %>%
  do(find_hull(.))

# Add polygons to the plot
p <- p + 
  geom_polygon(data = hulls, aes(x = CAP1, y = CAP2, fill = Location, group = Location), alpha = 0.02) 
p +
  scale_fill_manual(values = location_colours) + scale_color_manual(values = location_colours)+
  geom_segment(data = env_data, aes(x = 0, xend = CAP1, y = 0, yend = CAP2),
               arrow = arrow(length = unit(0.2, "cm")), color = "#000000", size = 0.8) +
  geom_text_repel(data = env_data, aes(x = CAP1, y = CAP2, label = name), 
                  color = "#000000", hjust =0, vjust = 0.0 ,force_pull = -0.06,
                  segment.size = 0.5,  # Controls the thickness of the lines
                  segment.linetype = "dashed", size = 4) + 
  theme_minimal() +
  theme(legend.position = "top") +
  labs(title = "db-RDA for ECM communities in LFP scats from 1993-1995", x = "CAP1", y = "CAP2") 



### 1.1.5 Cluster analysis -----

## Presence absence
ecm_pa <- decostand(otu_table(ecm),method="pa")
binary <- phyloseq(otu_table(ecm_pa, taxa_are_rows=TRUE), 
                   sample_data(ecm), 
                   tax_table(ecm))
otu_table_bin <- as(t(otu_table(binary)), "matrix") 
samp_table_bin <- data.frame(sample_data(binary))
otu_table_bin <- data.frame(otu_table_bin)
otu_table_bin$SampleID <- row.names(otu_table_bin)
samp_table_bin$SampleID <- row.names(samp_table_bin)
join <- left_join(samp_table_bin,otu_table_bin, by = "SampleID")

col_names_indices <- data.frame(Index = seq_along(colnames(join)), Column = colnames(join))
print(col_names_indices)
join<- join[-c(1,2,4:35)] # Remove columns that are not Location or OTUs
join

# Summarise OTUs by location
join <- join %>% group_by(Location)  %>% summarise(across(where(is.numeric), max)) %>% 
  column_to_rownames("Location")

raupcrick_location <- raupcrick(join)

site_clust <- hclust(raupcrick_location, method = "complete")
ecm_pa_clust <- plot(site_clust, main = "Hierarchical Clustering of Sites by ECM OTU _p/a ")


dend <- as.dendrogram(site_clust)
labels <- labels(dend)
group1 <- c("Bellbird", "Ellery One") #East Gippsland sites
group2 <- c("Riley", "West Buffalo") #North East Victoria sites

label_colors <- ifelse(labels %in% group1, "#5EA69D",ifelse(labels %in% group2, "#EAA48F", "black")) # Colour by Region 

dend <- dend %>% set("labels_col", label_colors) %>% set("labels_cex", 0.9)
par(font = 1)
plot(dend, main = "ECM OTU Clustering by Site - Presence/Absence")


## Relative abundance
otu_table <- as(t(otu_table(ecm)), "matrix") 
sample_data <- data.frame(sample_data(ecm))

otu_df <- data.frame(otu_table)
otu_df$SampleID <- row.names(otu_df)
sample_data$SampleID <- row.names(sample_data)
join <- left_join(sample_data,otu_df, by = "SampleID")

col_names_indices <- data.frame(Index = seq_along(colnames(join)), Column = colnames(join))
print(col_names_indices)

join<- join[-c(1,2,4:35)] # Remove columns that are not Location or OTUs
join


# Summarise OTUs by location
join <- join %>% group_by(Location)  %>% summarise(across(where(is.numeric), sum)) %>% 
  column_to_rownames("Location")


bray_distance_location <- vegdist(join, method = "bray")

site_clust <- hclust(bray_distance_location, method = "complete")
ecm_abund_clust <- plot(site_clust, main = "Hierarchical Clustering of Sites by ECM OTU ")


dend <- as.dendrogram(site_clust)
labels <- labels(dend)
group1 <- c("Bellbird", "Ellery One") #East Gippsland sites
group2 <- c("Riley", "West Buffalo") #North East Victoria sites
label_colors <- ifelse(labels %in% group1, "#5EA69D",ifelse(labels %in% group2, "#EAA48F", "black"))

dend <- dend %>% set("labels_col", label_colors) %>% set("labels_cex", 0.9)
par(font = 1)
plot(dend, main = "ECM OTU Clustering by Site - Relative Abundance")





### 1.1.6 Mantel correlation -----
site_coords <- read.csv("data/site_coords.csv")

EG_ecm <- prune_samples((sample_data(ecm)$Region %in% c("East Gippsland")), ecm) # EG = East Gippsland sub-population
NE_ecm <- prune_samples(!(sample_data(ecm)$Region %in% c("East Gippsland")), ecm) # NE = Northeast Victoria sub-population
 
region <- EG_ecm # swap to test EG or NE sites

sample_df <- data.frame(sample_data(region))
sample_df$SampleID <- rownames(sample_df)

# add coords to sample df
sample_df_location <- sample_df %>% left_join(site_coords, by = c("Location" = "data"))

rownames(sample_df_location) <- sample_df_location$SampleID
sample_data(region) <- sample_df_location

otu <- t(otu_table(region))

unique((sample_data(region))$Location)
# Jaccard for p-a data
#otu_pa <- (otu > 0) * 1
#dist <- vegdist(otu_pa, method = "jaccard", binary = TRUE)

# Bray for relative abundance data
dist <- vegdist(otu, method = "bray", binary = TRUE)


#geo distances
coords <- sample_df_location[, c("Longitude", "Latitude")]
dist_matrix <- distm(coords, fun = distHaversine)

distdist <- as.dist(dist_matrix)

unique(distdist)
unique(sample_df_location$Location)


mantel_ecm <- mantel.correlog(
  dist,
  distdist,
  n.class = 2,
  r.type = "spearman",
  nperm = 999
)

plot(mantel_ecm)
mantel_ecm$mantel.res




## 1.2. Truffle-like Ectomycorrhizal fungi ----

data <- prune_samples(sample_data(ITSrel.count2)$Year < 1996, ITSrel.count2)
data <- prune_samples(!(sample_data(data)$Location %in% c("Ellery Two","Puggaree Rd")), data) # Remove sites with < 5 samples

ecm.list=row.names(guild_table)[which(guild_table$primary_lifestyle=="ectomycorrhizal")] # Filter for ECM taxa
data_ecm=prune_taxa(ecm.list, data)


hypo.list=row.names(guild_table)[which(guild_table$Fruitbody_type=="gasteroid-hypogeous")]
data=prune_taxa(hypo.list, data_ecm)

data_truffle_ecm = prune_taxa(taxa_sums(data) > 0, data) # Remove OTUs with zero count
truffle_ecm <- prune_samples(sample_sums(data_truffle_ecm) > 0, data_truffle_ecm) # Remove samples with zero OTUs 

### 1.2.1 Ordinations  ----


iDist <- distance((truffle_ecm), method="bray")
iMDS  <- ordinate(truffle_ecm, "NMDS", distance=iDist)

data.scores <- as.data.frame(scores(iMDS))
sort(data.scores$NMDS1) 

truffle_ecm <- subset_samples(truffle_ecm, !agrf. %in% c("1", "52","18")) # Remove outliers that obscure patterns in ordination

iDist <- distance((truffle_ecm), method="bray")
iMDS  <- ordinate(truffle_ecm, "NMDS", distance=iDist)

# Ordinate samples by Location
ordination_df <- as.data.frame(scores(iMDS, display = "sites"))
ord_truffle_ecm <- cbind(ordination_df, sample_data(truffle_ecm))

location_colours <- c("#D1B08F", "#95A58E","#EAA48F","#92c9c4","#5EA69D")


ordplot <- ggplot(ord_truffle_ecm, aes(x = NMDS1, y = NMDS2, color = Location)) +
  geom_point(size = 3) +
  scale_colour_manual(values = location_colours)+
  theme_light() +
  theme(aspect.ratio = 1) +
  ggtitle("Spatial Differences in Community Composition of Truffle-like \nECM in LFP scats collected in 1993-1995") +
  theme(plot.title = element_text(size = 11))
ordplot

find_hull <- function(df) df[chull(df$NMDS1, df$NMDS2), ]

# Calculate the polygons for each site
hulls <- ord_truffle_ecm %>%
  group_by(Location) %>%
  do(find_hull(.))

# Add polygons to ordination plot
ordplot + 
  geom_polygon(data = hulls, aes(x = NMDS1, y = NMDS2, fill = Location, group = Location), alpha = 0.05)  +
  theme(plot.title = element_text(size = 13)) +scale_fill_manual(values = location_colours) + theme_light()




# Ordinate samples by Region
ordination_df <- as.data.frame(scores(iMDS, display = "sites"))
ord_truffle_ecm <- cbind(ordination_df, sample_data(truffle_ecm))

location_colours <- c("#D1B08F", "#95A58E","#EAA48F","#92c9c4","#5EA69D")


ordplot <- ggplot(ord_truffle_ecm, aes(x = NMDS1, y = NMDS2, color = Region)) +
  geom_point(size = 3) +
  scale_colour_manual(values = location_colours)+
  theme_light() +
  theme(aspect.ratio = 1) +
  ggtitle("Regional Differences in Community Composition of \nTruffle-like ECM in LFP scats collected in 1993-1995") +
  theme(plot.title = element_text(size = 13))
ordplot

find_hull <- function(df) df[chull(df$NMDS1, df$NMDS2), ]

# Calculate the polygons for each site
hulls <- ord_truffle_ecm %>%
  group_by(Region) %>%
  do(find_hull(.))

# Add polygons to ordination plot
ordplot + 
  geom_polygon(data = hulls, aes(x = NMDS1, y = NMDS2, fill = Region, group = Region), alpha = 0.05)  +
  theme(plot.title = element_text(size = 13)) +scale_fill_manual(values = location_colours) + theme_light()



### 1.2.2 PERMANOVAs  ----

sampledf <- data.frame(sample_data(truffle_ecm))
otu_table <- as(t(otu_table(truffle_ecm)), "matrix")
bray_dist <- vegdist(otu_table, method = "bray")



tax <- tax_table(truffle_ecm)
tax <- data.frame(tax)

taxcount <- tax %>% group_by(Genus, .drop = FALSE) %>% count()
taxcount
# PERMANOVA to test for differences in community composition between locations and seasons

permanova_result <- adonis2(bray_dist ~ Location, data = sampledf)
print(permanova_result)  

permanova_result <- adonis2(bray_dist ~ Location * Season, data = sampledf)
print(permanova_result) # Significant interaction between location and season 


# Pairwise differences in community composition between locations 
pairwise.adonis(bray_dist, phyloseq::sample_data(truffle_ecm)$Location)

# Test for difference in group dispersion between location
beta_location <- betadisper(bray_dist, sampledf$Location)
permutest(beta_location) # Not significant
plot(beta_location, hull=FALSE, ellipse=TRUE) 
boxplot(beta_location)

# Test for difference in group dispersion between season
beta_season <- betadisper(bray_dist, sampledf$Season)
permutest(beta_season) # Not significant
plot(beta_season, hull=FALSE, ellipse=TRUE) 
boxplot(beta_season)



### 1.2.3 Indicspecies ----

# Test for genera that are indicative of particular sample groupings

sampledf <- data.frame(sample_data(truffle_ecm))
data_genus <- tax_glom(truffle_ecm, taxrank = "Genus") # Group taxa at the level of genus

otu <- as.data.frame(t(otu_table(data_genus)))
tax <- as.data.frame(t(tax_table(data_genus))) %>% rownames_to_column(var = "Taxa")
genus <- tax %>% filter(Taxa == "Genus")
genus <- genus[-1]
genus_otu <- as.data.frame(rbind(genus, otu))
colnames(genus_otu) <- genus_otu[1,]
genus_otu <- genus_otu[-1,]
genus_otu <- genus_otu %>%
  mutate(across(everything(), as.numeric))

genus_otu_bind <- cbind(sampledf, genus_otu)
Location <- (sample_data(genus_otu_bind)$Location) # Replace with $Season or $Region to test for genera that are indicative of differences between these groupings 


indicator_taxa <- multipatt(genus_otu, Location, func = "r", 
                            control = how(nperm=999)) 
summary(indicator_taxa)


### 1.2.4 DB-RDA ----

data_genus <- tax_glom(truffle_ecm, taxrank = "Genus") # Group taxa at the level of genus
otu_table <- as(t(otu_table(data_genus)), "matrix") 
sample_data <- data.frame(sample_data(data_genus))

col_names_indices <- data.frame(Index = seq_along(colnames(sample_data)), Column = colnames(sample_data))
print(col_names_indices)

sample_data[9:29] <- lapply(sample_data[9:29], function(x) as.numeric(as.character(x)))
sample_data[9:29]  <- scale(sample_data[9:29]) # Scale climate variables


bray_distance <- vegdist(otu_table, method = "bray")

#use ordistep to chose model by permutation tests.

null<- capscale(bray_distance ~ 1, data = sample_data, distance = "bray")

# Prepare a model with all non correlated climate variables (selected in script 2.alpha_diversity.R , section 1.2.4)


full_model <- capscale(bray_distance ~ Location + Region + Season + week_precip + mon_maxt_av + three_mon_mint_av + three_mon_AHMI + twelve_mon_precip + twelve_mon_mint_av + mon_precip + mon_AHMI, data = sample_data, distance = "bray")

full_model <- capscale(bray_distance ~ Season + Location + three_mon_mint_av, data = sample_data, distance = "bray")

# Select model
model_select <- ordiR2step(null, scope = formula(full_model), R2scope = TRUE, pstep = 0.05) # remove twelve_mon_mint_av (high VIF), run model again

model_select <- capscale(bray_distance ~ Season + Location + three_mon_mint_av , data = sample_data, distance = "bray")


formula(model_select)
print(model_select)
summary(model_select)

# Check for VIF
vif_values <- vif.cca(model_select)
print(vif_values) # Remove variables with VIF > 10

marginal_effects <- anova(model_select, by = "margin")
marginal_effects

anova_overall <- anova.cca(model_select, permutations = 999)
print(anova_overall)
db_rda_scores <- scores(model_select, display = "sites")
db_rda_env <- scores(model_select, display = "bp")
scores <- as.data.frame(db_rda_env)

# Combine the scores with the sample data
sites_data <- data.frame(db_rda_scores, sample_data) 
env_data <- data.frame(db_rda_env)

location_colours <- c("#D1B08F", "#95A58E","#EAA48F","#92c9c4","#5EA69D")

rownames(env_data)
env_data$name <- c("Spring","Summer","Winter", "Ellery One", "Riley", "West Buffalo")

p <- ggplot(sites_data, aes(x = CAP1, y = CAP2, color = Location)) +
  geom_point(size = 3) 

find_hull <- function(df) df[chull(df$CAP1, df$CAP2), ]

# Calculate polygon for each site
hulls <- sites_data %>%
  group_by(Location) %>%
  do(find_hull(.))

# Add polygons to the plot
p <- p + 
  geom_polygon(data = hulls, aes(x = CAP1, y = CAP2, fill = Location, group = Location), alpha = 0.02) 
p +
  scale_fill_manual(values = location_colours) + scale_color_manual(values = location_colours)+
  geom_segment(data = env_data, aes(x = 0, xend = CAP1, y = 0, yend = CAP2),
               arrow = arrow(length = unit(0.2, "cm")), color = "#000000", size = 0.8) +
  geom_text_repel(data = env_data, aes(x = CAP1, y = CAP2, label = name), 
                  color = "#000000", hjust =0, vjust = 0.0 ,force_pull = -0.06,
                  segment.size = 0.5,  # Controls the thickness of the lines
                  segment.linetype = "dashed", size = 4) + 
  theme_minimal() +
  theme(legend.position = "top") +
  labs(title = "db-RDA for Truffle-like ECM communities \nin LFP scats from 1993-1995", x = "CAP1", y = "CAP2") 



### 1.2.5 Cluster analysis -----

## Presence absence
truffle_ecm_pa <- decostand(otu_table(truffle_ecm),method="pa")
binary <- phyloseq(otu_table(truffle_ecm_pa, taxa_are_rows=TRUE), 
                   sample_data(truffle_ecm), 
                   tax_table(truffle_ecm))
otu_table_bin <- as(t(otu_table(binary)), "matrix") 
samp_table_bin <- data.frame(sample_data(binary))
otu_table_bin <- data.frame(otu_table_bin)
otu_table_bin$SampleID <- row.names(otu_table_bin)
samp_table_bin$SampleID <- row.names(samp_table_bin)
join <- left_join(samp_table_bin,otu_table_bin, by = "SampleID")

col_names_indices <- data.frame(Index = seq_along(colnames(join)), Column = colnames(join))
print(col_names_indices)
join<- join[-c(1,2,4:35)] # Remove columns that are not Location or OTUs
join

# Summarise OTUs by location
join <- join %>% group_by(Location)  %>% summarise(across(where(is.numeric), max)) %>% 
  column_to_rownames("Location")

raupcrick_location <- raupcrick(join)

site_clust <- hclust(raupcrick_location, method = "complete")
truffle_ecm_pa_clust <- plot(site_clust, main = "Hierarchical Clustering of Sites by truffle_ecm OTU _p/a ")


dend <- as.dendrogram(site_clust)
labels <- labels(dend)
group1 <- c("Bellbird", "Ellery One") #East Gippsland sites
group2 <- c("Riley", "West Buffalo") #North East Victoria sites

label_colors <- ifelse(labels %in% group1, "#5EA69D",ifelse(labels %in% group2, "#EAA48F", "black")) # Colour by Region 

dend <- dend %>% set("labels_col", label_colors) %>% set("labels_cex", 0.9)
par(font = 1)
plot(dend, main = "Truffle-like ECM OTU Clustering by Site - Presence/Absence")


## Relative abundance
otu_table <- as(t(otu_table(truffle_ecm)), "matrix") 
sample_data <- data.frame(sample_data(truffle_ecm))

otu_df <- data.frame(otu_table)
otu_df$SampleID <- row.names(otu_df)
sample_data$SampleID <- row.names(sample_data)
join <- left_join(sample_data,otu_df, by = "SampleID")

col_names_indices <- data.frame(Index = seq_along(colnames(join)), Column = colnames(join))
print(col_names_indices)

join<- join[-c(1,2,4:35)] # Remove columns that are not Location or OTUs
join


# Summarise OTUs by location
join <- join %>% group_by(Location)  %>% summarise(across(where(is.numeric), sum)) %>% 
  column_to_rownames("Location")


bray_distance_location <- vegdist(join, method = "bray")

site_clust <- hclust(bray_distance_location, method = "complete")
truffle_ecm_abund_clust <- plot(site_clust, main = "Hierarchical Clustering of Sites by truffle_ecm OTU ")


dend <- as.dendrogram(site_clust)
labels <- labels(dend)
group1 <- c("Bellbird", "Ellery One") #East Gippsland sites
group2 <- c("Riley", "West Buffalo") #North East Victoria sites
label_colors <- ifelse(labels %in% group1, "#5EA69D",ifelse(labels %in% group2, "#EAA48F", "black"))

dend <- dend %>% set("labels_col", label_colors) %>% set("labels_cex", 0.9)
par(font = 1)
plot(dend, main = "Truffle-like ECM OTU Clustering by Site - Relative Abundance")




### 1.2.6 Mantel correlation -----
site_coords <- read.csv("data/site_coords.csv")

EG_ecm <- prune_samples((sample_data(truffle_ecm)$Region %in% c("East Gippsland")), truffle_ecm) # EG = East Gippsland sub-population
NE_ecm <- prune_samples(!(sample_data(truffle_ecm)$Region %in% c("East Gippsland")), truffle_ecm) # NE = Northeast Victoria sub-population

region <- EG_ecm # swap to test EG or NE sites

sample_df <- data.frame(sample_data(region))
sample_df$SampleID <- rownames(sample_df)

# add coords to sample df
sample_df_location <- sample_df %>% left_join(site_coords, by = c("Location" = "data"))

rownames(sample_df_location) <- sample_df_location$SampleID
sample_data(region) <- sample_df_location

otu <- t(otu_table(region))

unique((sample_data(region))$Location)
# Jaccard for p-a data
#otu_pa <- (otu > 0) * 1
#dist <- vegdist(otu_pa, method = "jaccard", binary = TRUE)

# Bray for relative abundance data
dist <- vegdist(otu, method = "bray", binary = TRUE)


#geo distances
coords <- sample_df_location[, c("Longitude", "Latitude")]
dist_matrix <- distm(coords, fun = distHaversine)

distdist <- as.dist(dist_matrix)

unique(distdist)
unique(sample_df_location$Location)


mantel_ecm_truffle <- mantel.correlog(
  dist,
  distdist,
  n.class = 2,
  r.type = "spearman",
  nperm = 999
)

plot(mantel_ecm_truffle)
mantel_ecm_truffle$mantel.res


# Plot east Gippsland mantel correlations for both ECM and truffle-like ECM 

eastgipps_truffle <- as.data.frame(mantel_ecm_truffle$mantel.res)
eastgipps_truffle$taxa <- "Truffle-like ECM"
eastgipps_ecm <- as.data.frame(mantel_ecm$mantel.res)
eastgipps_ecm$taxa <- "ECM"
combined <- rbind(eastgipps_truffle, eastgipps_ecm)



ggplot(combined, aes(x = class.index, y = Mantel.cor, color = taxa, group = taxa)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line() +
  geom_point(aes(shape = 'Pr(corrected)'), size = 2,show.legend = FALSE) + #add points to figure
  scale_color_manual(values = c( "#5EA69D", "#EAA48F")) +
  scale_x_continuous(breaks = combined$class.index,labels = scales::label_number(accuracy = 1, big.mark = ""))+ # label distance class on x axis 
  labs(x = "Midpoint of distance class (m)", y = "Mantel correlation",
       color = NULL) + theme_minimal()


### 1.2.7  Regional occurence of genera ---- 

genus <- tax_glom(truffle_ecm, taxrank = "Genus") # glom at level of Genus
genus_pa <- transform_sample_counts(genus, function(x) ifelse(x > 0, 1, 0)) # convert to presence-absence

region_presence <- merge_samples(genus_pa, "Region")
region_presence <- transform_sample_counts(region_presence, function(x) ifelse(x > 0, 1, 0))


otu <- as.data.frame(t(otu_table(region_presence)))
tax <- as.data.frame(tax_table(region_presence))

otu$Genus <- tax$Genus

genus_long <- pivot_longer(otu, cols = c("East Gippsland", "North East Victoria"), 
                        names_to = "Region", values_to = "Presence")

genus_long<- genus_long %>%
  mutate(Region = if_else(Region == "North East Victoria", "Northeast Victoria", Region))


genus_long$Genus <- factor(genus_long$Genus, levels = rev(sort(unique(genus_long$Genus)))) #sorts alphabetically

ggplot(genus_long, aes(x = Region, y = Genus, fill = factor(Presence))) +
  geom_tile(color = "black")  +
    scale_fill_manual(values = c("0" = "white", "1" = "#778773"), name = "Presence", labels = c("Absent", "Present"))+
  theme_minimal() +
  labs(title = "Presence of Truffle-like ECM genera between regions",
       x = "Region",
       y = NULL) +
  theme(axis.text.y = element_text(size = 9), legend.title = element_blank(),axis.title.x = element_text(size = 10))



# 2. Samples at Watchmaker and Bellbird between 1998-2005 ----


## 2.1. Ectomycorrhizal fungi ----

ecm.list=row.names(guild_table)[which(guild_table$primary_lifestyle=="ectomycorrhizal")]
data=prune_taxa(ecm.list, ITSrel.count2)
data <- prune_samples(sample_data(data)$Year >= 1998 & sample_data(data)$Year <= 2005, data) # Select years where samples were collected at both Bellbird & Watchmaker
data_ecm <- subset_samples(data, Location %in% "Bellbird" | Location %in% "Watchmaker" )


data_ecm = prune_taxa(taxa_sums(data_ecm) > 0, data_ecm) # Remove OTUs with zero count
ecm <- prune_samples(sample_sums(data_ecm) > 0, data_ecm) # Remove samples with zero OTUs 

### 2.1.1 Ordinations  ----


iDist <- distance((ecm), method="bray")
iMDS  <- ordinate(ecm, "NMDS", distance=iDist)

data.scores <- as.data.frame(scores(iMDS))
sort(data.scores$NMDS1) 

iDist <- distance((ecm), method="bray")
iMDS  <- ordinate(ecm, "NMDS", distance=iDist)

# Ordinate samples by Location
ordination_df <- as.data.frame(scores(iMDS, display = "sites"))
ord_ecm <- cbind(ordination_df, sample_data(ecm))

location_colours  <- c("#EAA48F","#92c9c4")


ordplot <- ggplot(ord_ecm, aes(x = NMDS1, y = NMDS2, color = Location)) +
  geom_point(size = 3) +
  scale_colour_manual(values = location_colours)+
  theme_light() +
  theme(aspect.ratio = 1) +
  ggtitle("Spatial Differences in Community Composition of ECM in LFP scats collected at Watchmaker and Bellbird") +
  theme(plot.title = element_text(size = 16))
ordplot

find_hull <- function(df) df[chull(df$NMDS1, df$NMDS2), ]

# Calculate the polygons for each site
hulls <- ord_ecm %>%
  group_by(Location) %>%
  do(find_hull(.))

# Add polygons to ordination plot
ordplot + 
  geom_polygon(data = hulls, aes(x = NMDS1, y = NMDS2, fill = Location, group = Location), alpha = 0.05)  +
  theme(plot.title = element_text(size = 13)) +scale_fill_manual(values = location_colours) + theme_light()



### 2.1.2 PERMANOVAs  ----

sampledf <- data.frame(sample_data(ecm))
otu_table <- as(t(otu_table(ecm)), "matrix")
bray_dist <- vegdist(otu_table, method = "bray")

# PERMANOVA to test for differences in community composition between locations and seasons

permanova_result <- adonis2(bray_dist ~ Location, data = sampledf)
print(permanova_result)  

permanova_result <- adonis2(bray_dist ~ Location * Season, data = sampledf)
print(permanova_result) # Significant interaction between location and season 

# Test for difference in group dispersion between location
beta_location <- betadisper(bray_dist, sampledf$Location)
permutest(beta_location) # Dispersion differs significantly between locations
plot(beta_location, hull=FALSE, ellipse=TRUE) 
boxplot(beta_location)

# Test for difference in group dispersion between season
beta_season <- betadisper(bray_dist, sampledf$Season)
permutest(beta_season) # Dispersion differs significantly between seasons
plot(beta_season, hull=FALSE, ellipse=TRUE) 
boxplot(beta_season)



### 2.1.3 Indicspecies ----

# Test for genera that are indicative of particular sample groupings

sampledf <- data.frame(sample_data(ecm))
data_genus <- tax_glom(ecm, taxrank = "Genus") # Group taxa at the level of genus

otu <- as.data.frame(t(otu_table(data_genus)))
tax <- as.data.frame(t(tax_table(data_genus))) %>% rownames_to_column(var = "Taxa")
genus <- tax %>% filter(Taxa == "Genus")
genus <- genus[-1]
genus_otu <- as.data.frame(rbind(genus, otu))
colnames(genus_otu) <- genus_otu[1,]
genus_otu <- genus_otu[-1,]
genus_otu <- genus_otu %>%
  mutate(across(everything(), as.numeric))

genus_otu_bind <- cbind(sampledf, genus_otu)
Location <- (sample_data(genus_otu_bind)$Location) # Replace with $Season to test for genera that are indicative of differences between these groupings 


indicator_taxa <- multipatt(genus_otu, Location, func = "r", 
                            control = how(nperm=999)) 
summary(indicator_taxa)



### 2.1.4 DB-RDA ----

data_genus <- tax_glom(ecm, taxrank = "Genus") # Group taxa at the level of genus
otu_table <- as(t(otu_table(data_genus)), "matrix") 
sample_data <- data.frame(sample_data(data_genus))

col_names_indices <- data.frame(Index = seq_along(colnames(sample_data)), Column = colnames(sample_data))
print(col_names_indices)

sample_data[9:29] <- lapply(sample_data[9:29], function(x) as.numeric(as.character(x)))
sample_data[9:29]  <- scale(sample_data[9:29]) # Scale climate variables

## Identify highly correlated variables

clim_vars <- sample_data[, 9:29] # Select climate variables

cor_mat <- cor(clim_vars, use = "pairwise.complete.obs") # create the correlation matrix
high_cor <- findCorrelation(cor_mat, cutoff = 0.70, verbose = TRUE, names = TRUE) # Find high correlated pairs (correlation > 0.7)
non_cor_vars <- setdiff(names(clim_vars), high_cor) # List the climate variables that are not highly correlated

# create distance matrix
bray_distance <- vegdist(otu_table, method = "bray")

# use ordistep to chose model by permutation tests.
null<- capscale(bray_distance ~ 1, data = sample_data, distance = "bray")

# Prepare a model with all non correlated climate variables 

full_model <- capscale(bray_distance ~ Location + Season + week_precip + three_month_precip + wk_maxt_av + three_mon_mint_av + twelve_mon_precip + twelve_mon_mint_av + twelve_mon_av_temp + mon_precip + mon_AHMI + fdsi_annual, data = sample_data, distance = "bray")

# Select model
model_select <- ordiR2step(null, scope = formula(full_model), R2scope = TRUE, pstep = 0.05)

formula(model_select)
print(model_select)
summary(model_select)

# Check for VIF
vif_values <- vif.cca(model_select)
print(vif_values)

marginal_effects <- anova(model_select, by = "margin")
marginal_effects

anova_overall <- anova.cca(model_select, permutations = 999)
print(anova_overall)
db_rda_scores <- scores(model_select, display = "sites")
db_rda_env <- scores(model_select, display = "bp")
scores <- as.data.frame(db_rda_env)

# Combine the scores with the sample data
sites_data <- data.frame(db_rda_scores, sample_data) 
env_data <- data.frame(db_rda_env)

location_colours <- c("#EAA48F","#92c9c4")

rownames(env_data)
env_data$name <- c("Spring","Summer","Winter", "Watchmaker", "Quarterly MINT", "Quarterly Precip", "Monthly AHMI")


p <- ggplot(sites_data, aes(x = CAP1, y = CAP2, color = Location)) +
  geom_point(aes(shape = Season), size = 3) +
  theme_minimal() 


find_hull <- function(df) df[chull(df$CAP1, df$CAP2), ]

# Calculate polygon for each site
hulls <- sites_data %>%
  group_by(Location) %>%
  do(find_hull(.))

# Add polygons to the plot
p <- p + 
  geom_polygon(data = hulls, aes(x = CAP1, y = CAP2, fill = Location, group = Location), alpha = 0.02) 
wm_bb_ecm_dbrda <-p +
  scale_fill_manual(values = location_colours) + scale_color_manual(values = location_colours)+
  geom_segment(data = env_data, aes(x = 0, xend = CAP1, y = 0, yend = CAP2),
               arrow = arrow(length = unit(0.2, "cm")), color = "#000000", size = 0.8) +
  geom_text_repel(data = env_data, aes(x = CAP1, y = CAP2, label = name), 
                  color = "#000000", hjust =0, vjust = 0.0 ,force_pull = -0.06,
                  segment.size = 0.5,  # Controls the thickness of the lines
                  segment.linetype = "dashed", size = 3.7) + 
  theme_minimal() +
  theme(legend.position = "top") +
  labs(title = "db-RDA for ECM communities in LFP scats at Watchmaker and Bellbird", x = "CAP1", y = "CAP2") 
wm_bb_ecm_dbrda




## 2.2. Truffle-like Ectomycorrhizal fungi ----

ecm.list=row.names(guild_table)[which(guild_table$primary_lifestyle=="ectomycorrhizal")]
data=prune_taxa(ecm.list, ITSrel.count2)
hypo.list=row.names(guild_table)[which(guild_table$Fruitbody_type=="gasteroid-hypogeous")]
data=prune_taxa(hypo.list, data)

data <- prune_samples(sample_data(data)$Year >= 1998 & sample_data(data)$Year <= 2005, data) # Select years where samples were collected at both Bellbird & Watchmaker
data_truffle <- subset_samples(data, Location %in% "Bellbird" | Location %in% "Watchmaker" )


data_truffle = prune_taxa(taxa_sums(data_truffle) > 0, data_truffle) # Remove OTUs with zero count
truffle_ecm <- prune_samples(sample_sums(data_truffle) > 0, data_truffle) # Remove samples with zero OTUs 

### 2.2.1 Ordinations  ----


iDist <- distance((truffle_ecm), method="bray")
iMDS  <- ordinate(truffle_ecm, "NMDS", distance=iDist)

data.scores <- as.data.frame(scores(iMDS))
sort(data.scores$NMDS1) 

iDist <- distance((truffle_ecm), method="bray")
iMDS  <- ordinate(truffle_ecm, "NMDS", distance=iDist)

# Ordinate samples by Location
ordination_df <- as.data.frame(scores(iMDS, display = "sites"))
ord_truffle_ecm <- cbind(ordination_df, sample_data(truffle_ecm))

location_colours  <- c("#EAA48F","#92c9c4")


ordplot <- ggplot(ord_truffle_ecm, aes(x = NMDS1, y = NMDS2, color = Location)) +
  geom_point(size = 3) +
  scale_colour_manual(values = location_colours)+
  theme_light() +
  theme(aspect.ratio = 1) +
  ggtitle("Spatial Differences in Community Composition of Truffle-like ECM fungi in LFP scats collected at Watchmaker and Bellbird") +
  theme(plot.title = element_text(size = 16))
ordplot

find_hull <- function(df) df[chull(df$NMDS1, df$NMDS2), ]

# Calculate the polygons for each site
hulls <- ord_truffle_ecm %>%
  group_by(Location) %>%
  do(find_hull(.))

# Add polygons to ordination plot
ordplot + 
  geom_polygon(data = hulls, aes(x = NMDS1, y = NMDS2, fill = Location, group = Location), alpha = 0.05)  +
  theme(plot.title = element_text(size = 13)) +scale_fill_manual(values = location_colours) + theme_light()



### 2.2.2 PERMANOVAs  ----

sampledf <- data.frame(sample_data(truffle_ecm))
otu_table <- as(t(otu_table(truffle_ecm)), "matrix")
bray_dist <- vegdist(otu_table, method = "bray")

# PERMANOVA to test for differences in community composition between locations and seasons

permanova_result <- adonis2(bray_dist ~ Location, data = sampledf)
print(permanova_result)  

permanova_result <- adonis2(bray_dist ~ Location * Season, data = sampledf)
print(permanova_result) # Significant interaction between location and season 

# Test for difference in group dispersion between location
beta_location <- betadisper(bray_dist, sampledf$Location)
permutest(beta_location) # Dispersion does not differ significantly between locations
plot(beta_location, hull=FALSE, ellipse=TRUE) 
boxplot(beta_location)

# Test for difference in group dispersion between season
beta_season <- betadisper(bray_dist, sampledf$Season)
permutest(beta_season) # Dispersion does not differ significantly between seasons
plot(beta_season, hull=FALSE, ellipse=TRUE) 
boxplot(beta_season)



### 2.2.3 Indicspecies ----

# Test for genera that are indicative of particular sample groupings

sampledf <- data.frame(sample_data(truffle_ecm))
data_genus <- tax_glom(truffle_ecm, taxrank = "Genus") # Group taxa at the level of genus

otu <- as.data.frame(t(otu_table(data_genus)))
tax <- as.data.frame(t(tax_table(data_genus))) %>% rownames_to_column(var = "Taxa")
genus <- tax %>% filter(Taxa == "Genus")
genus <- genus[-1]
genus_otu <- as.data.frame(rbind(genus, otu))
colnames(genus_otu) <- genus_otu[1,]
genus_otu <- genus_otu[-1,]
genus_otu <- genus_otu %>%
  mutate(across(everything(), as.numeric))

genus_otu_bind <- cbind(sampledf, genus_otu)
Location <- (sample_data(genus_otu_bind)$Location) # Replace with $Season to test for genera that are indicative of differences between these groupings 


indicator_taxa <- multipatt(genus_otu, Location, func = "r", 
                            control = how(nperm=999)) 
summary(indicator_taxa)



### 2.2.4 DB-RDA ----

data_genus <- tax_glom(truffle_ecm, taxrank = "Genus") # Group taxa at the level of genus
otu_table <- as(t(otu_table(data_genus)), "matrix") 
sample_data <- data.frame(sample_data(data_genus))

col_names_indices <- data.frame(Index = seq_along(colnames(sample_data)), Column = colnames(sample_data))
print(col_names_indices)

sample_data[9:29] <- lapply(sample_data[9:29], function(x) as.numeric(as.character(x)))
sample_data[9:29]  <- scale(sample_data[9:29]) # Scale climate variables

## Identify highly correlated variables

clim_vars <- sample_data[, 9:29] # Select climate variables

cor_mat <- cor(clim_vars, use = "pairwise.complete.obs") # create the correlation matrix
high_cor <- findCorrelation(cor_mat, cutoff = 0.70, verbose = TRUE, names = TRUE) # Find high correlated pairs (correlation > 0.7)
non_cor_vars <- setdiff(names(clim_vars), high_cor) # List the climate variables that are not highly correlated

# create distance matrix
bray_distance <- vegdist(otu_table, method = "bray")

# use ordistep to chose model by permutation tests.
null<- capscale(bray_distance ~ 1, data = sample_data, distance = "bray")

# Prepare a model with all non correlated climate variables 

full_model <- capscale(bray_distance ~ Location + Season + week_precip + three_month_precip + three_mon_mint_av + twelve_mon_precip + twelve_mon_mint_av + twelve_mon_av_temp + mon_precip + mon_AHMI + fdsi_annual, data = sample_data, distance = "bray")

# Select model
model_select <- ordiR2step(null, scope = formula(full_model), R2scope = TRUE, pstep = 0.05)

formula(model_select)
print(model_select)
summary(model_select)

# Check for VIF
vif_values <- vif.cca(model_select)
print(vif_values)

marginal_effects <- anova(model_select, by = "margin")
marginal_effects

anova_overall <- anova.cca(model_select, permutations = 999)
print(anova_overall)
db_rda_scores <- scores(model_select, display = "sites")
db_rda_env <- scores(model_select, display = "bp")
scores <- as.data.frame(db_rda_env)

# Combine the scores with the sample data
sites_data <- data.frame(db_rda_scores, sample_data) 
env_data <- data.frame(db_rda_env)

location_colours  <- c("#EAA48F","#92c9c4")

rownames(env_data)
env_data$name <- c("Spring","Summer","Winter",  "Quarterly MINT","Watchmaker", "Annual MINT")


p <- ggplot(sites_data, aes(x = CAP1, y = CAP2, color = Location)) +
  geom_point(aes(shape = Season), size = 3) +
  theme_minimal() 


find_hull <- function(df) df[chull(df$CAP1, df$CAP2), ]

# Calculate polygon for each site
hulls <- sites_data %>%
  group_by(Location) %>%
  do(find_hull(.))

# Add polygons to the plot
p <- p + 
  geom_polygon(data = hulls, aes(x = CAP1, y = CAP2, fill = Location, group = Location), alpha = 0.02) 
wm_bb_ecm_dbrda_truffle <- p +
  scale_fill_manual(values = location_colours) + scale_color_manual(values = location_colours)+
  geom_segment(data = env_data, aes(x = 0, xend = CAP1, y = 0, yend = CAP2),
               arrow = arrow(length = unit(0.2, "cm")), color = "#000000", size = 0.8) +
  geom_text_repel(data = env_data, aes(x = CAP1, y = CAP2, label = name), 
                  color = "#000000", hjust =0, vjust = 0.0 ,force_pull = -0.06,
                  segment.size = 0.5,  # Controls the thickness of the lines
                  segment.linetype = "dashed", size = 3.7) + 
  theme_minimal() +
  theme(legend.position = "top") +
  labs(title = "db-RDA for Truffle-like ECM fungal communities in LFP scats at Watchmaker and Bellbird", x = "CAP1", y = "CAP2") 

# plot both ecm and truffle-like ecm DBRDA

wm_bb_ecm_dbrda
wm_bb_ecm_dbrda_truffle

wm_bb_ecm_dbrda_plot <- wm_bb_ecm_dbrda + 
  theme_light() +labs(title = NULL) + theme(legend.position = "none")

wm_bb_ecm_dbrda_truffle_plot <- wm_bb_ecm_dbrda_truffle + 
  theme_light() +labs(title = NULL) 

combined_plot <- ((wm_bb_ecm_dbrda_plot + wm_bb_ecm_dbrda_truffle_plot))   +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = "a")

combined_plot



# 3. Samples at Watchmaker to test LFP traits ----


## 3.1. Ectomycorrhizal fungi ----


data <- prune_samples(sample_data(ITSrel.count2)$Year >= 1998 & sample_data(ITSrel.count)$Year <= 2005, ITSrel.count)
data <- subset_samples(data, Location %in% "Watchmaker" )

ecm.list=row.names(guild_table)[which(guild_table$primary_lifestyle=="ectomycorrhizal")] # Filter for ECM taxa
data_ecm=prune_taxa(ecm.list, data)

data_ecm = prune_taxa(taxa_sums(data_ecm) > 0, data_ecm) # Remove OTUs with zero count
ecm <- prune_samples(sample_sums(data_ecm) > 0, data_ecm) # Remove samples with zero OTUs 


### 3.1.1. Ordinations ----


# Add trait data 

traits  <- read.csv('data/LFP_traits.csv',stringsAsFactors = FALSE) 

watchmaker_traits <- subset(traits, select = c(Sample, Sex, Mass))
ecm_df <- data.frame(sample_data(ecm))
ecm_df$Sample <- row.names(ecm_df)
watchmaker_traits$Sample

traits_sampledata <- ecm_df %>%   left_join(watchmaker_traits, by = c("Sample")) # Join physeq sample data with trait data

traits_sampledata <- traits_sampledata %>% filter(Mass != "NA") #Remove row with NA
ecm_cleaned <- prune_samples(sample_data(ecm)$agrf. %in% traits_sampledata$agrf., ecm) # remove rows in physeq that are not in traits_sampledata 

row.names(traits_sampledata) <- traits_sampledata$Sample #add row names back in order to merge with the physeq
sampledata <- sample_data(traits_sampledata)
ecm_trait = merge_phyloseq(ecm_cleaned, sampledata) # final physeq with LFP traits
 

iDist <- distance((ecm_trait), method="bray")
iMDS  <- ordinate(ecm_trait, "NMDS", distance=iDist)

data.scores <- as.data.frame(scores(iMDS))
sort(data.scores$NMDS1) 


# Ordinate samples by Location
ordination_df <- as.data.frame(scores(iMDS, display = "sites"))
ord_ecm_trait <- cbind(ordination_df, sample_data(ecm_trait))


ord_ecm_trait$Mass <- as.numeric(ord_ecm_trait$Mass)

ord_ecm_trait$Mass_bins <- (cut_number(ord_ecm_trait$Mass,3, labels = c("Low", "Medium", "High")))

ord_ecm_trait <- ord_ecm_trait %>% mutate(Sex = case_when(Sex == "F" ~ "Female",Sex == "M" ~ "Male"))
sex_colours <- c("#EAA48F","#92c9c4")



ordplot_sex <- ggplot(ord_ecm_trait, aes(x = NMDS1, y = NMDS2, color = Sex, shape = Mass_bins)) +
  geom_point(size = 3) +
  scale_colour_manual(values = sex_colours)+
  theme_light() +
  theme(aspect.ratio = 1) +
  ggtitle("Differences in Community Composition of ECM fungi in LFP Scats according to LFP Sex") +
  theme(plot.title = element_text(size = 16)) + scale_shape_discrete(name = "LFP Mass Class")
ordplot_sex

find_hull <- function(df) df[chull(df$NMDS1, df$NMDS2), ]

# Calculate the polygons for each site
hulls <- ord_ecm_trait %>%
  group_by(Sex) %>%
  do(find_hull(.))

# Add polygons to ordination plot
ordplot_sex <- ordplot_sex + 
  geom_polygon(data = hulls, aes(x = NMDS1, y = NMDS2, fill = Sex, group = Sex), alpha = 0.05)  +
  theme(plot.title = element_text(size = 13)) +scale_fill_manual(values = sex_colours) + theme_light()

ordplot_sex

ordplot_mass_ecm <- ggplot(ord_ecm_trait, aes(x = NMDS1, y = NMDS2, color = Mass_bins, shape = Sex)) +
  geom_point(size = 3) +
  theme_light() +
  theme(aspect.ratio = 1) +
  ggtitle("Differences in Community Composition of ECM fungi in LFP Scats according to LFP Mass") +
  theme(plot.title = element_text(size = 13)) + scale_color_manual(
    name = "LFP mass (g)",
    values = c("#C9E5E2", "#61888B", "#004148"))
ordplot_mass_ecm





## Plot relative abundance of ECM genera between males and females


colour_theme <- c("#C2CFD0" 
                 , "#92c9c4", "#90D8D2", "#809B97", "#95A58E","#4a8780", "#006973",  "#93715c",
                 "#9a8678", "grey",  "#E0D1C3","#D1B08F","#A67048", "#778773", "#4C2B16","#B27440",
                 "#CF866F", "#EAA48F", "#e3ad9a","#F8D8D2"
)
ps_transformed <- ecm_trait %>%  tax_glom(taxrank = "Genus") %>%    # Glom at the genus level
  transform_sample_counts(function(x) { x / sum(x) }) %>% # Transform to relative abundance
  psmelt()                                              

# Step 2: Rank the genera by total abundance across all Sexs
genus_ranks <- ps_transformed %>%
  group_by(Genus) %>%
  summarize(Total_Abundance = sum(Abundance, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(Total_Abundance)) 

genus_ranks <- genus_ranks%>%
  mutate(Rank = row_number()) %>%
  mutate(Genus = ifelse(Rank > 19, "Other", Genus)) 

ps_phylum <- ps_transformed %>%
  left_join(genus_ranks, by = "Genus") %>%
  mutate(Genus = ifelse(is.na(Rank), "Other", Genus)) %>%
  group_by(Sex, Genus) %>%
  summarize(Abundance = sum(Abundance, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(Sex, Genus) %>%
  summarize(Total_Abundance = sum(Abundance, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(Sex) %>%
  mutate(Relative_Abundance = Total_Abundance / sum(Total_Abundance, na.rm = TRUE)) %>%
  ungroup() 


ps_phylum_all <- ps_phylum

ps_phylum_all <- ps_phylum_all %>% mutate(Sex = case_when(Sex == "F" ~ "Female",Sex == "M" ~ "Male"))
sex_wm <-  ggplot(ps_phylum_all, aes(x = Sex, y = Relative_Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colour_theme) +
  labs(x = "LFP Sex",
       y = "Relative Abundance",
       title = NULL) +
  scale_y_continuous(labels = scales::percent) +        # Ensure y-axis is in percentage
  theme_light()

# Plot with ordination

sex_wm_plot <- sex_wm + 
  theme_light() +labs(title = NULL) + theme(axis.title.x = element_blank())

ordplot_ecm_wm_plot <- ordplot_sex + 
  theme_light() +labs(title = NULL) 

combined_plot <- ((ordplot_ecm_wm_plot + sex_wm_plot))   +
  plot_layout(widths = c(1, 0.5)) +
  plot_annotation(tag_levels = "a")

combined_plot




### 3.1.2 PERMANOVAs  ----

sampledf <- data.frame(sample_data(ecm_trait))
otu_table <- as(t(otu_table(ecm_trait)), "matrix")
bray_dist <- vegdist(otu_table, method = "bray")

# PERMANOVA to test for differences in community composition between mass and sex

permanova_result <- adonis2(bray_dist ~ Sex, data = sampledf)
print(permanova_result)   # Composition differs significantly between Male and Female LFP

permanova_result <- adonis2(bray_dist ~ Mass, data = sampledf)
print(permanova_result)  # Composition differs significantly between mass of LFP

permanova_result <- adonis2(bray_dist ~ Sex * Mass, data = sampledf)
print(permanova_result) # No interaction between Sex and Mass 

permanova_result <- adonis2(bray_dist ~ Sex + Season, data = sampledf)
print(permanova_result) # No interaction between Sex and Mass 

# Test for difference in group dispersion between Sex
beta_location <- betadisper(bray_dist, sampledf$Location)
permutest(beta_location) # Dispersion does not differ significantly between locations
plot(beta_location, hull=FALSE, ellipse=TRUE) 
boxplot(beta_location)

# Test for difference in group dispersion between season
beta_season <- betadisper(bray_dist, sampledf$Season)
permutest(beta_season) # Dispersion does not differ significantly between seasons
plot(beta_season, hull=FALSE, ellipse=TRUE) 
boxplot(beta_season)



### 3.1.3 Indicspecies ----

# Test for genera that are indicative of particular sample groupings

sampledf <- data.frame(sample_data(ecm_trait))
data_genus <- tax_glom(ecm_trait, taxrank = "Genus") # Group taxa at the level of genus

otu <- as.data.frame(t(otu_table(data_genus)))
tax <- as.data.frame(t(tax_table(data_genus))) %>% rownames_to_column(var = "Taxa")
genus <- tax %>% filter(Taxa == "Genus")
genus <- genus[-1]
genus_otu <- as.data.frame(rbind(genus, otu))
colnames(genus_otu) <- genus_otu[1,]
genus_otu <- genus_otu[-1,]
genus_otu <- genus_otu %>%
  mutate(across(everything(), as.numeric))


# Sex
genus_otu_bind <- cbind(sampledf, genus_otu)
Sex <- (sample_data(genus_otu_bind)$Sex) 

indicator_taxa <- multipatt(genus_otu, Sex, func = "r", 
                            control = how(nperm=999)) 
summary(indicator_taxa)

# Mass

# Break mass into 3 bins - low, medium, and high
sampledf$Mass <- as.numeric(sampledf$Mass)

sampledf$Mass_bins <- (cut_number(sampledf$Mass,3, labels = c("Low", "Medium", "High")))


genus_otu_bind <- cbind(sampledf, genus_otu)
Mass <- (sample_data(genus_otu_bind)$Mass_bins) 

indicator_taxa <- multipatt(genus_otu, Mass, func = "r", 
                            control = how(nperm=999)) 
summary(indicator_taxa)

# Visualise how mass is distributed between bins 

ggplot(sampledf, aes(x = Mass_bins, y = Mass, fill = Mass_bins)) +
  geom_boxplot() +scale_fill_manual(
    name = "LFP mass (g)",
    values = c("#C9E5E2", "#61888B", "#004148")
  ) +
  theme_light() +
  labs(
    x = "LFP Mass Category",
    y = "Body Mass (g)",
    title = "Distribution of LFP Mass across Mass Categories"
  ) 

sampledf_trait_ecm <- sampledf

## 3.2. Truffle-like Ectomycorrhizal fungi ----


data <- prune_samples(sample_data(ITSrel.count2)$Year >= 1998 & sample_data(ITSrel.count)$Year <= 2005, ITSrel.count)
data <- subset_samples(data, Location %in% "Watchmaker" )

ecm.list=row.names(guild_table)[which(guild_table$primary_lifestyle=="ectomycorrhizal")] # Filter for ECM taxa
data_ecm=prune_taxa(ecm.list, data)

hypo.list=row.names(guild_table)[which(guild_table$Fruitbody_type=="gasteroid-hypogeous")]
data_truffle=prune_taxa(hypo.list, data_ecm)

data_truffle = prune_taxa(taxa_sums(data_truffle) > 0, data_truffle) # Remove OTUs with zero count
truffle_ecm <- prune_samples(sample_sums(data_truffle) > 0, data_truffle) # Remove samples with zero OTUs 


### 3.2.1. Ordinations ----


# Add trait data 

traits  <- read.csv('data/LFP_traits.csv',stringsAsFactors = FALSE) 

watchmaker_traits <- subset(traits, select = c(Sample, Sex, Mass))
truffle_ecmdf <- data.frame(sample_data(truffle_ecm))
truffle_ecmdf$Sample <- row.names(truffle_ecmdf)
watchmaker_traits$Sample

traits_sampledata <- truffle_ecmdf %>%   left_join(watchmaker_traits, by = c("Sample")) # Join physeq sample data with trait data

traits_sampledata <- traits_sampledata %>% filter(Mass != "NA") #Remove row with NA
#Add binned categories for mass
traits_sampledata$Mass <- as.numeric(traits_sampledata$Mass)
traits_sampledata$Mass_bins <- (cut_number(traits_sampledata$Mass,3, labels = c("Low", "Medium", "High")))


truffle_ecm_cleaned <- prune_samples(sample_data(truffle_ecm)$agrf. %in% traits_sampledata$agrf., truffle_ecm) # remove rows in physeq that are not in traits_sampledata 

row.names(traits_sampledata) <- traits_sampledata$Sample #add row names back in order to merge with the physeq
sampledata <- sample_data(traits_sampledata)
truffle_ecm_trait = merge_phyloseq(truffle_ecm_cleaned, sampledata) # final physeq with LFP traits


iDist <- distance((truffle_ecm_trait), method="bray")
iMDS  <- ordinate(truffle_ecm_trait, "NMDS", distance=iDist)

data.scores <- as.data.frame(scores(iMDS))
sort(data.scores$NMDS1) 


# Ordinate samples by Location
ordination_df <- as.data.frame(scores(iMDS, display = "sites"))
ord_truffle_ecm_trait <- cbind(ordination_df, sample_data(truffle_ecm_trait))

ord_truffle_ecm_trait <- ord_truffle_ecm_trait %>% mutate(Sex = case_when(Sex == "F" ~ "Female",Sex == "M" ~ "Male"))
sex_colours <- c("#EAA48F","#92c9c4")


ordplot_sex <- ggplot(ord_truffle_ecm_trait, aes(x = NMDS1, y = NMDS2, color = Sex)) +
  geom_point(size = 3) +
  scale_colour_manual(values = sex_colours)+
  theme_light() +
  theme(aspect.ratio = 1) +
  ggtitle("Differences in Community Composition of Truffle-like ECM fungi in LFP Scats according to LFP Sex") +
  theme(plot.title = element_text(size = 16))
ordplot_sex

find_hull <- function(df) df[chull(df$NMDS1, df$NMDS2), ]

# Calculate the polygons for each site
hulls <- ord_truffle_ecm_trait %>%
  group_by(Sex) %>%
  do(find_hull(.))

# Add polygons to ordination plot
ordplot_sex + 
  geom_polygon(data = hulls, aes(x = NMDS1, y = NMDS2, fill = Sex, group = Sex), alpha = 0.05)  +
  theme(plot.title = element_text(size = 13)) +scale_fill_manual(values = sex_colours) + theme_light()


ordplot_mass_truffle_ecm <- ggplot(ord_truffle_ecm_trait, aes(x = NMDS1, y = NMDS2, color = Mass_bins, shape = Sex)) +
  geom_point(size = 3) +
  theme_light() +
  theme(aspect.ratio = 1) +
  ggtitle("Differences in Community Composition of ECM fungi in LFP Scats according to LFP Mass") +
  theme(plot.title = element_text(size = 13)) + scale_color_manual(
    name = "LFP Mass Class",
    values = c("#C9E5E2", "#61888B", "#004148")) + scale_shape_discrete(name = "LFP Sex")

ordplot_mass_truffle_ecm

## Plot ecm and truffle-like ecm mass and sex ordinations
ordplot_mass_ecm
ordplot_mass_truffle_ecm


combined_plot <- (ordplot_mass_ecm + ggtitle(NULL) +theme(legend.position = "none")) +
  (ordplot_mass_truffle_ecm + ggtitle(NULL)) +  plot_annotation(tag_levels = "a")

combined_plot




 ### 3.2.2 PERMANOVAs  ----

sampledf <- data.frame(sample_data(truffle_ecm_trait))
otu_table <- as(t(otu_table(truffle_ecm_trait)), "matrix")
bray_dist <- vegdist(otu_table, method = "bray")

# PERMANOVA to test for differences in community composition between mass and sex

permanova_result <- adonis2(bray_dist ~ Sex, data = sampledf)
print(permanova_result)   # Composition does not differ significantly between Male and Female LFP

permanova_result <- adonis2(bray_dist ~ Mass, data = sampledf)
print(permanova_result)  # Composition differs significantly between mass of LFP

permanova_result <- adonis2(bray_dist ~ Sex * Mass, data = sampledf)
print(permanova_result) # No interaction between Sex and Mass 

permanova_result <- adonis2(bray_dist ~ Mass + Season, data = sampledf)
print(permanova_result) #Season explains the most variation

# Test for difference in group dispersion between Sex
beta_location <- betadisper(bray_dist, sampledf$Location)
permutest(beta_location) # Dispersion does not differ significantly between locations
plot(beta_location, hull=FALSE, ellipse=TRUE) 
boxplot(beta_location)

# Test for difference in group dispersion between season
beta_season <- betadisper(bray_dist, sampledf$Season)
permutest(beta_season) # Dispersion does not differ significantly between seasons
plot(beta_season, hull=FALSE, ellipse=TRUE) 
boxplot(beta_season)



### 3.2.3 Indicspecies ----

# Test for genera that are indicative of particular sample groupings

sampledf <- data.frame(sample_data(truffle_ecm_trait))
data_genus <- tax_glom(truffle_ecm_trait, taxrank = "Genus") # Group taxa at the level of genus

otu <- as.data.frame(t(otu_table(data_genus)))
tax <- as.data.frame(t(tax_table(data_genus))) %>% rownames_to_column(var = "Taxa")
genus <- tax %>% filter(Taxa == "Genus")
genus <- genus[-1]
genus_otu <- as.data.frame(rbind(genus, otu))
colnames(genus_otu) <- genus_otu[1,]
genus_otu <- genus_otu[-1,]
genus_otu <- genus_otu %>%
  mutate(across(everything(), as.numeric))


# Sex
genus_otu_bind <- cbind(sampledf, genus_otu)
Sex <- (sample_data(genus_otu_bind)$Sex) 

indicator_taxa <- multipatt(genus_otu, Sex, func = "r", 
                            control = how(nperm=999)) 
summary(indicator_taxa)

# Mass

genus_otu_bind <- cbind(sampledf, genus_otu)
Mass <- (sample_data(genus_otu_bind)$Mass_bins) 

indicator_taxa <- multipatt(genus_otu, Mass, func = "r", 
                            control = how(nperm=999)) 
summary(indicator_taxa)


# Visualise how mass is distributed between bins 






ggplot(sampledf, aes(x = Mass_bins, y = Mass, fill = Mass_bins)) +
  geom_boxplot() +scale_fill_manual(
    name = "LFP mass (g)",
    values = c("#C9E5E2", "#61888B", "#004148")
  ) +
  theme_light() +
  labs(
    x = "LFP Mass Category",
    y = "Body Mass (g)",
    title = "Distribution of LFP Mass across Mass Categories"
  ) 


sampledf_trait_truffle_ecm <- sampledf

# Plot with the all ECM data.
# sampledf_trait_ecm and sampledf_trait_truffle_ecm

all_datasets <- bind_rows(
  sampledf_trait_ecm   %>% mutate(Dataset = "All ECM Fungi"),
  sampledf_trait_truffle_ecm %>% mutate(Dataset = "Truffle-like ECM Fungi")
)

rary(dplyr)

# summarise max and min of each group to plot on figure
range <- all_datasets %>%
  filter(!is.na(Mass)) %>%
  group_by(Mass_bins, Dataset) %>%
  summarise(ymin = min(Mass), ymax = max(Mass))

ggplot(all_datasets, aes(x = Mass_bins, y = Mass, fill = Dataset)) +
  geom_boxplot(width = 0.6,
               position = position_dodge(width = 0.7),
               color = "black") + 
  geom_text(data = range, aes(x = Mass_bins, y = ymin, label = ymin, color = Dataset),
    position = position_dodge(width = 0.7),vjust = 1.6, size = 3) +
  geom_text(data = range, aes(x = Mass_bins, y = ymax, label = ymax, color = Dataset),
            position = position_dodge(width = 0.7),vjust = -0.5, size = 3 ) +
  scale_y_continuous(limits = c(1200, 2200)) + labs(x = "LFP Mass Categories", y ="LFP Mass (g)")+
  scale_fill_manual(values = c("All ECM Fungi" = "#2F7D78","Truffle-like ECM Fungi" = "#C9E5E2"),name = "Dataset") +
  scale_color_manual(values = c("All ECM Fungi" = "#004148","Truffle-like ECM Fungi" = "#2F7D78"), guide = "none") +
theme_light()


# Plot presence of genera between Male and Female LFPs
genus <- tax_glom(truffle_ecm_trait, taxrank = "Genus") # glom at level of Genus
genus_pa <- transform_sample_counts(genus, function(x) ifelse(x > 0, 1, 0)) # convert to presence-absence

sex_presence <- merge_samples(genus_pa, "Sex")
sex_presence <- transform_sample_counts(sex_presence, function(x) ifelse(x > 0, 1, 0))


otu <- as.data.frame(t(otu_table(sex_presence)))
tax <- as.data.frame(tax_table(sex_presence))

otu$Genus <- tax$Genus

genus_long <- pivot_longer(otu, cols = c("F", "M"), 
                           names_to = "Sex", values_to = "Presence")

genus_long<- genus_long %>%
  mutate(Sex = case_when(Sex == "F" ~ "Female",Sex == "M"~ "Male" , TRUE ~ Sex))

genus_long$Genus <- factor(genus_long$Genus, levels = rev(sort(unique(genus_long$Genus)))) #sorts alphabetically

ggplot(genus_long, aes(x = Sex, y = Genus, fill = factor(Presence))) +
  geom_tile(color = "black")  +
  scale_fill_manual(values = c("0" = "white", "1" = "#778773"), name = "Presence", labels = c("Absent", "Present"))+
  theme_minimal() +
  labs(title = "Presence of Truffle-like ECM genera between Male and Female LFPs",
       x = "Sex",
       y = NULL) +
  theme(axis.text.y = element_text(size = 9), legend.title = element_blank(),axis.title.x = element_text(size = 10))



library(MicEco)
genus <- tax_glom(truffle_ecm_trait, taxrank = "Genus") # glom at level of Genus
ps_euler(genus, group = "Sex", shape = "ellipse", plot = TRUE, quantities = list(type=c("counts")),  col = "#7C644C")




# 4. Samples at Bellbird to test topography ----


## 4.1. Ectomycorrhizal fungi ----

data <- subset_samples(ITSrel.count2, Location %in% "Bellbird" )

ecm.list=row.names(guild_table)[which(guild_table$primary_lifestyle=="ectomycorrhizal")] # Filter for ECM taxa
data_ecm=prune_taxa(ecm.list, data)

data_ecm = prune_taxa(taxa_sums(data_ecm) > 0, data_ecm) # Remove OTUs with zero count
ecm <- prune_samples(sample_sums(data_ecm) > 0, data_ecm) # Remove samples with zero OTUs 


### 4.1.1. Ordinations ----


# Add topography data  

topography <- read.csv('data/topography.csv') # Read in data with topographic position of LFP scat collection

ecm_sampledf <- data.frame(sample_data(ecm))



topo_sampledata <- ecm_sampledf %>%   left_join(topography, by = c("Sample")) # Join physeq sample data with topography data

topo_sampledata <- topo_sampledata %>% filter(TWI != 0) # Remove records missing data

# convert aspect to eastness and northness

topo_sampledata$eastness <- sin(topo_sampledata$Aspect * (pi/180))
topo_sampledata$northness <- cos(topo_sampledata$Aspect * (pi/180))




ecm_cleaned <- prune_samples(sample_data(ecm)$agrf. %in% topo_sampledata$agrf., ecm) # remove rows in physeq that are not in topo_sampledata 

row.names(topo_sampledata) <- topo_sampledata$Sample #add row names back in order to merge with the physeq
sampledata <- sample_data(topo_sampledata)
ecm_topo = merge_phyloseq(ecm_cleaned, sampledata) # final physeq with LFP traits


iDist <- distance((ecm_topo), method="bray")
iMDS  <- ordinate(ecm_topo, "NMDS", distance=iDist)

data.scores <- as.data.frame(scores(iMDS))
sort(data.scores$NMDS1) 

ecm_topo <- subset_samples(ecm_topo, !agrf. %in% c("1")) # Remove outliers that obscure patterns in ordination

iDist <- distance((ecm_topo), method="bray")
iMDS  <- ordinate(ecm_topo, "NMDS", distance=iDist)

# Ordinate samples by Location
ordination_df <- as.data.frame(scores(iMDS, display = "sites"))
ord_ecm_topo <- cbind(ordination_df, sample_data(ecm_topo))

ordplot_topo <- ggplot(ord_ecm_topo, aes(x = NMDS1, y = NMDS2, color = Aspect)) + #Change to visualise TWI/Slope/Aspect
  geom_point(size = 3) +
  theme_light() +
  theme(aspect.ratio = 1) +
  ggtitle("Differences in Community Composition of ECM fungi in LFP Scats across a topographic gradient") +
  theme(plot.title = element_text(size = 13)) + scale_color_gradient(low="#C9E5E2", high="#004148")
ordplot_topo


### 4.1.2 PERMANOVAs  ----

sampledf <- data.frame(sample_data(ecm_topo))
otu_table <- as(t(otu_table(ecm_topo)), "matrix")
bray_dist <- vegdist(otu_table, method = "bray")

# PERMANOVA to test for differences in community composition across topography (change variable to test TWI, Slope, Aspect)

permanova_result <- adonis2(bray_dist ~ northness, data = sampledf)
print(permanova_result)   # Composition does not differ significantly across topographic gradient






## 4.2. Truffle-like Ectomycorrhizal fungi ----

data <- subset_samples(ITSrel.count2, Location %in% "Bellbird" )

ecm.list=row.names(guild_table)[which(guild_table$primary_lifestyle=="ectomycorrhizal")] # Filter for ECM taxa
data_ecm=prune_taxa(ecm.list, data)

hypo.list=row.names(guild_table)[which(guild_table$Fruitbody_type=="gasteroid-hypogeous")]
data_truffle=prune_taxa(hypo.list, data_ecm)

data_truffle = prune_taxa(taxa_sums(data_truffle) > 0, data_truffle) # Remove OTUs with zero count
truffle_ecm <- prune_samples(sample_sums(data_truffle) > 0, data_truffle) # Remove samples with zero OTUs 


### 4.2.1. Ordinations ----


# Add topography data  

topography <- read.csv('data/topography.csv') # Read in data with topographic position of LFP scat collection

sampledata <- data.frame(sample_data(truffle_ecm))
sampledata$Sample


topo_sampledata <- sampledata %>%   left_join(topography, by = c("Sample")) # Join physeq sample data with topography data

topo_sampledata <- topo_sampledata %>% filter(TWI != 0) # Remove records missing data


# convert aspect to eastness and northness

topo_sampledata$eastness <- sin(topo_sampledata$Aspect * (pi/180))
topo_sampledata$northness <- cos(topo_sampledata$Aspect * (pi/180))


truffle_ecm_cleaned <- prune_samples(sample_data(truffle_ecm)$agrf. %in% topo_sampledata$agrf., truffle_ecm) # remove rows in physeq that are not in topo_sampledata 

row.names(topo_sampledata) <- topo_sampledata$Sample #add row names back in order to merge with the physeq
sampledata <- sample_data(topo_sampledata)
truffle_ecm_topo = merge_phyloseq(truffle_ecm_cleaned, sampledata) # final physeq with LFP traits



iDist <- distance((truffle_ecm_topo), method="bray")
iMDS  <- ordinate(truffle_ecm_topo, "NMDS", distance=iDist)

data.scores <- as.data.frame(scores(iMDS))
sort(data.scores$NMDS1) 

truffle_ecm_topo <- subset_samples(truffle_ecm_topo, !agrf. %in% c("1","18","286")) # Remove outliers that obscure patterns in ordination


iDist <- distance((truffle_ecm_topo), method="bray")
iMDS  <- ordinate(truffle_ecm_topo, "NMDS", distance=iDist)

# Ordinate samples by Location
ordination_df <- as.data.frame(scores(iMDS, display = "sites"))
ord_truffle_ecm_topo <- cbind(ordination_df, sample_data(truffle_ecm_topo))

ordplot_topo <- ggplot(ord_truffle_ecm_topo, aes(x = NMDS1, y = NMDS2, color = TWI)) + #Change to visualise TWI/Slope/Aspect
  geom_point(size = 3) +
  theme_light() +
  theme(aspect.ratio = 1) +
  ggtitle("Differences in Community Composition of Truffle-like ECM fungi in LFP Scats across a topographic gradient") +
  theme(plot.title = element_text(size = 13)) + scale_color_gradient(low="#C9E5E2", high="#004148")
ordplot_topo


### 4.2.2 PERMANOVAs  ----

sampledf <- data.frame(sample_data(truffle_ecm_topo))
otu_table <- as(t(otu_table(truffle_ecm_topo)), "matrix")
bray_dist <- vegdist(otu_table, method = "bray")

# PERMANOVA to test for differences in community composition across topography (change variable to test TWI, Slope, Aspect)


permanova_result <- adonis2(bray_dist ~ northness , data = sampledf)
print(permanova_result)   # Composition does not differ differ significantly across topographic gradient






