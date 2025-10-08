# Purpose: Statistical analysis of differences in truffle-like fungal richness between seasons, geography, and LFP traits

# Load packages 


library(phyloseq)
library(ggplot2)
library(vegan)
library(ggpubr)
library(patchwork)
library(pscl)
library(caret)
library(MASS)
library(MuMIn)
library(car)
library(emmeans)
library(glmmTMB)

# 0. Visualise sample collections ----

# 3 datasets:
## 1: Samples from 5 sites between 1993-95
## 2: Samples from Watchmaker and Bellbird between 1998-2005
## 3: Samples from Bellbird with data on elevation of scat collection

# Filter for each dataset and combine into one physeq object 

data <- prune_samples(sample_data(ITSrel.count2)$Year < 1996, ITSrel.count2)
data <- prune_samples(!(sample_data(data)$Location %in% c("Ellery Two","Puggaree Rd")), data) # Remove sites with < 5 records
dataset1 <- data


data <- prune_samples(sample_data(ITSrel.count2)$Year >= 1998 & sample_data(ITSrel.count2)$Year <= 2005, ITSrel.count2)
data <- subset_samples(data, Location %in% "Bellbird" | Location %in% "Watchmaker" )
dataset2 <- data

dataset_1and2 <- merge_phyloseq(dataset1, dataset2) #merge dataset 1 and 2


topography <- read.csv('data/topography.csv')
data <- subset_samples(ITSrel.count2, Location %in% "Bellbird")
sample_data_df <- data.frame(sample_data(data))
merged_df <- merge(sample_data_df, topography, by = "Sample") # Merge physeq sample data with data on topographic position of sample collection
row.names(merged_df) <- merged_df$Sample
sample_data(data) <- merged_df
non_zero <- merged_df$Elevation != 0
data <- prune_samples(non_zero, data)
dataset3 <- data

# remove duplicates of dataset3 from dataset_1and2 before merging

sample_ID_dataset1and2 <- sample_data(dataset_1and2)$textbefore
sample_ID_dataset3 <- sample_data(dataset3)$textbefore
duplicates <- intersect(sample_ID_dataset1and2, sample_ID_dataset3)

dataset3_pruned <- prune_samples(!(sample_data(dataset3)$textbefore %in% duplicates), dataset3) # removes duplicates 

# Merge dataset 1&2 with dataset3
all_data <- merge_phyloseq(dataset_1and2, dataset3_pruned)

# Merge the Location columns
sample_data <- data.frame(sample_data(all_data))
sample_data$Sample <- rownames(sample_data)
sample_data <- sample_data %>%
  mutate(Location = coalesce(Location, Location.x, Location.y))
rownames(sample_data) <- sample_data$Sample
sample_data(all_data) <- sample_data(sample_data)

all_data <- prune_taxa(taxa_sums(all_data) > 0, all_data) ## remove OTUs with zero count 
all_data<- prune_samples(sample_sums(all_data) > 0, all_data) ## remove samples with zero OTUs

location_colours <- c("#AFD1CE", "#95AA90", "#D1B08F", "#CF866F", "#9D7963")

sample_data_df <- data.frame(sample_data(all_data))

# group data to sum number of samples per site per year
location_year_count <- sample_data_df %>%
  group_by(Location, Year) %>% #
  summarize(Sample_Count = n()) %>%
  ungroup()

location_year_count$Year <- as.factor(location_year_count$Year)

ggplot(location_year_count, aes(x = Year, y = Sample_Count, fill = Location)) +
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = location_colours,name = "Scat collection site") +
  labs(#title = "Total Number of Scat Samples Each Year",
    x = "Year",
    y = "Number of Samples") + 
  theme_light() +
  scale_x_discrete(breaks = c(1993,1995,1997,1999,2001,2003,2005,2007,2009,2011,2013)) + 
  theme(    legend.key.size = unit(0.5, "cm"))


# 1. Samples at all sites between 1993-1995 ----


data <- prune_samples(sample_data(ITSrel.count)$Year < 1996, ITSrel.count)
data <- prune_samples(!(sample_data(data)$Location %in% c("Ellery Two","Puggaree Rd")), data) # Remove sites with < 5 samples


## 1.1. Ectomycorrhizal fungi ----

ecm.list=row.names(guild_table)[which(guild_table$primary_lifestyle=="ectomycorrhizal")] # Filter for ECM taxa
data_ecm=prune_taxa(ecm.list, data)


data_ecm = prune_taxa(taxa_sums(data_ecm) > 0, data_ecm) # Remove OTUs with zero count
ecm <- prune_samples(sample_sums(data_ecm) > 0, data_ecm) # Remove samples with zero OTUs 


otu_table_df <- as.data.frame(otu_table(ecm))
otu_table_df <- t(otu_table_df)

spec_accum <- specaccum(otu_table_df)
spec_accum_df <- data.frame(Samples = spec_accum$sites,Richness = spec_accum$richness,SD = spec_accum$sd)

species_accumulation <- ggplot(spec_accum_df, aes(x = Samples, y = Richness)) +
  geom_line() +
  geom_ribbon(aes(ymin = Richness - SD, ymax = Richness + SD), alpha = 0.2) +
  labs(title = "ECM OTU Accumulation Curve",
       x = "Number of Samples",
       y = "OTU Richness") +
  theme_light()

species_accumulation


### 1.1.1. Calculate richness ----


meta = data.frame(sample_data(ecm))
filtered_reads <- sample_sums((ecm))
meta <- cbind(meta, filtered_reads)
diversity <-estimate_richness(ecm, measures=c("Observed", "Chao1", "Shannon", "Simpson"))
alphadiv_ecm <- cbind(diversity, meta)

### 1.1.2. Richness between seasons ----

# Negative binomial GLM 
nb_model_season <- glm.nb(Observed ~ Season , data = alphadiv_ecm)
summary(nb_model_season)

# Test differences with emmeans
emm <- emmeans(nb_model_season, ~ Season)
pairwise_results <- pairs(emm)
plot(emm)


# Plot seasonal differences

## Plot differences in richness

alphadiv_ecm$Season <- factor(alphadiv_ecm$Season, 
                                levels = c("Summer", "Autumn", "Winter", "Spring"))


season_colors <- c("Summer"= "#EAA48F","Autumn" = "#D1B08F","Winter" = "#92c9c4","Spring" = "#95A58E")                   

ecm_richness_season_plot <-ggplot(alphadiv_ecm, aes(x = Season, y = Observed)) +
  geom_boxplot(aes(fill  = Season)) +
  theme_minimal() +
  labs(title = "Seasonal Differences in Diversity of ECM ",
       x = "Season",
       y = "Observed Diversity") +
  stat_compare_means(comparisons = list(c("Winter","Spring")), method = "t.test", label = "p.signif") + # Based on emmeans results
  theme(
    axis.title.y = element_text(size = 13),
    axis.title.x = element_text(size = 13),
    text = element_text(size = 15), legend.position = "none") +
  scale_fill_manual(values = season_colors)
ecm_richness_season_plot


### 1.1.3. Richness between locations ----

# Negative binomial GLM - Regional differences
nb_model_region <- glm.nb(Observed ~ Region , data = alphadiv_ecm)
summary(nb_model_region)

# Negative binomial GLM - Site differences
nb_model_site <- glm.nb(Observed ~ Location , data = alphadiv_ecm)
summary(nb_model_site)



### 1.1.4. GLMs ----

## Identify highly correlated variables

clim_vars <- alphadiv_ecm[, 14:34] # Select climate variables

cor_mat <- cor(clim_vars, use = "pairwise.complete.obs") # create the correlation matrix
high_cor <- findCorrelation(cor_mat, cutoff = 0.70, verbose = TRUE, names = TRUE) # Find high correlated pairs (correlation > 0.7)
non_cor_vars <- setdiff(names(clim_vars), high_cor) # List the climate variables that are not highly correlated

# Select best model

# Scale climate variables

alphadiv_ecm <- alphadiv_ecm %>% mutate_at(vars(14:34), as.numeric)
alphadiv_ecm[, c(14:34)] <- scale(alphadiv_ecm[, c(14:34)])

# Use dredge to identify non-correlated climate variables that best predict observed richness of truffle-like ECM, using negative binomial regression


model <-glm.nb(Observed ~ week_precip + three_month_precip + three_mon_AHMI + twelve_mon_precip + twelve_mon_mint_av + wk_av_temp + mon_precip + mon_AHMI, data = alphadiv_ecm, na.action = "na.fail")


glm_nb <-dredge(model, beta= "sd", evaluate = TRUE, rank = "AIC",
                   m.lim = c(0,2),fixed = NULL ,    # Maximum of two terms 
                   extra=c("adjR^2", "BIC", "Cp"),
                   trace=TRUE)
View(glm_nb)

# Best model (lowest AIC)
nb_model <- glm.nb(Observed ~  twelve_mon_precip + week_precip, data = alphadiv_ecm)
summary(nb_model) #Variables are not significant

# Compare to null model
nb_model_null<- glm.nb(Observed ~  1 , data = alphadiv_ecm)
summary(nb_model_null)
AIC(nb_model,nb_model_null)
anova(nb_model,nb_model_null) # Selected model not significantly better than the null model




## 1.2. Truffle-like Ectomycorrhizal fungi ----

ecm.list=row.names(guild_table)[which(guild_table$primary_lifestyle=="ectomycorrhizal")] # Filter for ECM taxa
data_ecm=prune_taxa(ecm.list, data)


hypo.list=row.names(guild_table)[which(guild_table$Fruitbody_type=="gasteroid-hypogeous")]
data_truffle=prune_taxa(hypo.list, data_ecm)

data_truffle = prune_taxa(taxa_sums(data_truffle) > 0, data_truffle) # Remove OTUs with zero count
truffle_ecm <- prune_samples(sample_sums(data_truffle) > 0, data_truffle) # Remove samples with zero OTUs 


otu_table_df <- as.data.frame(otu_table(truffle_ecm))
otu_table_df <- t(otu_table_df)

spec_accum <- specaccum(otu_table_df)
spec_accum_df <- data.frame(Samples = spec_accum$sites,Richness = spec_accum$richness,SD = spec_accum$sd)

species_accumulation <- ggplot(spec_accum_df, aes(x = Samples, y = Richness)) +
  geom_line() +
  geom_ribbon(aes(ymin = Richness - SD, ymax = Richness + SD), alpha = 0.2) +
  labs(title = "Truffle-like ECM OTU Accumulation Curve",
       x = "Number of Samples",
       y = "OTU Richness") +
  theme_light()

species_accumulation


### 1.2.1. Calculate richness ----


meta = data.frame(sample_data(truffle_ecm))
filtered_reads <- sample_sums((truffle_ecm))
meta <- cbind(meta, filtered_reads)
diversity <-estimate_richness(truffle_ecm, measures=c("Observed", "Chao1", "Shannon", "Simpson"))
alphadiv_truffle_ecm <- cbind(diversity, meta)

### 1.2.2. Richness between seasons ----

# Negative binomial GLM 
nb_model_season <- glm.nb(Observed ~ Season , data = alphadiv_truffle_ecm)
summary(nb_model_season)


# Plot seasonal differences

## Plot differences in richness

alphadiv_truffle_ecm$Season <- factor(alphadiv_truffle_ecm$Season, 
                              levels = c("Summer", "Autumn", "Winter", "Spring"))


season_colors <- c("Summer"= "#EAA48F","Autumn" = "#D1B08F","Winter" = "#92c9c4","Spring" = "#95A58E")                   

truffle_ecm_richness_season_plot <-ggplot(alphadiv_truffle_ecm, aes(x = Season, y = Observed)) +
  geom_boxplot(aes(fill  = Season)) +
  theme_minimal() +
  labs(title = "Seasonal Differences in Diversity of Truffle-like ECM ",
       x = "Season",
       y = "Observed Diversity") +

  theme(
    axis.title.y = element_text(size = 13),
    axis.title.x = element_text(size = 13),
    text = element_text(size = 15), legend.position = "none") +
  scale_fill_manual(values = season_colors)
truffle_ecm_richness_season_plot


### 1.2.3. Richness between locations ----

# Negative binomial GLM - Regional differences
nb_model_region <- glm.nb(Observed ~ Region , data = alphadiv_truffle_ecm)
summary(nb_model_region)

# Negative binomial GLM - Site differences
nb_model_site <- glm.nb(Observed ~ Location , data = alphadiv_truffle_ecm)
summary(nb_model_site)



### 1.2.4. GLMs ----

## Identify highly correlated variables

clim_vars <- alphadiv_truffle_ecm[, 14:34] # Select climate variables

cor_mat <- cor(clim_vars, use = "pairwise.complete.obs") # create the correlation matrix
high_cor <- findCorrelation(cor_mat, cutoff = 0.70, verbose = TRUE, names = TRUE) # Find high correlated pairs (correlation > 0.7)
non_cor_vars <- setdiff(names(clim_vars), high_cor) # List the climate variables that are not highly correlated

# Select best model

# Scale climate variables

alphadiv_truffle_ecm <- alphadiv_truffle_ecm %>% mutate_at(vars(14:34), as.numeric)
alphadiv_truffle_ecm[, c(14:34)] <- scale(alphadiv_truffle_ecm[, c(14:34)])

# Use dredge to identify non-correlated climate variables that best predict observed richness of truffle-like truffle_ecm, using negative binomial regression


model <-glm.nb(Observed ~ week_precip + three_month_precip + three_mon_AHMI + twelve_mon_precip + twelve_mon_mint_av + wk_av_temp + mon_precip + mon_AHMI, data = alphadiv_truffle_ecm, na.action = "na.fail")


glm_nb <-dredge(model, beta= "sd", evaluate = TRUE, rank = "AIC",
                m.lim = c(0,2),fixed = NULL ,    # Maximum of two terms 
                extra=c("adjR^2", "BIC", "Cp"),
                trace=TRUE)
View(glm_nb)

# Best model (lowest AIC)
nb_model <- glm.nb(Observed ~  twelve_mon_precip + three_month_precip, data = alphadiv_truffle_ecm)
summary(nb_model) 

# Compare to null model
nb_model_null<- glm.nb(Observed ~  1 , data = alphadiv_truffle_ecm)
summary(nb_model_null)
AIC(nb_model,nb_model_null)
anova(nb_model,nb_model_null) # Selected model not significantly better than the null model





# 2. Samples at Watchmaker and Bellbird between 1998-2005 ----

data <- prune_samples(sample_data(ITSrel.count)$Year >= 1998 & sample_data(ITSrel.count)$Year <= 2005, ITSrel.count)
data <- subset_samples(data, Location %in% "Bellbird" | Location %in% "Watchmaker" )


## 2.1. Ectomycorrhizal fungi ----

ecm.list=row.names(guild_table)[which(guild_table$primary_lifestyle=="ectomycorrhizal")] # Filter for ECM taxa
data_ecm=prune_taxa(ecm.list, data)


data_ecm = prune_taxa(taxa_sums(data_ecm) > 0, data_ecm) # Remove OTUs with zero count
ecm <- prune_samples(sample_sums(data_ecm) > 0, data_ecm) # Remove samples with zero OTUs 


otu_table_df <- as.data.frame(otu_table(ecm))
otu_table_df <- t(otu_table_df)


### 2.1.1. Calculate richness ----


meta = data.frame(sample_data(ecm))
filtered_reads <- sample_sums((ecm))
meta <- cbind(meta, filtered_reads)
diversity <-estimate_richness(ecm, measures=c("Observed", "Chao1", "Shannon", "Simpson"))
alphadiv_ecm <- cbind(diversity, meta)

### 2.1.2. Richness between seasons ----

# Negative binomial GLM 
nb_model_season <- glm.nb(Observed ~ Season , data = alphadiv_ecm)
summary(nb_model_season)

# Test differences with emmeans
emm <- emmeans(nb_model_season, ~ Season)
pairwise_results <- pairs(emm)
plot(emm)



### 2.1.3. Richness between locations ----

# Negative binomial GLM - Site differences
nb_model_site <- glm.nb(Observed ~ Location , data = alphadiv_ecm)
summary(nb_model_site)




## 2.2. Truffle-like Ectomycorrhizal fungi ----

ecm.list=row.names(guild_table)[which(guild_table$primary_lifestyle=="ectomycorrhizal")] # Filter for ECM taxa
data_ecm=prune_taxa(ecm.list, data)


hypo.list=row.names(guild_table)[which(guild_table$Fruitbody_type=="gasteroid-hypogeous")]
data_truffle=prune_taxa(hypo.list, data_ecm)

data_truffle = prune_taxa(taxa_sums(data_truffle) > 0, data_truffle) # Remove OTUs with zero count
truffle_ecm <- prune_samples(sample_sums(data_truffle) > 0, data_truffle) # Remove samples with zero OTUs 



### 2.2.1. Calculate richness ----


meta = data.frame(sample_data(truffle_ecm))
filtered_reads <- sample_sums((truffle_ecm))
meta <- cbind(meta, filtered_reads)
diversity <-estimate_richness(truffle_ecm, measures=c("Observed", "Chao1", "Shannon", "Simpson"))
alphadiv_truffle_ecm <- cbind(diversity, meta)

### 2.2.2. Richness between seasons ----

# Negative binomial GLM 
nb_model_season <- glm.nb(Observed ~ Season , data = alphadiv_truffle_ecm)
summary(nb_model_season)


### 2.2.3. Richness between locations ----


# Negative binomial GLM - Site differences
nb_model_site <- glm.nb(Observed ~ Location , data = alphadiv_truffle_ecm)
summary(nb_model_site)



# 3. Samples at Watchmaker to test LFP traits ----


data <- prune_samples(sample_data(ITSrel.count)$Year >= 1998 & sample_data(ITSrel.count)$Year <= 2005, ITSrel.count)
data <- subset_samples(data, Location %in% "Watchmaker" )

## 3.1. Ectomycorrhizal fungi ----

ecm.list=row.names(guild_table)[which(guild_table$primary_lifestyle=="ectomycorrhizal")] # Filter for ECM taxa
data_ecm=prune_taxa(ecm.list, data)

data_ecm = prune_taxa(taxa_sums(data_ecm) > 0, data_ecm) # Remove OTUs with zero count
ecm <- prune_samples(sample_sums(data_ecm) > 0, data_ecm) # Remove samples with zero OTUs 


### 3.1.1. Calculate richness ----


meta = data.frame(sample_data(ecm))
filtered_reads <- sample_sums((ecm))
meta <- cbind(meta, filtered_reads)
diversity <-estimate_richness(ecm, measures=c("Observed", "Chao1", "Shannon", "Simpson"))
alphadiv_ecm <- cbind(diversity, meta)


traits  <- read.csv('data/LFP_traits.csv',stringsAsFactors = FALSE) 
traits <- subset(traits, select = c(Sample, Sex, Mass))

traits_watchmaker <- alphadiv_ecm %>%
  left_join(traits, by = c("Sample"))

na_positions <- is.na(traits_watchmaker)
print(na_positions)
na_indices <- which(na_positions, arr.ind = TRUE)
print(na_indices)
traits_watchmaker_ecm <- traits_watchmaker[-29,] #remove row with no data




### 3.1.2 Richness between sex and mass ----

# Negative binomial GLM for sex
nb_model_season <- glm.nb(Observed ~ Sex , data = traits_watchmaker_ecm)
summary(nb_model_season)

# Plot
ggplot(traits_watchmaker_ecm, aes(x = Sex, y = Observed, fill = Sex )) +
  geom_boxplot() + labs(title = "Differences in ECM richness between Male and Female LFPs", x = "LFP Sex", y = "ECM Richness")

# Negative binomial GLM for mass
nb_model_mass <- glm.nb(Observed ~ Mass, data = traits_watchmaker_ecm)
summary(nb_model_mass)

# Plot
ggplot(traits_watchmaker_ecm, aes(x = Mass, y = Observed )) +
  geom_point() + geom_smooth(method = "lm")+ labs(title = "Differences in ECM richness across LFP size", x = "LFP Mass", y = "ECM Richness")

nb_model_mass <- glm.nb(Observed ~ Mass + I(Mass^2), data = traits_watchmaker_ecm) #Appears unimodal on plot - try quadratic term
summary(nb_model_mass) #not significant



## 3.2. Truffle-like Ectomycorrhizal fungi ----

ecm.list=row.names(guild_table)[which(guild_table$primary_lifestyle=="ectomycorrhizal")] # Filter for ECM taxa
data_ecm=prune_taxa(ecm.list, data)

hypo.list=row.names(guild_table)[which(guild_table$Fruitbody_type=="gasteroid-hypogeous")]
data_truffle=prune_taxa(hypo.list, data_ecm)

data_truffle = prune_taxa(taxa_sums(data_truffle) > 0, data_truffle) # Remove OTUs with zero count
truffle_ecm <- prune_samples(sample_sums(data_truffle) > 0, data_truffle) # Remove samples with zero OTUs 


### 3.2.1. Calculate richness ----


meta = data.frame(sample_data(truffle_ecm))
filtered_reads <- sample_sums((truffle_ecm))
meta <- cbind(meta, filtered_reads)
diversity <-estimate_richness(truffle_ecm, measures=c("Observed", "Chao1", "Shannon", "Simpson"))
alphadiv_truffle_ecm <- cbind(diversity, meta)


traits  <- read.csv('data/LFP_traits.csv',stringsAsFactors = FALSE) 
traits <- subset(traits, select = c(Sample, Sex, Mass))

traits_watchmaker <- alphadiv_truffle_ecm %>%
  left_join(traits, by = c("Sample"))

na_positions <- is.na(traits_watchmaker)
print(na_positions)
na_indices <- which(na_positions, arr.ind = TRUE)
print(na_indices)
traits_watchmaker_truffle_ecm <- traits_watchmaker[-23,] #remove row with no data




### 3.2.2 Richness between sex and mass ----

# Negative binomial GLM for sex
nb_model_season <- glm.nb(Observed ~ Sex , data = traits_watchmaker_truffle_ecm)
summary(nb_model_season)

# Plot
ggplot(traits_watchmaker_truffle_ecm, aes(x = Sex, y = Observed, fill = Sex )) +
  geom_boxplot() + labs(title = "Differences in Truffle-like ECM fungal richness between Male and Female LFPs", x = "LFP Sex", y = "Truffle-like ECM fungal richness")

# Negative binomial GLM for mass
nb_model_mass <- glm.nb(Observed ~ Mass, data = traits_watchmaker_truffle_ecm)
summary(nb_model_mass)

# Plot
ggplot(traits_watchmaker_truffle_ecm, aes(x = Mass, y = Observed )) +
  geom_point() + geom_smooth(method = "lm")+ labs(title = "Differences in Truffle-like ECM fungal richness across LFP size", x = "LFP Mass", y = "Truffle-like ECM fungal richness")







# 4. Samples at Bellbird to test topography ----


data <- subset_samples(ITSrel.count, Location %in% "Bellbird" )

## 4.1. Ectomycorrhizal fungi ----

ecm.list=row.names(guild_table)[which(guild_table$primary_lifestyle=="ectomycorrhizal")] # Filter for ECM taxa
data_ecm=prune_taxa(ecm.list, data)

data_ecm = prune_taxa(taxa_sums(data_ecm) > 0, data_ecm) # Remove OTUs with zero count
ecm <- prune_samples(sample_sums(data_ecm) > 0, data_ecm) # Remove samples with zero OTUs 

### 4.1.1. Calculate richness ----


meta = data.frame(sample_data(ecm))
filtered_reads <- sample_sums((ecm))
meta <- cbind(meta, filtered_reads)
diversity <-estimate_richness(ecm, measures=c("Observed", "Chao1", "Shannon", "Simpson"))
alphadiv_ecm <- cbind(diversity, meta)

# Read in file with elevation of each scat sample collection at Bellbird
topography <- read.csv('data/topography.csv')
topography <- topography %>% rename(Elevation = Elevation.of.trap.across.trapping.grid..m.a.s.l.)

# Merge dataframe with alpha diversity data 
ecm_topography <- merge(alphadiv_ecm, topography, by = "Sample")
ecm_topography$Elevation <- as.numeric(ecm_topography$Elevation) 
ecm_topography <- ecm_topography %>% filter(Elevation != 0) # Remove records missing an elevation

### 4.1.2. GLM ----


nb_model <- glm.nb(Observed ~  Elevation , data = ecm_topography) 
summary(nb_model)

season_nb_glm  <- glmmTMB(Observed ~ Elevation + (1 | Season), family = nbinom2(link = "log"), data = ecm_topography) # Test the effect of Elevation when season is included as random effect
summary(season_nb_glm)

# Plot
ggplot(ecm_topography, aes(x = Elevation, y = Observed )) +
  geom_point() + geom_smooth(method = "lm")+ labs(title = "ECM richness in scats collected across a topographic gradient", x = "Elevation (m.a.s.l)", y = "ECM Richness")




## 4.2. Truffle-like Ectomycorrhizal fungi ----

ecm.list=row.names(guild_table)[which(guild_table$primary_lifestyle=="ectomycorrhizal")] # Filter for ECM taxa
data_ecm=prune_taxa(ecm.list, data)

hypo.list=row.names(guild_table)[which(guild_table$Fruitbody_type=="gasteroid-hypogeous")]
data_truffle=prune_taxa(hypo.list, data_ecm)

data_truffle = prune_taxa(taxa_sums(data_truffle) > 0, data_truffle) # Remove OTUs with zero count
truffle_ecm <- prune_samples(sample_sums(data_truffle) > 0, data_truffle) # Remove samples with zero OTUs 


### 4.2.1. Calculate richness ----


meta = data.frame(sample_data(truffle_ecm))
filtered_reads <- sample_sums((truffle_ecm))
meta <- cbind(meta, filtered_reads)
diversity <-estimate_richness(truffle_ecm, measures=c("Observed", "Chao1", "Shannon", "Simpson"))
alphadiv_truffle_ecm <- cbind(diversity, meta)

# Read in file with elevation of each scat sample collection at Bellbird
topography <- read.csv('data/topography.csv')
topography <- topography %>% rename(Elevation = Elevation.of.trap.across.trapping.grid..m.a.s.l.)

# Merge dataframe with alpha diversity data 
truffle_ecm_topography <- merge(alphadiv_truffle_ecm, topography, by = "Sample")
truffle_ecm_topography$Elevation <- as.numeric(truffle_ecm_topography$Elevation) 
truffle_ecm_topography <- truffle_ecm_topography %>% filter(Elevation != 0) # Remove records missing an elevation

### 4.2.2. GLM ----

nb_model <- glm.nb(Observed ~  Elevation , data = truffle_ecm_topography) 
summary(nb_model)

season_nb_glm  <- glmmTMB(Observed ~ Elevation + (1 | Season), family = nbinom2(link = "log"), data = truffle_ecm_topography) # Test the effect of Elevation when season is included as random effect
summary(season_nb_glm)

# Plot
ggplot(truffle_ecm_topography, aes(x = Elevation, y = Observed )) +
  geom_point() + geom_smooth(method = "lm")+ labs(title = "truffle_ecm richness in scats collected across a topographic gradient", x = "Elevation (m.a.s.l)", y = "truffle_ecm Richness")

