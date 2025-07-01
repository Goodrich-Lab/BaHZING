## code to prepare `iHMP microbiome and diet` dataset
library(phyloseq)
library(readr)
library(tidyverse)

Sys.setenv("VROOM_CONNECTION_SIZE" = 10000000)


# The following data is from https://ibdmdb.org/ -----------
# Complete metadata is from here: https://ibdmdb.org/
# This publication:  DOI: 10.1038/s41586-019-1237-9

# Read in metadata
ihmp_complete_metadata <- read_csv("https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/metadata/hmp2_metadata_2018-08-20.csv")

# Read in microbiome data, downloaded from here: https://ibdmdb.org/downloads/html/products_MGX_2017-08-12.html
ihmp_species_full <- read_tsv("https://raw.githubusercontent.com/borenstein-lab/microbiome-metabolome-curated-data/refs/heads/main/data/processed_data/iHMP_IBDMDB_2019/species.counts.tsv")

# Modify metadata names
ihmp_complete_metadata1 <- ihmp_complete_metadata |>
  janitor::clean_names() |>
  tidylog::rename_all(~str_remove(., "_etc")) |>
  tidylog::rename(Sample = external_id) |>
  tidylog::select(-project) |>
  tidylog::filter(Sample %in% ihmp_species_full$Sample,
                  data_type == "metagenomics") |>
  tidylog::rename(soft_drinks = soft_drinks_tea_or_coffee_with_sugar_corn_syrup_maple_syrup_cane_sugar,
                  diet_soft_drinks = diet_soft_drinks_tea_or_coffee_with_sugar_stevia_equal_splenda,
                  fruit_juice = fruit_juice_orange_apple_cranberry_prune,
                  alcohol = alcohol_beer_brandy_spirits_hard_liquor_wine_aperitif,
                  yogurt = yogurt_or_other_foods_containing_active_bacterial_cultures_kefir_sauerkraut,
                  dairy = dairy_milk_cream_ice_cream_cheese_cream_cheese,
                  probiotic = probiotic,
                  fruits_no_juice = fruits_no_juice_apples_raisins_bananas_oranges_strawberries_blueberries,
                  vegetables = vegetables_salad_tomatoes_onions_greens_carrots_peppers_green_beans,
                  beans_soy = beans_tofu_soy_soy_burgers_lentils_mexican_beans_lima_beans,
                  whole_grains = whole_grains_wheat_oats_brown_rice_rye_quinoa_wheat_bread_wheat_pasta,
                  starch = starch_white_rice_bread_pizza_potatoes_yams_cereals_pancakes,
                  eggs = eggs,
                  processed_meat = processed_meat_other_red_or_white_meat_such_as_lunch_meat_ham_salami_bologna,
                  red_meat = red_meat_beef_hamburger_pork_lamb,
                  white_meat = white_meat_chicken_turkey,
                  shellfish = shellfish_shrimp_lobster_scallops,
                  fish = fish_fish_nuggets_breaded_fish_fish_cakes_salmon_tuna,
                  sweets = sweets_pies_jam_chocolate_cake_cookies)

# Filter to only include the first visit for each participant
metadata <- ihmp_complete_metadata1 |>
  tidylog::filter(!is.na(water)) |>
  group_by(participant_id) |>
  tidylog::filter(visit_num == min(visit_num)) |>
  ungroup()

# Filter species dataframe
ihmp_species <- ihmp_species_full |>
  tidylog::filter(Sample %in% metadata$Sample)


# 1. Diet data ------

## a) Format diet data ------
# get diet data column names:
diet_var_name <- metadata |>
  tidylog::select(soft_drinks:sweets,
                  -x1_in_the_past_2_weeks_have_you_received_any_of_the_following_medications) |>
  colnames()

# Create numeric formatted diet data
metadata <- metadata |>
  mutate_at(
    .vars = diet_var_name,
    .funs = list("dietnum" = ~case_when(. == "No, I did not consume these products in the last 7 days" ~ 0,
                                        . == "Within the past 4 to 7 days" ~ 1,
                                        . == "Within the past 2 to 3 days" ~ 2,
                                        . == "Yesterday, 1 to 2 times"     ~ 3,
                                        . == "Yesterday, 3 or more times"  ~ 4))) |>
  mutate_at(.vars = metadata |>
              tidylog::select(ends_with("_dietnum")) |> colnames(),
            .funs = list("fct" = ~as.factor(.)))

# Select only diet
diet_dat <- metadata |>
  tidylog::select(ends_with("_dietnum"))

diet_dat_fct <- metadata |>
  tidylog::select(ends_with("_dietnum_fct"))

## b) Get summary table of diet data ####
# diet_dat |>
#   gtsummary::tbl_summary() |>
#   as_kable() |>
#   kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
#                             full_width = F,
#                             position = "center")

# Plot histogram of diet data
# diet_dat |>
#   gather(cols, value) %>%
#   ggplot(aes(x = value)) +
#   geom_histogram() +
#   facet_wrap(.~cols)

# 2. Covariates ---------

## a. Connvert categorical variables to indicator variables
education <- paste0("educ",0:7)
metadata <- metadata %>%
  mutate(sex=ifelse(sex=="Female",0,1)) %>%
  mutate(race0=ifelse(race=="White",1,0),
         race1=ifelse(race=="Black or African American",1,0),
         race2=ifelse(race=="American Indian or Alaska Native",1,0),
         race3=ifelse(race=="More than one race",1,0)) %>%
  mutate(educ0=ifelse(education_level=="7th grade or less",1,0),
         educ1=ifelse(education_level=="Some high school",1,0),
         educ2=ifelse(education_level=="High school graduate or GED",1,0),
         educ3=ifelse(education_level=="Some college, no degree",1,0),
         educ4=ifelse(education_level=="Bachelor's degree",1,0),
         educ5=ifelse(education_level=="Master's degree",1,0),
         educ6=ifelse(education_level=="Professional/Doctoral degree",1,0),
         educ7=ifelse(education_level=="Unknown/Not Reported",1,0)) %>%
  mutate(across(all_of(education),~ifelse(is.na(.),0,.)))

#FIx NAs in diet data
diet <- colnames(diet_dat)
metadata <- metadata %>%
  mutate(across(all_of(diet),~ifelse(is.na(.),0,.)))

#Fix NA in age data
metadata <- metadata %>%
  mutate(consent_age=ifelse(is.na(consent_age),0,consent_age))

# 3. Get taxa, otu, and sample tables ----------

## a. Get otu table --------
otumat <- ihmp_species |>
  column_to_rownames("Sample") |>
  as.matrix() |>
  t()

# Get percent of zeros
otumat_0s <- otumat == 0
pct0 <- rowSums(otumat_0s)/ncol(otumat_0s)

otumat <- otumat[pct0 != 1,]


## a) Get taxa table --------
# Get microbiome column names
species_metadata <- data.frame(full = rownames(otumat))

# Seperate data by ";"
species_metadata <- species_metadata |>
  tidyr::separate_wider_delim(col = full,
                              delim = ";",
                              names = c("Domain", "Phylum", "Class", "Order",
                                        "Family", "Genus", "Species"),
                              cols_remove = FALSE)

taxmat <- species_metadata |>
  column_to_rownames("full") |>
  as.matrix()


## c. Get sample table --------
SAMPDAT <- metadata |>
  column_to_rownames("Sample") |>
  sample_data()

## d. Convert to phyloseq object ------
library("phyloseq")
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)

sample_names(SAMPDAT)
sample_names(OTU)

iHMP = phyloseq(OTU, TAX, SAMPDAT)


# Save dataset
usethis::use_data(iHMP, overwrite = TRUE)


#4. Create reduced Dataset for testing -----
filterTaxaByPrevalence <- function(ps, percentSamplesPresentIn){
  #define threshold
  prevalenceThreshold <- percentSamplesPresentIn * nsamples(ps)
  #apply a function to all samples that determines if a taxa is fully absent from a particular sample
  toKeep <- apply(data.frame(otu_table(ps)), 2, function(taxa) return(sum(taxa > 0) > prevalenceThreshold))
  #this is placed into toKeep which is a TRUE/FALSE vector of whether or not a that taxa is present in enough samples
  ps_filt <- prune_taxa(toKeep, ps)
  return(ps_filt)
}

#Transpose phyloseq object to taxa are columns
iHMP <- t(iHMP)
#Reduce dataset
iHMP_Reduced <- filterTaxaByPrevalence(iHMP, 0.97)
#31961 Taxa before filtering
#2913 taxa after filtering

# Save dataset
usethis::use_data(iHMP_Reduced, overwrite = TRUE)

