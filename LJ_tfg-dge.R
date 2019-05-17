# Analysis of differential expression for temperature and Fe metatranscriptomics.

library(magrittr)
library(ggplot2)
library(beepr)
library(dplyr)
library(ggfortify)
library(testthat)
library(reshape2)
library(tidyverse)
library(edgeR)
library(gridExtra)

# load group normalized data
#LJ can also do this for MCL data.. Import the MCL data, combine the cluster# and taxa to make clusters unique
#then change the name of the taxa.group column to grpnorm_taxgrp to fit in with the rest of this code. 
#grp_norm <- read.csv("allMCL.csv")
#grp_norm$orf_id <- paste (grp_norm$cluster, grp_norm$Taxa.group)
#colnames(grp_norm)[1] <- "grpnorm_taxgrp"

grp_norm <- read.csv("annotation_allTFG.grpnorm_mmetsp_fc_pn_reclassified.edgeR.csv")



# read in experiment key with associated treatments
tfg_key <- read.csv("tfg_key.csv")
tfg_key2 <- tfg_key %>% 
  dplyr::filter(time == 5 | time == 0, 
                b12 != 1 | is.na(b12))

# subset large original dataframe
grp_norm_sub <- grp_norm %>% dplyr::select('orf_id',# 'group', 
                                           'grpnorm_taxgrp', 'TFG_t5_1A', 'TFG_t5_1B', 'TFG_t5_1C', 
                                           'TFG_t5_3A', 'TFG_t5_3B', 'TFG_t5_3C', 'TFG_t5_4A', 
                                           'TFG_t5_4B', 'TFG_t5_4C', 'TFG_t5_6A', 'TFG_t5_6B', 
                                           'TFG_t5_6C', 'TFG_t5_7A', 'TFG_t5_7B', 'TFG_t5_7C', 
                                           'TFG_t5_9A', 'TFG_t5_9B', 'TFG_t5_9C', 'TFG_t0_A', 
                                           'TFG_t0_B', 'TFG_t0_C')

grp_norm_sub_not0 <- grp_norm %>% dplyr::select('orf_id',# 'group', 
                                           'grpnorm_taxgrp', 'TFG_t5_1A', 'TFG_t5_1B', 'TFG_t5_1C', 
                                           'TFG_t5_3A', 'TFG_t5_3B', 'TFG_t5_3C', 'TFG_t5_4A', 
                                           'TFG_t5_4B', 'TFG_t5_4C', 'TFG_t5_6A', 'TFG_t5_6B', 
                                           'TFG_t5_6C', 'TFG_t5_7A', 'TFG_t5_7B', 'TFG_t5_7C', 
                                           'TFG_t5_9A', 'TFG_t5_9B', 'TFG_t5_9C')

# change the row names to be the ORF ID
rownames(grp_norm_sub) <- grp_norm_sub$orf_id

# combine the expression values all into one column, with each referring to a sample ID
grp_norm_t <- melt(grp_norm_sub, variable.name = 'sample_key', 
                   value.name = 'expression_value') %>% 
              inner_join(y = tfg_key, by = 'sample_key')

# read in normalization factors
norm_factors <- read.csv("norm_factors_mmetsp_fc_pn_reclassified.csv")

norm_factors_sub <- norm_factors %>% dplyr::select('group', 'TFG_t5_1A', 'TFG_t5_1B', 'TFG_t5_1C', 
                                                   'TFG_t5_3A', 'TFG_t5_3B', 'TFG_t5_3C', 'TFG_t5_4A', 
                                                   'TFG_t5_4B', 'TFG_t5_4C', 'TFG_t5_6A', 'TFG_t5_6B', 
                                                   'TFG_t5_6C', 'TFG_t5_7A', 'TFG_t5_7B', 'TFG_t5_7C', 
                                                   'TFG_t5_9A', 'TFG_t5_9B', 'TFG_t5_9C')

# making a dataframe of taxon specific normalization factors
tax_specific_abundance_norms <- data.frame(grpnorm_taxgrp = norm_factors$group,
                                           norm_factors_group = rowMeans(norm_factors_sub[,c(2:ncol(norm_factors_sub))]))


# grouping of treatments
samples_present <- inner_join(x = data.frame(sample_key = names(grp_norm_sub_not0)[3:23]), 
                              y = tfg_key, by  = 'sample_key')
treatment_groups <- paste(samples_present$temp, samples_present$fe)

## multiplying abundance scaling factors
grp_norm_sub_not0 <- inner_join(x = grp_norm_sub_not0, 
                                y = tax_specific_abundance_norms, 
                                by = 'grpnorm_taxgrp')
grp_norm_sub_scaled <- grp_norm_sub_not0[c(3:20)]*grp_norm_sub_not0$norm_factors_group


# making a 'DGE List' for edgeR differential expression analysis
dge_list_tfg <- DGEList(counts = grp_norm_sub_scaled,
                        group = treatment_groups,
                        remove.zeros = FALSE, 
                        genes = grp_norm_sub_not0$orf_id)

## used for old analysis with traditional model formulation
# fe <- as.logical(samples_present$fe)
# temperature_factor <- as.factor(temperature)

groups <- paste('fe', samples_present$fe, 'temp', samples_present$temp, sep = '_')
groups <- gsub(pattern = "-", replacement = "", x = groups)

# model matrix with each individual group as a unique factor
tfg_mm_factors <- model.matrix(~0 + groups) 

# tfg_mm_contin is the old model formulation
# tfg_mm_contin <- model.matrix(~fe*temperature_factor) 

# estimating dispersion parameters
tfg_disp_factors <- estimateGLMCommonDisp(dge_list_tfg, 
                         design = tfg_mm_factors)
# tfg_disp_contin <- estimateGLMCommonDisp(dge_list_tfg, 
#                          design = tfg_mm_contin)

# fitting gene-wise models
tfg_fit_factors <- glmQLFit(tfg_disp_factors, design = tfg_mm_factors, robust = TRUE)
# tfg_fit_contin <- glmQLFit(tfg_disp_contin, design = tfg_mm_contin, robust = TRUE)

# quasi likelihood DE test

# setting up contrasts to look at fe, temp, and interactions. See these links for guidance:
# for contrast formulation of temp_3vstemp_05:
# https://support.bioconductor.org/p/80224/
# for guidance on setting up model contrasts, see section 'Complicated Contrasts'
#https://bioconductor.org/packages/release/workflows/vignettes/RnaSeqGeneEdgeRQL/inst/doc/edgeRQL.html
# for discussion about using contrasts over traditional model formulation:
# https://support.bioconductor.org/p/115795/#115831

# the first set of contrasts are pairwise contrasts.
my_contrasts <- makeContrasts(fe0_temp3vstemp05 = groupsfe_0_temp_3 - groupsfe_0_temp_0.5, 
                              fe0_temp6vstemp05 =  groupsfe_0_temp_6 - groupsfe_0_temp_0.5,
                              fe0_temp6vstemp3 = groupsfe_0_temp_6 - groupsfe_0_temp_3,
                              fe2_temp3vstemp05 = groupsfe_2_temp_3 - groupsfe_2_temp_0.5, 
                              fe2_temp6vstemp05 = groupsfe_2_temp_6 - groupsfe_2_temp_0.5, 
                              fe2_temp6vstemp3 = groupsfe_2_temp_6 - groupsfe_2_temp_3,
                              # Fe pairwise comparisons
                              temp05_fe2vsfe0 = groupsfe_2_temp_0.5 - groupsfe_0_temp_0.5,
                              temp3_fe2vsfe0 = groupsfe_2_temp_3 - groupsfe_0_temp_3,
                              temp6_fe2vsfe0 = groupsfe_2_temp_6 - groupsfe_0_temp_6,
                              # fe_all contrast: all the coefficients associated with iron compared with those without iron. testing overall expression of +Fe treatments
                              fe_all = (groupsfe_2_temp_0.5 + groupsfe_2_temp_3 + groupsfe_2_temp_6 - groupsfe_0_temp_0.5 - groupsfe_0_temp_3 - groupsfe_0_temp_6)/3, 
                              # testing overall expression of Temp 3 treatments compared with temp 0.5
                              temp_3vstemp_05 = ((groupsfe_2_temp_3 - groupsfe_2_temp_0.5) + (groupsfe_0_temp_3 - groupsfe_0_temp_0.5))/2,
                              # testing overall expression of Temp 6 treatments compared with temp 0.5
                              temp_6vstemp_05 = ((groupsfe_2_temp_6 - groupsfe_2_temp_0.5) + (groupsfe_0_temp_6 - groupsfe_0_temp_0.5))/2,
                              # testing overall expression fo Temp 6 treatments compared with temp 3
                              temp_6vstemp_3 = ((groupsfe_2_temp_6 - groupsfe_2_temp_3) + (groupsfe_0_temp_6 - groupsfe_0_temp_3))/2,
                              # testing interaction between temp and iron, between temps 6 and 0.5
                              tempfe_interaction_temp_6vstemp_05 = (groupsfe_2_temp_6 - groupsfe_0_temp_6) - (groupsfe_2_temp_0.5 - groupsfe_0_temp_0.5),
                              # testing interaction between temp and iron, between temps 6 and 0.5
                              tempfe_interaction_temp_3vstemp_05 = (groupsfe_2_temp_3 - groupsfe_0_temp_3) - (groupsfe_2_temp_0.5 - groupsfe_0_temp_0.5),
                              levels = tfg_mm_factors)


# pairwise comparisons

# temperature pairwise comparisons
# low iron
fe0_temp3vstemp05_test <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'fe0_temp3vstemp05'])
fe0_temp6vstemp05_test <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'fe0_temp6vstemp05'])
fe0_temp6vstemp3_test <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'fe0_temp6vstemp3'])
# high iron
fe2_temp3vstemp05_test <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'fe2_temp3vstemp05'])
fe2_temp6vstemp05_test <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'fe2_temp6vstemp05'])
fe2_temp6vstemp3_test <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'fe2_temp6vstemp3'])

# fe pairwise comparisons
temp05_fe2vsfe0_test <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'temp05_fe2vsfe0'])
temp3_fe2vsfe0_test <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'temp3_fe2vsfe0'])
temp6_fe2vsfe0_test <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'temp6_fe2vsfe0'])





# overall Fe DE values testing
fe_all_test <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'fe_all'])
# temperature comparison values
temp_3vstemp_05_test <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'temp_3vstemp_05'])
temp_6vstemp_05_test <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'temp_6vstemp_05'])
temp_6vstemp_3_test <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'temp_6vstemp_3'])
#temperature and Fe 'interaction' values
tempfe_interaction_temp_6vstemp_05_test <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'tempfe_interaction_temp_6vstemp_05'])
tempfe_interaction_temp_3vstemp_05_test <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'tempfe_interaction_temp_3vstemp_05'])

# function for processing topTags() output that spits out the logFC and adjusted p value
out_toptag <- function(qlf_test_output){
  # qlf_test_output <- fe_all_test
  top_out <- topTags(qlf_test_output, n = Inf, adjust.method = 'BH', sort.by = 'none')
  # get name of input
  nm <- deparse(substitute(qlf_test_output))
  nm_adj <- gsub(pattern = "_test", replacement = "", x = nm)
  
  # making the foldchange and p value columns
  foldchange <- top_out[[1]]$logFC
  pvalue <- top_out[[1]]$FDR
  
  result_df <- data.frame(foldchange, pvalue)
  
  #adjusting names so that they have the test name appended with 'pvalue' and 'foldchange'
  names(result_df) <- paste0(nm_adj, "_", names(result_df))
  
  return(result_df)
}


## up / down / none decisions about DE genes
# pairwise comparisons
# temperature comparisons:
# low iron
fe0_temp3vstemp05_decide <- decideTestsDGE(fe0_temp3vstemp05_test)
fe0_temp3vstemp05_info <- out_toptag(fe0_temp3vstemp05_test)

fe0_temp6vstemp05_decide <- decideTestsDGE(fe0_temp6vstemp05_test)
fe0_temp6vstemp05_info <- out_toptag(fe0_temp6vstemp05_test)

fe0_temp6vstemp3_decide <- decideTestsDGE(fe0_temp6vstemp3_test)
fe0_temp6vstemp3_info <- out_toptag(fe0_temp6vstemp3_test)

# high iron
fe2_temp3vstemp05_decide <- decideTestsDGE(fe2_temp3vstemp05_test)
fe2_temp3vstemp05_info <- out_toptag(fe2_temp3vstemp05_test)

fe2_temp6vstemp05_decide <- decideTestsDGE(fe2_temp6vstemp05_test)
fe2_temp6vstemp05_info <- out_toptag(fe2_temp6vstemp05_test)

fe2_temp6vstemp3_decide <- decideTestsDGE(fe2_temp6vstemp3_test)
fe2_temp6vstemp3_info <- out_toptag(fe2_temp6vstemp3_test)

# fe pairwise comparisons
temp05_fe2vsfe0_decide <- decideTestsDGE(temp05_fe2vsfe0_test)
temp05_fe2vsfe0_info <- out_toptag(temp05_fe2vsfe0_test)

temp3_fe2vsfe0_decide <- decideTestsDGE(temp3_fe2vsfe0_test)
temp3_fe2vsfe0_info <- out_toptag(temp3_fe2vsfe0_test)

temp6_fe2vsfe0_decide <- decideTestsDGE(temp6_fe2vsfe0_test)
temp6_fe2vsfe0_info <- out_toptag(temp6_fe2vsfe0_test)


# Fe decision
fe_all_decide <- decideTestsDGE(fe_all_test)
fe_all_info <- out_toptag(fe_all_test)

# Temperature decision
temp_3vstemp_05_decide <- decideTestsDGE(temp_3vstemp_05_test)
temp_3vstemp_05_info <- out_toptag(temp_3vstemp_05_test)

temp_6vstemp_05_decide <- decideTestsDGE(temp_6vstemp_05_test)
temp_6vstemp_05_info <- out_toptag(temp_6vstemp_05_test)

temp_6vstemp_3_decide <- decideTestsDGE(temp_6vstemp_3_test)
temp_6vstemp_3_info <- out_toptag(temp_6vstemp_3_test)

# interaction DE decision
tempfe_interaction_temp_6vstemp_05_decide <- decideTestsDGE(tempfe_interaction_temp_6vstemp_05_test)
tempfe_interaction_temp_6vstemp_05_info <- out_toptag(tempfe_interaction_temp_6vstemp_05_test)

tempfe_interaction_temp_3vstemp_05_decide <- decideTestsDGE(tempfe_interaction_temp_3vstemp_05_test)
tempfe_interaction_temp_3vstemp_05_info <- out_toptag(tempfe_interaction_temp_3vstemp_05_test)

# checking total number of DE genes for each test
# fe_all_decide@.Data %>% abs() %>% sum()
# 
# temp_3vstemp_05_decide@.Data %>% abs() %>% sum()
# temp_6vstemp_05_decide@.Data %>% abs() %>% sum()
# temp_6vstemp_3_decide@.Data %>% abs() %>% sum()
# 
# tempfe_interaction_temp_3vstemp_05_decide@.Data %>% abs() %>% sum()
# tempfe_interaction_temp_6vstemp_05_decide@.Data %>% abs() %>% sum()

# getting only the rows which do not have na for grpnorm_taxgrp
#grp_norm_no_grp_na <- grp_norm[!is.na(grp_norm$grpnorm_taxgrp),]
#LJ changed the line above into this because I didn't have any NA values
grp_norm_no_grp_na <- dplyr::filter (grp_norm, !grepl('^$', grpnorm_taxgrp))

# dataframe with actual values used for DE analysis (transformed vals by tax-specific scaling factor)
grp_norm_sub_scaled_copy <- grp_norm_sub_scaled
names(grp_norm_sub_scaled_copy) <- paste0(names(grp_norm_sub_scaled_copy), '-scaled')


## making compiled dataframe
tfg_de_df <- data.frame(Fe_vs_noFe_de = fe_all_decide %>% as.vector(),
                        Fe_vs_noFe_pvalue = fe_all_info$fe_all_pvalue,
                        Fe_vs_noFe_foldchange = fe_all_info$fe_all_foldchange,
                        # temperature comparisons
                        temp_3C_vs_0C_de = temp_3vstemp_05_decide %>% as.vector(),
                        temp_3C_vs_0C_pvalue = temp_3vstemp_05_info$temp_3vstemp_05_pvalue,
                        temp_3C_vs_0C_foldchange = temp_3vstemp_05_info$temp_3vstemp_05_foldchange,
                        
                        temp_6C_vs_0C_de = temp_6vstemp_05_decide %>% as.vector(),
                        temp_6C_vs_0C_pvalue = temp_6vstemp_05_info$temp_6vstemp_05_pvalue,
                        temp_6C_vs_0C_foldchange = temp_6vstemp_05_info$temp_6vstemp_05_foldchange,
                        
                        temp_6C_vs_3C_de = temp_6vstemp_3_decide %>% as.vector(),
                        temp_6C_vs_3C_pvalue = temp_6vstemp_3_info$temp_6vstemp_3_pvalue,
                        temp_6C_vs_3C_foldchange = temp_6vstemp_3_info$temp_6vstemp_3_foldchange,
                        
                        int_fe_temp_6C_vs_0C_de = tempfe_interaction_temp_6vstemp_05_decide %>% as.vector(),
                        int_fe_temp_6C_vs_0C_pvalue = tempfe_interaction_temp_6vstemp_05_info$tempfe_interaction_temp_6vstemp_05_pvalue,
                        int_fe_temp_6C_vs_0C_foldchange = tempfe_interaction_temp_6vstemp_05_info$tempfe_interaction_temp_6vstemp_05_foldchange,
                        int_fe_temp_3C_vs_0C_de = tempfe_interaction_temp_3vstemp_05_decide %>% as.vector(),
                        int_fe_temp_3C_vs_0C_pvalue = tempfe_interaction_temp_3vstemp_05_info$tempfe_interaction_temp_3vstemp_05_pvalue,
                        int_fe_temp_3C_vs_0C_foldchange = tempfe_interaction_temp_3vstemp_05_info$tempfe_interaction_temp_3vstemp_05_foldchange,
                        # pairwise comparisons
                        noFe_3C_vs_noFe_0C_de = fe0_temp3vstemp05_decide %>% as.vector(),
                        noFe_3C_vs_noFe_0C_pvalue = fe0_temp3vstemp05_info$fe0_temp3vstemp05_pvalue,
                        noFe_3C_vs_noFe_0C_foldchange = fe0_temp3vstemp05_info$fe0_temp3vstemp05_foldchange,
                        
                        Fe_3C_vs_Fe_0C_de = fe2_temp3vstemp05_decide %>% as.vector(),
                        Fe_3C_vs_Fe_0C_pvalue = fe2_temp3vstemp05_info$fe2_temp3vstemp05_pvalue,
                        Fe_3C_vs_Fe_0C_foldchange = fe2_temp3vstemp05_info$fe2_temp3vstemp05_foldchange,
                        
                        noFe_6C_vs_noFe_0C_de = fe0_temp6vstemp05_decide %>% as.vector(),
                        noFe_6C_vs_noFe_0C_pvalue = fe0_temp6vstemp05_info$fe0_temp6vstemp05_pvalue,
                        noFe_6C_vs_noFe_0C_foldchange = fe0_temp6vstemp05_info$fe0_temp6vstemp05_foldchange,
                        
                        Fe_6C_vs_Fe_0C_de = fe2_temp6vstemp05_decide %>% as.vector(),
                        Fe_6C_vs_Fe_0C_pvalue = fe2_temp6vstemp05_info$fe2_temp6vstemp05_pvalue,
                        Fe_6C_vs_Fe_0C_foldchange = fe2_temp6vstemp05_info$fe2_temp6vstemp05_foldchange,
                        
                        noFe_6C_vs_noFe_3C_de = fe0_temp6vstemp3_decide %>% as.vector(),
                        noFe_6C_vs_noFe_3C_pvalue = fe0_temp6vstemp3_info$fe0_temp6vstemp3_pvalue,
                        noFe_6C_vs_noFe_3C_foldchange = fe0_temp6vstemp3_info$fe0_temp6vstemp3_foldchange,
                        
                        Fe_6C_vs_Fe_3C_de = fe2_temp6vstemp3_decide %>% as.vector(),
                        Fe_6C_vs_Fe_3C_pvalue = fe2_temp6vstemp3_info$fe2_temp6vstemp3_pvalue,
                        Fe_6C_vs_Fe_3C_foldchange = fe2_temp6vstemp3_info$fe2_temp6vstemp3_foldchange,
                        
                        Fe_0C_vs_noFe_0C_de = temp05_fe2vsfe0_decide %>% as.vector(),
                        Fe_0C_vs_noFe_0C_pvalue = temp05_fe2vsfe0_info$temp05_fe2vsfe0_pvalue,
                        Fe_0C_vs_noFe_0C_foldchange = temp05_fe2vsfe0_info$temp05_fe2vsfe0_foldchange,
                        
                        Fe_3C_vs_noFe_3C_de = temp3_fe2vsfe0_decide %>% as.vector(),
                        Fe_3C_vs_noFe_3C_pvalue = temp3_fe2vsfe0_info$temp3_fe2vsfe0_pvalue,
                        Fe_3C_vs_noFe_3C_foldchange = temp3_fe2vsfe0_info$temp3_fe2vsfe0_foldchange,
                        
                        Fe_6C_vs_noFe_6C_de = temp6_fe2vsfe0_decide %>% as.vector(),
                        Fe_6C_vs_noFe_6C_pvalue = temp6_fe2vsfe0_info$temp6_fe2vsfe0_pvalue,
                        Fe_6C_vs_noFe_6C_foldchange = temp6_fe2vsfe0_info$temp6_fe2vsfe0_foldchange,
                        
                        # including dataframe with annotations
                        grp_norm_no_grp_na,
                        # including scaled values,
                        grp_norm_sub_scaled_copy,
                        # columns for plotting DE analysis
                        fe_neg_de = ifelse(fe_all_decide == -1, yes = -1, no = 0) %>% as.vector(),
                        fe_pos_de = ifelse(fe_all_decide == 1, yes = 1, no = 0) %>% as.vector(),
                        # temperature is a bit harder than Fe. The gene has to be DE in either comparisons from 3vs0.5 *or* 6vs0.5 to be counted in this column. But if the gene is negatively expressed from 3vs0.5, and positively expressed from 6vs0.5, it would be excluded.
                        temp_neg_de = ifelse((temp_3vstemp_05_decide == -1 | temp_6vstemp_05_decide == -1) &
                                               (temp_3vstemp_05_decide != 1 | temp_6vstemp_05_decide != 1), 
                                             yes = -1, no = 0) %>% as.vector(),
                        
                        temp_pos_de = ifelse((temp_3vstemp_05_decide == 1 | temp_6vstemp_05_decide == 1) &
                                               (temp_3vstemp_05_decide != -1 | temp_6vstemp_05_decide != -1), 
                                             yes = 1, no = 0) %>% as.vector(),
                        tempfe_neg_de = ifelse((tempfe_interaction_temp_6vstemp_05_decide == -1 | tempfe_interaction_temp_3vstemp_05_decide == -1) & 
                                                 (tempfe_interaction_temp_6vstemp_05_decide != 1 & tempfe_interaction_temp_3vstemp_05_decide != 1), 
                                               yes = -1, no = 0) %>% as.vector(),
                        tempfe_pos_de = ifelse((tempfe_interaction_temp_6vstemp_05_decide == 1 | tempfe_interaction_temp_3vstemp_05_decide == 1) & 
                                                 (tempfe_interaction_temp_6vstemp_05_decide != -1 & tempfe_interaction_temp_3vstemp_05_decide != -1), 
                                               yes = 1, no = 0) %>% as.vector(),
                        # abs_fe is for ranking taxonomic groups for plotting later on.
                        abs_fe = ifelse(abs(fe_all_decide) == 1, yes = 1, no = 0) %>% as.vector())


# old analysis

# fe0_temp05vs3 <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'fe0_temp05vstemp3'])
# fe0_temp05vs6 <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'fe0_temp05vstemp6'])
# fe0_temp3vs6 <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'fe0_temp3vstemp6'])
# 
# fe2_temp05vs3 <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'fe2_temp05vstemp3'])
# fe2_temp05vs6 <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'fe2_temp05vstemp6'])
# fe2_temp3vs6 <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'fe2_temp3vstemp6'])
# 
# 
# fe0_temp05vs3_decide <- decideTestsDGE(object = fe0_temp05vs3)
# fe0_temp05vs6_decide <- decideTestsDGE(object = fe0_temp05vs6)
# fe0_temp3vs6_decide <- decideTestsDGE(object = fe0_temp3vs6)
# 
# fe2_temp05vs3_decide <- decideTestsDGE(object = fe2_temp05vs3)
# fe2_temp05vs6_decide <- decideTestsDGE(object = fe2_temp05vs6)
# fe2_temp3vs6_decide <- decideTestsDGE(object = fe2_temp3vs6)
# 
# 
# fe2_temp05vs3_decide@.Data %>% abs() %>% sum()
# fe2_temp05vs6_decide@.Data %>% abs() %>% sum()
# fe0_temp05vs3_decide@.Data %>% abs() %>% sum()
# fe0_temp05vs6_decide@.Data %>% abs() %>% sum()


## getting significance tests of old model formulation

tfg_inter <- glmQLFTest(tfg_fit_contin, coef = 1)
tfg_fe <- glmQLFTest(tfg_fit_contin, coef = 2)
tfg_temp3 <- glmQLFTest(tfg_fit_contin, coef = 3)
tfg_temp6 <- glmQLFTest(tfg_fit_contin, coef = 4)
tfg_fetemp3 <- glmQLFTest(tfg_fit_contin, coef = 5)
tfg_fetemp6 <- glmQLFTest(tfg_fit_contin, coef = 6)
# 
# tfg_inter_decide <- decideTestsDGE(tfg_inter)
# tfg_fe_decide <- decideTestsDGE(tfg_fe)
# tfg_temp3_decide <- decideTestsDGE(tfg_temp3)
# tfg_temp6_decide <- decideTestsDGE(tfg_temp6)
# tfg_fetemp3_decide <- decideTestsDGE(tfg_fetemp3)
# tfg_fetemp6_decide <- decideTestsDGE(tfg_fetemp6)


### old version writing final dataframe

# tfg_de_df <- data.frame(tfg_inter_decide, tfg_fe_decide, tfg_temp3_decide,
#            tfg_temp6_decide, tfg_fetemp3_decide, tfg_fetemp6_decide,
#            grp_norm_sub_not0,
#           abs_fe = ifelse(abs(tfg_fe_decide) == 1, yes = 1, no = 0) %>% as.vector(),
#           abs_temp = ifelse(abs(tfg_temp3_decide) == 1 | abs(tfg_temp6_decide) == 1, yes = 1, no = 0) %>% as.vector(),
#           abs_tempfe = ifelse(abs(tfg_fetemp3_decide) == 1 | abs(tfg_fetemp6_decide) == 1, yes = 1, no = 0) %>% as.vector(),
#           # summaries of negative DE by treatment *group*
#           fe_neg_de = ifelse((tfg_fe_decide == -1 | tfg_fe_decide == -1) & 
#                                  (tfg_fe_decide != 1 | tfg_fe_decide!= 1), yes = -1, no = 0) %>% as.vector(),
#           temp_neg_de = ifelse((tfg_temp3_decide == -1 | tfg_temp6_decide == -1) & 
#                                  (tfg_temp3_decide != 1 | tfg_temp6_decide!= 1), yes = -1, no = 0) %>% as.vector(),
#           fetemp_neg_de = ifelse((tfg_fetemp3_decide == -1 | tfg_fetemp6_decide == -1) & 
#                                  (tfg_fetemp3_decide != 1 | tfg_fetemp6_decide!= 1), yes = -1, no = 0) %>% as.vector(),
#           
#           # summaries of positive DE by treatment *group*
#           fe_pos_de = ifelse((tfg_fe_decide == 1 | tfg_fe_decide == 1) & 
#                                  (tfg_fe_decide != -1 | tfg_fe_decide!= -1), yes = 1, no = 0) %>% as.vector(),
#           temp_pos_de = ifelse((tfg_temp3_decide == 1 | tfg_temp6_decide == 1) & 
#                                  (tfg_temp3_decide != -1 | tfg_temp6_decide!= -1), yes = 1, no = 0) %>% as.vector(),
#           fetemp_pos_de = ifelse((tfg_fetemp3_decide == 1 | tfg_fetemp6_decide == 1) & 
#                                  (tfg_fetemp3_decide != -1 | tfg_fetemp6_decide!= -1), yes = 1, no = 0) %>% as.vector())

write.csv(tfg_de_df, file = 'LJ_tfg_de_annotation_allTFG.grpnorm_mmetsp_fc_pn_reclassified.csv')

### plotting de results by taxon

# blank plot just used for the axis.
blank_p <- tfg_de_df %>% 
  group_by(grpnorm_taxgrp) %>% 
  summarise(fe_grp = sum(abs_fe)) %>%  
  ggplot(aes(x = fe_grp, 
             y = fct_reorder(grpnorm_taxgrp, fe_grp))) + 
  # geom_point() +
  theme_bw() +
  xlab('') + 
  theme(axis.text.x = element_text(colour = 'white'), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), 
        axis.ticks = element_blank()) +
  ylab('');blank_p

# old plots that show total DE, not subdivided into positive and negative DE:

# p1_fe <- tfg_de_df %>% 
#   group_by(grpnorm_taxgrp) %>% 
#   summarise(fe_grp = sum(abs_fe)) %>%  
#   ggplot(aes(x = fe_grp, y = fct_reorder(grpnorm_taxgrp, fe_grp))) + 
#   geom_point(size = 3, alpha = 0.5) +
#   theme_bw() +
#   xlab('Iron-Related DE') + 
#   ylab('') +
#   theme(axis.text.y = element_blank());p1_fe
# 
# p2_temp <- tfg_de_df %>% 
#   group_by(grpnorm_taxgrp) %>% 
#   summarise(fe_grp = sum(abs_fe),
#             temp_grp = sum(abs_temp)) %>% 
#   ggplot(aes(x = temp_grp, y = fct_reorder(grpnorm_taxgrp, fe_grp))) + 
#   geom_point(size = 3, alpha = 0.5) +
#   theme_bw() + 
#   ylab('') + 
#   xlab('Temperature-Related DE') +
#   theme(axis.text.y = element_blank());p2_temp
# 
# p3_tempfe <- tfg_de_df %>% 
#   group_by(grpnorm_taxgrp) %>% 
#   summarise(fe_grp = sum(abs_fe),
#             fetemp_grp = sum(abs_tempfe)) %>% 
#   ggplot(aes(x = fetemp_grp, y = fct_reorder(grpnorm_taxgrp, fe_grp))) + 
#   geom_point(size = 3, alpha = 0.5) +
#   theme_bw() + 
#   ylab('') + 
#   xlab('Temperature-Iron-Related DE') +
#   theme(axis.text.y = element_blank());p3_tempfe
# 
# 
# grid.arrange(blank_p, p1_fe, p2_temp, p3_tempfe, nrow = 1)


# DE plot with positive and negative differences shown.

pn_fe <- tfg_de_df %>% 
  group_by(grpnorm_taxgrp) %>% 
  summarise(fe_grp = sum(abs_fe),
            fe_grp_neg = sum(fe_neg_de),
            fe_grp_pos = sum(fe_pos_de)) %>%  
  ggplot(aes(x = fe_grp_neg, y = fct_reorder(grpnorm_taxgrp, fe_grp))) + 
  geom_point(size = 3, alpha = 0.5, colour = 'blue') +
  geom_point(aes(x = fe_grp_pos, y = fct_reorder(grpnorm_taxgrp, fe_grp)),
             size = 3, alpha = 0.5, colour = 'red') +
  geom_segment(aes(x = fe_grp_neg, xend = 0, 
                   y = fct_reorder(grpnorm_taxgrp, fe_grp), 
                   yend = fct_reorder(grpnorm_taxgrp, fe_grp)), colour = 'blue', alpha = 0.2, lwd = 2) +
  geom_segment(aes(x = fe_grp_pos, xend = 0, 
                   y = fct_reorder(grpnorm_taxgrp, fe_grp), 
                   yend = fct_reorder(grpnorm_taxgrp, fe_grp)), colour = 'red', alpha = 0.2, lwd = 2) +
  theme_bw() +
  xlab('Iron-Related DE') + 
  ylab('') +
  theme(axis.text.y = element_blank());pn_fe


pn_temp <- tfg_de_df %>% 
  group_by(grpnorm_taxgrp) %>% 
  summarise(fe_grp = sum(abs_fe),
            temp_grp_neg = sum(temp_neg_de),
            temp_grp_pos = sum(temp_pos_de)) %>%  
  ggplot(aes(x = temp_grp_neg, y = fct_reorder(grpnorm_taxgrp, fe_grp))) + 
  geom_point(size = 3, alpha = 0.5, colour = 'blue') +
  geom_point(aes(x = temp_grp_pos, y = fct_reorder(grpnorm_taxgrp, fe_grp)),
             size = 3, alpha = 0.5, colour = 'red') +
  geom_segment(aes(x = temp_grp_neg, xend = 0, 
                   y = fct_reorder(grpnorm_taxgrp, fe_grp), 
                   yend = fct_reorder(grpnorm_taxgrp, fe_grp)), colour = 'blue', alpha = 0.2, lwd = 2) +
  geom_segment(aes(x = temp_grp_pos, xend = 0, 
                   y = fct_reorder(grpnorm_taxgrp, fe_grp), 
                   yend = fct_reorder(grpnorm_taxgrp, fe_grp)), colour = 'red', alpha = 0.2, lwd = 2) +
  theme_bw() +
  xlab('Temperature-Related DE') + 
  ylab('') +
  theme(axis.text.y = element_blank());pn_temp

pn_tempfe <- tfg_de_df %>% 
  group_by(grpnorm_taxgrp) %>% 
  summarise(fe_grp = sum(abs_fe),
            fetemp_grp_neg = sum(tempfe_neg_de),
            fetemp_grp_pos = sum(tempfe_pos_de)) %>%  
  ggplot(aes(x = fetemp_grp_neg, y = fct_reorder(grpnorm_taxgrp, fe_grp))) + 
  geom_point(size = 3, alpha = 0.5, colour = 'blue') +
  geom_point(aes(x = fetemp_grp_pos, y = fct_reorder(grpnorm_taxgrp, fe_grp)),
             size = 3, alpha = 0.5, colour = 'red') +
  geom_segment(aes(x = fetemp_grp_neg, xend = 0, 
                   y = fct_reorder(grpnorm_taxgrp, fe_grp), 
                   yend = fct_reorder(grpnorm_taxgrp, fe_grp)), colour = 'blue', alpha = 0.2, lwd = 2) +
  geom_segment(aes(x = fetemp_grp_pos, xend = 0, 
                   y = fct_reorder(grpnorm_taxgrp, fe_grp), 
                   yend = fct_reorder(grpnorm_taxgrp, fe_grp)), colour = 'red', alpha = 0.2, lwd = 2) +
  theme_bw() +
  xlab('Iron-Temperature-Related DE') + 
  ylab('') +
  theme(axis.text.y = element_blank());pn_tempfe

grid.arrange(blank_p, pn_fe, pn_temp, pn_tempfe, nrow = 1)




##`````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````##
##-------------------------------------------------------------------------------------------------------------------------------------##
#this does the same thing as above, but for MCL data


# Analysis of differential expression for temperature and Fe metatranscriptomics.

library(magrittr)
library(ggplot2)
library(beepr)
library(dplyr)
library(ggfortify)
library(testthat)
library(reshape2)
library(tidyverse)
library(edgeR)
library(gridExtra)

# load group normalized data
#LJ can also do this for MCL data.. Import the MCL data, combine the cluster# and taxa to make clusters unique
#then change the name of the taxa.group column to grpnorm_taxgrp to fit in with the rest of this code. 
grp_norm <- read.csv("allMCL.csv")
grp_norm$orf_id <- paste (grp_norm$cluster, grp_norm$Taxa.group)
colnames(grp_norm)[1] <- "grpnorm_taxgrp"

#grp_norm <- read.csv("annotation_allTFG.grpnorm_mmetsp_fc_pn_reclassified.edgeR.csv")



# read in experiment key with associated treatments
tfg_key <- read.csv("tfg_key.csv")
tfg_key2 <- tfg_key %>% 
  dplyr::filter(time == 5 | time == 0, 
                b12 != 1 | is.na(b12))

# subset large original dataframe
grp_norm_sub <- grp_norm %>% dplyr::select('orf_id',# 'group', 
                                           'grpnorm_taxgrp', 'TFG_t5_1A', 'TFG_t5_1B', 'TFG_t5_1C', 
                                           'TFG_t5_3A', 'TFG_t5_3B', 'TFG_t5_3C', 'TFG_t5_4A', 
                                           'TFG_t5_4B', 'TFG_t5_4C', 'TFG_t5_6A', 'TFG_t5_6B', 
                                           'TFG_t5_6C', 'TFG_t5_7A', 'TFG_t5_7B', 'TFG_t5_7C', 
                                           'TFG_t5_9A', 'TFG_t5_9B', 'TFG_t5_9C', 'TFG_t0_A', 
                                           'TFG_t0_B', 'TFG_t0_C')

grp_norm_sub_not0 <- grp_norm %>% dplyr::select('orf_id',# 'group', 
                                                'grpnorm_taxgrp', 'TFG_t5_1A', 'TFG_t5_1B', 'TFG_t5_1C', 
                                                'TFG_t5_3A', 'TFG_t5_3B', 'TFG_t5_3C', 'TFG_t5_4A', 
                                                'TFG_t5_4B', 'TFG_t5_4C', 'TFG_t5_6A', 'TFG_t5_6B', 
                                                'TFG_t5_6C', 'TFG_t5_7A', 'TFG_t5_7B', 'TFG_t5_7C', 
                                                'TFG_t5_9A', 'TFG_t5_9B', 'TFG_t5_9C')

# change the row names to be the ORF ID
rownames(grp_norm_sub) <- grp_norm_sub$orf_id

# combine the expression values all into one column, with each referring to a sample ID
grp_norm_t <- melt(grp_norm_sub, variable.name = 'sample_key', 
                   value.name = 'expression_value') %>% 
  inner_join(y = tfg_key, by = 'sample_key')

# read in normalization factors
norm_factors <- read.csv("norm_factors_mmetsp_fc_pn_reclassified.csv")

norm_factors_sub <- norm_factors %>% dplyr::select('group', 'TFG_t5_1A', 'TFG_t5_1B', 'TFG_t5_1C', 
                                                   'TFG_t5_3A', 'TFG_t5_3B', 'TFG_t5_3C', 'TFG_t5_4A', 
                                                   'TFG_t5_4B', 'TFG_t5_4C', 'TFG_t5_6A', 'TFG_t5_6B', 
                                                   'TFG_t5_6C', 'TFG_t5_7A', 'TFG_t5_7B', 'TFG_t5_7C', 
                                                   'TFG_t5_9A', 'TFG_t5_9B', 'TFG_t5_9C')

# making a dataframe of taxon specific normalization factors
tax_specific_abundance_norms <- data.frame(grpnorm_taxgrp = norm_factors$group,
                                           norm_factors_group = rowMeans(norm_factors_sub[,c(2:ncol(norm_factors_sub))]))


# grouping of treatments
samples_present <- inner_join(x = data.frame(sample_key = names(grp_norm_sub_not0)[3:23]), 
                              y = tfg_key, by  = 'sample_key')
treatment_groups <- paste(samples_present$temp, samples_present$fe)

## multiplying abundance scaling factors
grp_norm_sub_not0 <- inner_join(x = grp_norm_sub_not0, 
                                y = tax_specific_abundance_norms, 
                                by = 'grpnorm_taxgrp')
grp_norm_sub_scaled <- grp_norm_sub_not0[c(3:20)]*grp_norm_sub_not0$norm_factors_group


# making a 'DGE List' for edgeR differential expression analysis
dge_list_tfg <- DGEList(counts = grp_norm_sub_scaled,
                        group = treatment_groups,
                        remove.zeros = FALSE, 
                        genes = grp_norm_sub_not0$orf_id)

## used for old analysis with traditional model formulation
# fe <- as.logical(samples_present$fe)
# temperature_factor <- as.factor(temperature)

groups <- paste('fe', samples_present$fe, 'temp', samples_present$temp, sep = '_')
groups <- gsub(pattern = "-", replacement = "", x = groups)

# model matrix with each individual group as a unique factor
tfg_mm_factors <- model.matrix(~0 + groups) 

# tfg_mm_contin is the old model formulation
# tfg_mm_contin <- model.matrix(~fe*temperature_factor) 

# estimating dispersion parameters
tfg_disp_factors <- estimateGLMCommonDisp(dge_list_tfg, 
                                          design = tfg_mm_factors)
# tfg_disp_contin <- estimateGLMCommonDisp(dge_list_tfg, 
#                          design = tfg_mm_contin)

# fitting gene-wise models
tfg_fit_factors <- glmQLFit(tfg_disp_factors, design = tfg_mm_factors, robust = TRUE)
# tfg_fit_contin <- glmQLFit(tfg_disp_contin, design = tfg_mm_contin, robust = TRUE)

# quasi likelihood DE test

# setting up contrasts to look at fe, temp, and interactions. See these links for guidance:
# for contrast formulation of temp_3vstemp_05:

# the first set of contrasts are pairwise contrasts.
my_contrasts <- makeContrasts(fe0_temp3vstemp05 = groupsfe_0_temp_3 - groupsfe_0_temp_0.5, 
                              fe0_temp6vstemp05 =  groupsfe_0_temp_6 - groupsfe_0_temp_0.5,
                              fe0_temp6vstemp3 = groupsfe_0_temp_6 - groupsfe_0_temp_3,
                              fe2_temp3vstemp05 = groupsfe_2_temp_3 - groupsfe_2_temp_0.5, 
                              fe2_temp6vstemp05 = groupsfe_2_temp_6 - groupsfe_2_temp_0.5, 
                              fe2_temp6vstemp3 = groupsfe_2_temp_6 - groupsfe_2_temp_3,
                              # Fe pairwise comparisons
                              temp05_fe2vsfe0 = groupsfe_2_temp_0.5 - groupsfe_0_temp_0.5,
                              temp3_fe2vsfe0 = groupsfe_2_temp_3 - groupsfe_0_temp_3,
                              temp6_fe2vsfe0 = groupsfe_2_temp_6 - groupsfe_0_temp_6,
                              # fe_all contrast: all the coefficients associated with iron compared with those without iron. testing overall expression of +Fe treatments
                              fe_all = (groupsfe_2_temp_0.5 + groupsfe_2_temp_3 + groupsfe_2_temp_6 - groupsfe_0_temp_0.5 - groupsfe_0_temp_3 - groupsfe_0_temp_6)/3, 
                              # testing overall expression of Temp 3 treatments compared with temp 0.5
                              temp_3vstemp_05 = ((groupsfe_2_temp_3 - groupsfe_2_temp_0.5) + (groupsfe_0_temp_3 - groupsfe_0_temp_0.5))/2,
                              # testing overall expression of Temp 6 treatments compared with temp 0.5
                              temp_6vstemp_05 = ((groupsfe_2_temp_6 - groupsfe_2_temp_0.5) + (groupsfe_0_temp_6 - groupsfe_0_temp_0.5))/2,
                              # testing overall expression fo Temp 6 treatments compared with temp 3
                              temp_6vstemp_3 = ((groupsfe_2_temp_6 - groupsfe_2_temp_3) + (groupsfe_0_temp_6 - groupsfe_0_temp_3))/2,
                              # testing interaction between temp and iron, between temps 6 and 0.5
                              tempfe_interaction_temp_6vstemp_05 = (groupsfe_2_temp_6 - groupsfe_0_temp_6) - (groupsfe_2_temp_0.5 - groupsfe_0_temp_0.5),
                              # testing interaction between temp and iron, between temps 6 and 0.5
                              tempfe_interaction_temp_3vstemp_05 = (groupsfe_2_temp_3 - groupsfe_0_temp_3) - (groupsfe_2_temp_0.5 - groupsfe_0_temp_0.5),
                              levels = tfg_mm_factors)


# pairwise comparisons

# temperature pairwise comparisons
# low iron
fe0_temp3vstemp05_test <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'fe0_temp3vstemp05'])
fe0_temp6vstemp05_test <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'fe0_temp6vstemp05'])
fe0_temp6vstemp3_test <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'fe0_temp6vstemp3'])
# high iron
fe2_temp3vstemp05_test <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'fe2_temp3vstemp05'])
fe2_temp6vstemp05_test <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'fe2_temp6vstemp05'])
fe2_temp6vstemp3_test <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'fe2_temp6vstemp3'])

# fe pairwise comparisons
temp05_fe2vsfe0_test <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'temp05_fe2vsfe0'])
temp3_fe2vsfe0_test <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'temp3_fe2vsfe0'])
temp6_fe2vsfe0_test <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'temp6_fe2vsfe0'])





# overall Fe DE values testing
fe_all_test <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'fe_all'])
# temperature comparison values
temp_3vstemp_05_test <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'temp_3vstemp_05'])
temp_6vstemp_05_test <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'temp_6vstemp_05'])
temp_6vstemp_3_test <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'temp_6vstemp_3'])
#temperature and Fe 'interaction' values
tempfe_interaction_temp_6vstemp_05_test <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'tempfe_interaction_temp_6vstemp_05'])
tempfe_interaction_temp_3vstemp_05_test <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'tempfe_interaction_temp_3vstemp_05'])

# function for processing topTags() output that spits out the logFC and adjusted p value
out_toptag <- function(qlf_test_output){
  # qlf_test_output <- fe_all_test
  top_out <- topTags(qlf_test_output, n = Inf, adjust.method = 'BH', sort.by = 'none')
  # get name of input
  nm <- deparse(substitute(qlf_test_output))
  nm_adj <- gsub(pattern = "_test", replacement = "", x = nm)
  
  # making the foldchange and p value columns
  foldchange <- top_out[[1]]$logFC
  pvalue <- top_out[[1]]$FDR
  
  result_df <- data.frame(foldchange, pvalue)
  
  #adjusting names so that they have the test name appended with 'pvalue' and 'foldchange'
  names(result_df) <- paste0(nm_adj, "_", names(result_df))
  
  return(result_df)
}


## up / down / none decisions about DE genes
# pairwise comparisons
# temperature comparisons:
# low iron
fe0_temp3vstemp05_decide <- decideTestsDGE(fe0_temp3vstemp05_test)
fe0_temp3vstemp05_info <- out_toptag(fe0_temp3vstemp05_test)

fe0_temp6vstemp05_decide <- decideTestsDGE(fe0_temp6vstemp05_test)
fe0_temp6vstemp05_info <- out_toptag(fe0_temp6vstemp05_test)

fe0_temp6vstemp3_decide <- decideTestsDGE(fe0_temp6vstemp3_test)
fe0_temp6vstemp3_info <- out_toptag(fe0_temp6vstemp3_test)

# high iron
fe2_temp3vstemp05_decide <- decideTestsDGE(fe2_temp3vstemp05_test)
fe2_temp3vstemp05_info <- out_toptag(fe2_temp3vstemp05_test)

fe2_temp6vstemp05_decide <- decideTestsDGE(fe2_temp6vstemp05_test)
fe2_temp6vstemp05_info <- out_toptag(fe2_temp6vstemp05_test)

fe2_temp6vstemp3_decide <- decideTestsDGE(fe2_temp6vstemp3_test)
fe2_temp6vstemp3_info <- out_toptag(fe2_temp6vstemp3_test)

# fe pairwise comparisons
temp05_fe2vsfe0_decide <- decideTestsDGE(temp05_fe2vsfe0_test)
temp05_fe2vsfe0_info <- out_toptag(temp05_fe2vsfe0_test)

temp3_fe2vsfe0_decide <- decideTestsDGE(temp3_fe2vsfe0_test)
temp3_fe2vsfe0_info <- out_toptag(temp3_fe2vsfe0_test)

temp6_fe2vsfe0_decide <- decideTestsDGE(temp6_fe2vsfe0_test)
temp6_fe2vsfe0_info <- out_toptag(temp6_fe2vsfe0_test)


# Fe decision
fe_all_decide <- decideTestsDGE(fe_all_test)
fe_all_info <- out_toptag(fe_all_test)

# Temperature decision
temp_3vstemp_05_decide <- decideTestsDGE(temp_3vstemp_05_test)
temp_3vstemp_05_info <- out_toptag(temp_3vstemp_05_test)

temp_6vstemp_05_decide <- decideTestsDGE(temp_6vstemp_05_test)
temp_6vstemp_05_info <- out_toptag(temp_6vstemp_05_test)

temp_6vstemp_3_decide <- decideTestsDGE(temp_6vstemp_3_test)
temp_6vstemp_3_info <- out_toptag(temp_6vstemp_3_test)

# interaction DE decision
tempfe_interaction_temp_6vstemp_05_decide <- decideTestsDGE(tempfe_interaction_temp_6vstemp_05_test)
tempfe_interaction_temp_6vstemp_05_info <- out_toptag(tempfe_interaction_temp_6vstemp_05_test)

tempfe_interaction_temp_3vstemp_05_decide <- decideTestsDGE(tempfe_interaction_temp_3vstemp_05_test)
tempfe_interaction_temp_3vstemp_05_info <- out_toptag(tempfe_interaction_temp_3vstemp_05_test)

# checking total number of DE genes for each test
# fe_all_decide@.Data %>% abs() %>% sum()
# 
# temp_3vstemp_05_decide@.Data %>% abs() %>% sum()
# temp_6vstemp_05_decide@.Data %>% abs() %>% sum()
# temp_6vstemp_3_decide@.Data %>% abs() %>% sum()
# 
# tempfe_interaction_temp_3vstemp_05_decide@.Data %>% abs() %>% sum()
# tempfe_interaction_temp_6vstemp_05_decide@.Data %>% abs() %>% sum()

# getting only the rows which do not have na for grpnorm_taxgrp
#grp_norm_no_grp_na <- grp_norm[!is.na(grp_norm$grpnorm_taxgrp),]
#LJ changed the line above into this because I didn't have any NA values
grp_norm_no_grp_na <- dplyr::filter (grp_norm, !grepl('^$', grpnorm_taxgrp))

# dataframe with actual values used for DE analysis (transformed vals by tax-specific scaling factor)
grp_norm_sub_scaled_copy <- grp_norm_sub_scaled
names(grp_norm_sub_scaled_copy) <- paste0(names(grp_norm_sub_scaled_copy), '-scaled')


## making compiled dataframe
tfg_de_df <- data.frame(Fe_vs_noFe_de = fe_all_decide %>% as.vector(),
                        Fe_vs_noFe_pvalue = fe_all_info$fe_all_pvalue,
                        Fe_vs_noFe_foldchange = fe_all_info$fe_all_foldchange,
                        # temperature comparisons
                        temp_3C_vs_0C_de = temp_3vstemp_05_decide %>% as.vector(),
                        temp_3C_vs_0C_pvalue = temp_3vstemp_05_info$temp_3vstemp_05_pvalue,
                        temp_3C_vs_0C_foldchange = temp_3vstemp_05_info$temp_3vstemp_05_foldchange,
                        
                        temp_6C_vs_0C_de = temp_6vstemp_05_decide %>% as.vector(),
                        temp_6C_vs_0C_pvalue = temp_6vstemp_05_info$temp_6vstemp_05_pvalue,
                        temp_6C_vs_0C_foldchange = temp_6vstemp_05_info$temp_6vstemp_05_foldchange,
                        
                        temp_6C_vs_3C_de = temp_6vstemp_3_decide %>% as.vector(),
                        temp_6C_vs_3C_pvalue = temp_6vstemp_3_info$temp_6vstemp_3_pvalue,
                        temp_6C_vs_3C_foldchange = temp_6vstemp_3_info$temp_6vstemp_3_foldchange,
                        
                        int_fe_temp_6C_vs_0C_de = tempfe_interaction_temp_6vstemp_05_decide %>% as.vector(),
                        int_fe_temp_6C_vs_0C_pvalue = tempfe_interaction_temp_6vstemp_05_info$tempfe_interaction_temp_6vstemp_05_pvalue,
                        int_fe_temp_6C_vs_0C_foldchange = tempfe_interaction_temp_6vstemp_05_info$tempfe_interaction_temp_6vstemp_05_foldchange,
                        int_fe_temp_3C_vs_0C_de = tempfe_interaction_temp_3vstemp_05_decide %>% as.vector(),
                        int_fe_temp_3C_vs_0C_pvalue = tempfe_interaction_temp_3vstemp_05_info$tempfe_interaction_temp_3vstemp_05_pvalue,
                        int_fe_temp_3C_vs_0C_foldchange = tempfe_interaction_temp_3vstemp_05_info$tempfe_interaction_temp_3vstemp_05_foldchange,
                        # pairwise comparisons
                        noFe_3C_vs_noFe_0C_de = fe0_temp3vstemp05_decide %>% as.vector(),
                        noFe_3C_vs_noFe_0C_pvalue = fe0_temp3vstemp05_info$fe0_temp3vstemp05_pvalue,
                        noFe_3C_vs_noFe_0C_foldchange = fe0_temp3vstemp05_info$fe0_temp3vstemp05_foldchange,
                        
                        Fe_3C_vs_Fe_0C_de = fe2_temp3vstemp05_decide %>% as.vector(),
                        Fe_3C_vs_Fe_0C_pvalue = fe2_temp3vstemp05_info$fe2_temp3vstemp05_pvalue,
                        Fe_3C_vs_Fe_0C_foldchange = fe2_temp3vstemp05_info$fe2_temp3vstemp05_foldchange,
                        
                        noFe_6C_vs_noFe_0C_de = fe0_temp6vstemp05_decide %>% as.vector(),
                        noFe_6C_vs_noFe_0C_pvalue = fe0_temp6vstemp05_info$fe0_temp6vstemp05_pvalue,
                        noFe_6C_vs_noFe_0C_foldchange = fe0_temp6vstemp05_info$fe0_temp6vstemp05_foldchange,
                        
                        Fe_6C_vs_Fe_0C_de = fe2_temp6vstemp05_decide %>% as.vector(),
                        Fe_6C_vs_Fe_0C_pvalue = fe2_temp6vstemp05_info$fe2_temp6vstemp05_pvalue,
                        Fe_6C_vs_Fe_0C_foldchange = fe2_temp6vstemp05_info$fe2_temp6vstemp05_foldchange,
                        
                        noFe_6C_vs_noFe_3C_de = fe0_temp6vstemp3_decide %>% as.vector(),
                        noFe_6C_vs_noFe_3C_pvalue = fe0_temp6vstemp3_info$fe0_temp6vstemp3_pvalue,
                        noFe_6C_vs_noFe_3C_foldchange = fe0_temp6vstemp3_info$fe0_temp6vstemp3_foldchange,
                        
                        Fe_6C_vs_Fe_3C_de = fe2_temp6vstemp3_decide %>% as.vector(),
                        Fe_6C_vs_Fe_3C_pvalue = fe2_temp6vstemp3_info$fe2_temp6vstemp3_pvalue,
                        Fe_6C_vs_Fe_3C_foldchange = fe2_temp6vstemp3_info$fe2_temp6vstemp3_foldchange,
                        
                        Fe_0C_vs_noFe_0C_de = temp05_fe2vsfe0_decide %>% as.vector(),
                        Fe_0C_vs_noFe_0C_pvalue = temp05_fe2vsfe0_info$temp05_fe2vsfe0_pvalue,
                        Fe_0C_vs_noFe_0C_foldchange = temp05_fe2vsfe0_info$temp05_fe2vsfe0_foldchange,
                        
                        Fe_3C_vs_noFe_3C_de = temp3_fe2vsfe0_decide %>% as.vector(),
                        Fe_3C_vs_noFe_3C_pvalue = temp3_fe2vsfe0_info$temp3_fe2vsfe0_pvalue,
                        Fe_3C_vs_noFe_3C_foldchange = temp3_fe2vsfe0_info$temp3_fe2vsfe0_foldchange,
                        
                        Fe_6C_vs_noFe_6C_de = temp6_fe2vsfe0_decide %>% as.vector(),
                        Fe_6C_vs_noFe_6C_pvalue = temp6_fe2vsfe0_info$temp6_fe2vsfe0_pvalue,
                        Fe_6C_vs_noFe_6C_foldchange = temp6_fe2vsfe0_info$temp6_fe2vsfe0_foldchange,
                        
                        # including dataframe with annotations
                        grp_norm_no_grp_na,
                        # including scaled values,
                        grp_norm_sub_scaled_copy,
                        # columns for plotting DE analysis
                        fe_neg_de = ifelse(fe_all_decide == -1, yes = -1, no = 0) %>% as.vector(),
                        fe_pos_de = ifelse(fe_all_decide == 1, yes = 1, no = 0) %>% as.vector(),
                        # temperature is a bit harder than Fe. The gene has to be DE in either comparisons from 3vs0.5 *or* 6vs0.5 to be counted in this column. But if the gene is negatively expressed from 3vs0.5, and positively expressed from 6vs0.5, it would be excluded.
                        temp_neg_de = ifelse((temp_3vstemp_05_decide == -1 | temp_6vstemp_05_decide == -1) &
                                               (temp_3vstemp_05_decide != 1 | temp_6vstemp_05_decide != 1), 
                                             yes = -1, no = 0) %>% as.vector(),
                        
                        temp_pos_de = ifelse((temp_3vstemp_05_decide == 1 | temp_6vstemp_05_decide == 1) &
                                               (temp_3vstemp_05_decide != -1 | temp_6vstemp_05_decide != -1), 
                                             yes = 1, no = 0) %>% as.vector(),
                        tempfe_neg_de = ifelse((tempfe_interaction_temp_6vstemp_05_decide == -1 | tempfe_interaction_temp_3vstemp_05_decide == -1) & 
                                                 (tempfe_interaction_temp_6vstemp_05_decide != 1 & tempfe_interaction_temp_3vstemp_05_decide != 1), 
                                               yes = -1, no = 0) %>% as.vector(),
                        tempfe_pos_de = ifelse((tempfe_interaction_temp_6vstemp_05_decide == 1 | tempfe_interaction_temp_3vstemp_05_decide == 1) & 
                                                 (tempfe_interaction_temp_6vstemp_05_decide != -1 & tempfe_interaction_temp_3vstemp_05_decide != -1), 
                                               yes = 1, no = 0) %>% as.vector(),
                        # abs_fe is for ranking taxonomic groups for plotting later on.
                        abs_fe = ifelse(abs(fe_all_decide) == 1, yes = 1, no = 0) %>% as.vector())


# old analysis

# fe0_temp05vs3 <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'fe0_temp05vstemp3'])
# fe0_temp05vs6 <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'fe0_temp05vstemp6'])
# fe0_temp3vs6 <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'fe0_temp3vstemp6'])
# 
# fe2_temp05vs3 <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'fe2_temp05vstemp3'])
# fe2_temp05vs6 <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'fe2_temp05vstemp6'])
# fe2_temp3vs6 <- glmQLFTest(tfg_fit_factors, contrast = my_contrasts[,'fe2_temp3vstemp6'])
# 
# 
# fe0_temp05vs3_decide <- decideTestsDGE(object = fe0_temp05vs3)
# fe0_temp05vs6_decide <- decideTestsDGE(object = fe0_temp05vs6)
# fe0_temp3vs6_decide <- decideTestsDGE(object = fe0_temp3vs6)
# 
# fe2_temp05vs3_decide <- decideTestsDGE(object = fe2_temp05vs3)
# fe2_temp05vs6_decide <- decideTestsDGE(object = fe2_temp05vs6)
# fe2_temp3vs6_decide <- decideTestsDGE(object = fe2_temp3vs6)
# 
# 
# fe2_temp05vs3_decide@.Data %>% abs() %>% sum()
# fe2_temp05vs6_decide@.Data %>% abs() %>% sum()
# fe0_temp05vs3_decide@.Data %>% abs() %>% sum()
# fe0_temp05vs6_decide@.Data %>% abs() %>% sum()


## getting significance tests of old model formulation

tfg_inter <- glmQLFTest(tfg_fit_contin, coef = 1)
tfg_fe <- glmQLFTest(tfg_fit_contin, coef = 2)
tfg_temp3 <- glmQLFTest(tfg_fit_contin, coef = 3)
tfg_temp6 <- glmQLFTest(tfg_fit_contin, coef = 4)
tfg_fetemp3 <- glmQLFTest(tfg_fit_contin, coef = 5)
tfg_fetemp6 <- glmQLFTest(tfg_fit_contin, coef = 6)
# 
# tfg_inter_decide <- decideTestsDGE(tfg_inter)
# tfg_fe_decide <- decideTestsDGE(tfg_fe)
# tfg_temp3_decide <- decideTestsDGE(tfg_temp3)
# tfg_temp6_decide <- decideTestsDGE(tfg_temp6)
# tfg_fetemp3_decide <- decideTestsDGE(tfg_fetemp3)
# tfg_fetemp6_decide <- decideTestsDGE(tfg_fetemp6)


### old version writing final dataframe

# tfg_de_df <- data.frame(tfg_inter_decide, tfg_fe_decide, tfg_temp3_decide,
#            tfg_temp6_decide, tfg_fetemp3_decide, tfg_fetemp6_decide,
#            grp_norm_sub_not0,
#           abs_fe = ifelse(abs(tfg_fe_decide) == 1, yes = 1, no = 0) %>% as.vector(),
#           abs_temp = ifelse(abs(tfg_temp3_decide) == 1 | abs(tfg_temp6_decide) == 1, yes = 1, no = 0) %>% as.vector(),
#           abs_tempfe = ifelse(abs(tfg_fetemp3_decide) == 1 | abs(tfg_fetemp6_decide) == 1, yes = 1, no = 0) %>% as.vector(),
#           # summaries of negative DE by treatment *group*
#           fe_neg_de = ifelse((tfg_fe_decide == -1 | tfg_fe_decide == -1) & 
#                                  (tfg_fe_decide != 1 | tfg_fe_decide!= 1), yes = -1, no = 0) %>% as.vector(),
#           temp_neg_de = ifelse((tfg_temp3_decide == -1 | tfg_temp6_decide == -1) & 
#                                  (tfg_temp3_decide != 1 | tfg_temp6_decide!= 1), yes = -1, no = 0) %>% as.vector(),
#           fetemp_neg_de = ifelse((tfg_fetemp3_decide == -1 | tfg_fetemp6_decide == -1) & 
#                                  (tfg_fetemp3_decide != 1 | tfg_fetemp6_decide!= 1), yes = -1, no = 0) %>% as.vector(),
#           
#           # summaries of positive DE by treatment *group*
#           fe_pos_de = ifelse((tfg_fe_decide == 1 | tfg_fe_decide == 1) & 
#                                  (tfg_fe_decide != -1 | tfg_fe_decide!= -1), yes = 1, no = 0) %>% as.vector(),
#           temp_pos_de = ifelse((tfg_temp3_decide == 1 | tfg_temp6_decide == 1) & 
#                                  (tfg_temp3_decide != -1 | tfg_temp6_decide!= -1), yes = 1, no = 0) %>% as.vector(),
#           fetemp_pos_de = ifelse((tfg_fetemp3_decide == 1 | tfg_fetemp6_decide == 1) & 
#                                  (tfg_fetemp3_decide != -1 | tfg_fetemp6_decide!= -1), yes = 1, no = 0) %>% as.vector())

write.csv(tfg_de_df, file = 'LJ_MCL_tfg_de_annotation_allTFG.grpnorm_mmetsp_fc_pn_reclassified.csv')

### plotting de results by taxon

# blank plot just used for the axis.
blank_p <- tfg_de_df %>% 
  group_by(grpnorm_taxgrp) %>% 
  summarise(fe_grp = sum(abs_fe)) %>%  
  ggplot(aes(x = fe_grp, 
             y = fct_reorder(grpnorm_taxgrp, fe_grp))) + 
  # geom_point() +
  theme_bw() +
  xlab('') + 
  theme(axis.text.x = element_text(colour = 'white'), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), 
        axis.ticks = element_blank()) +
  ylab('');blank_p

# old plots that show total DE, not subdivided into positive and negative DE:

# p1_fe <- tfg_de_df %>% 
#   group_by(grpnorm_taxgrp) %>% 
#   summarise(fe_grp = sum(abs_fe)) %>%  
#   ggplot(aes(x = fe_grp, y = fct_reorder(grpnorm_taxgrp, fe_grp))) + 
#   geom_point(size = 3, alpha = 0.5) +
#   theme_bw() +
#   xlab('Iron-Related DE') + 
#   ylab('') +
#   theme(axis.text.y = element_blank());p1_fe
# 
# p2_temp <- tfg_de_df %>% 
#   group_by(grpnorm_taxgrp) %>% 
#   summarise(fe_grp = sum(abs_fe),
#             temp_grp = sum(abs_temp)) %>% 
#   ggplot(aes(x = temp_grp, y = fct_reorder(grpnorm_taxgrp, fe_grp))) + 
#   geom_point(size = 3, alpha = 0.5) +
#   theme_bw() + 
#   ylab('') + 
#   xlab('Temperature-Related DE') +
#   theme(axis.text.y = element_blank());p2_temp
# 
# p3_tempfe <- tfg_de_df %>% 
#   group_by(grpnorm_taxgrp) %>% 
#   summarise(fe_grp = sum(abs_fe),
#             fetemp_grp = sum(abs_tempfe)) %>% 
#   ggplot(aes(x = fetemp_grp, y = fct_reorder(grpnorm_taxgrp, fe_grp))) + 
#   geom_point(size = 3, alpha = 0.5) +
#   theme_bw() + 
#   ylab('') + 
#   xlab('Temperature-Iron-Related DE') +
#   theme(axis.text.y = element_blank());p3_tempfe
# 
# 
# grid.arrange(blank_p, p1_fe, p2_temp, p3_tempfe, nrow = 1)


# DE plot with positive and negative differences shown.

pn_fe <- tfg_de_df %>% 
  group_by(grpnorm_taxgrp) %>% 
  summarise(fe_grp = sum(abs_fe),
            fe_grp_neg = sum(fe_neg_de),
            fe_grp_pos = sum(fe_pos_de)) %>%  
  ggplot(aes(x = fe_grp_neg, y = fct_reorder(grpnorm_taxgrp, fe_grp))) + 
  geom_point(size = 3, alpha = 0.5, colour = 'blue') +
  geom_point(aes(x = fe_grp_pos, y = fct_reorder(grpnorm_taxgrp, fe_grp)),
             size = 3, alpha = 0.5, colour = 'red') +
  geom_segment(aes(x = fe_grp_neg, xend = 0, 
                   y = fct_reorder(grpnorm_taxgrp, fe_grp), 
                   yend = fct_reorder(grpnorm_taxgrp, fe_grp)), colour = 'blue', alpha = 0.2, lwd = 2) +
  geom_segment(aes(x = fe_grp_pos, xend = 0, 
                   y = fct_reorder(grpnorm_taxgrp, fe_grp), 
                   yend = fct_reorder(grpnorm_taxgrp, fe_grp)), colour = 'red', alpha = 0.2, lwd = 2) +
  theme_bw() +
  xlab('Iron-Related DE') + 
  ylab('') +
  theme(axis.text.y = element_blank());pn_fe


pn_temp <- tfg_de_df %>% 
  group_by(grpnorm_taxgrp) %>% 
  summarise(fe_grp = sum(abs_fe),
            temp_grp_neg = sum(temp_neg_de),
            temp_grp_pos = sum(temp_pos_de)) %>%  
  ggplot(aes(x = temp_grp_neg, y = fct_reorder(grpnorm_taxgrp, fe_grp))) + 
  geom_point(size = 3, alpha = 0.5, colour = 'blue') +
  geom_point(aes(x = temp_grp_pos, y = fct_reorder(grpnorm_taxgrp, fe_grp)),
             size = 3, alpha = 0.5, colour = 'red') +
  geom_segment(aes(x = temp_grp_neg, xend = 0, 
                   y = fct_reorder(grpnorm_taxgrp, fe_grp), 
                   yend = fct_reorder(grpnorm_taxgrp, fe_grp)), colour = 'blue', alpha = 0.2, lwd = 2) +
  geom_segment(aes(x = temp_grp_pos, xend = 0, 
                   y = fct_reorder(grpnorm_taxgrp, fe_grp), 
                   yend = fct_reorder(grpnorm_taxgrp, fe_grp)), colour = 'red', alpha = 0.2, lwd = 2) +
  theme_bw() +
  xlab('Temperature-Related DE') + 
  ylab('') +
  theme(axis.text.y = element_blank());pn_temp

pn_tempfe <- tfg_de_df %>% 
  group_by(grpnorm_taxgrp) %>% 
  summarise(fe_grp = sum(abs_fe),
            fetemp_grp_neg = sum(tempfe_neg_de),
            fetemp_grp_pos = sum(tempfe_pos_de)) %>%  
  ggplot(aes(x = fetemp_grp_neg, y = fct_reorder(grpnorm_taxgrp, fe_grp))) + 
  geom_point(size = 3, alpha = 0.5, colour = 'blue') +
  geom_point(aes(x = fetemp_grp_pos, y = fct_reorder(grpnorm_taxgrp, fe_grp)),
             size = 3, alpha = 0.5, colour = 'red') +
  geom_segment(aes(x = fetemp_grp_neg, xend = 0, 
                   y = fct_reorder(grpnorm_taxgrp, fe_grp), 
                   yend = fct_reorder(grpnorm_taxgrp, fe_grp)), colour = 'blue', alpha = 0.2, lwd = 2) +
  geom_segment(aes(x = fetemp_grp_pos, xend = 0, 
                   y = fct_reorder(grpnorm_taxgrp, fe_grp), 
                   yend = fct_reorder(grpnorm_taxgrp, fe_grp)), colour = 'red', alpha = 0.2, lwd = 2) +
  theme_bw() +
  xlab('Iron-Temperature-Related DE') + 
  ylab('') +
  theme(axis.text.y = element_blank());pn_tempfe

grid.arrange(blank_p, pn_fe, pn_temp, pn_tempfe, nrow = 1)













#MCL1 <- read.csv("MCL_tfg_de_annotation_allTFG.grpnorm_mmetsp_fc_pn_reclassified.csv")
#MCL2 <- read.csv("LJ_MCL_tfg_de_annotation_allTFG.grpnorm_mmetsp_fc_pn_reclassified.csv")

#all(MCL1$Fe_3C_vs_noFe_3C_foldchange == MCL2$Fe_3C_vs_noFe_3C_foldchange)

