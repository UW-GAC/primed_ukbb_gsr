# library packages
library(AnVIL, quietly = TRUE)
library(AnvilDataModels, quietly = TRUE)
library(R.utils, quietly = TRUE) # unzip the .bgz file
library(curl, quietly = TRUE) # read from public AWS
library(data.table, quietly = TRUE) # fread from publicly available AWS
library(tidyverse, quietly = TRUE, warn.conflicts = FALSE) # write_tsv function at the end
library(tools, quietly = TRUE) # get MD5 hash of the data file
library(argparser, quietly = TRUE) # execute script from console


# identify which phenotype you want to analyze
p <- arg_parser("Pan UKBB GSR Download")

p <- add_argument(parser = p,
                  arg = "--phenocode",
                  type = "character",
                  nargs = Inf,
                  help = "Identify the unique code of the phenotype(s) as described in the Pan-UK Biobank phenotype manifest.")

p <- add_argument(parser = p,
                  arg = "--population",
                  type = "character",
                  default = "all_available",
                  nargs = Inf,
                  help = "Identify the specific population(s) for which you want the Pan-UK Biobank GSR data.
                          If you list multiple phenotypes, then data will be restricted to these populations for every listed phenotype.
                          Confirm that the population(s) of interest have data for your phenotype(s) of interest using the Pan-UK Biobank phenotype manifest.")

argv <- parse_args(parser = p)
phenocode_list <- argv$phenocode
population_list <- argv$population

# set global timeout options
hour_timeout <- 1
options(timeout = 60*60*hour_timeout)
rm(list = c("hour_timeout"))


# write function for pulling compressed data directly from publicly available AWS
read_AWS <- function(url, save_location){
  # index the system time when function started
  start_time <- Sys.time()
  
  # identify the file name
  file_name <- gsub("\\..*", "", basename(url))
  
  # label file path
  file_path <- paste0(file_name, ".gz", collapse = "")
  
  # download compressed .tsv file from Pan-UK Biobank's file on AWS
  curl_download(url = url,
                destfile = file_path,
                mode = "wb",
                quiet = FALSE)
    
  # STATUS
  message("Downloaded compressed .tsv file.")
  
  # read in the decompressed file as tsv
  data <- fread(file = file_path,
                sep = "\t",
                header = TRUE,
                nrows = Inf,
                stringsAsFactors = FALSE)
  
  # STATUS
  message("Read in .tsv file.")
  
  # STATUS
  message("Completed in ", round(difftime(Sys.time(), start_time, units = "mins"), 0), " minutes.")
  
  return(data)
}


# prepare functions to extract key information from phenotype manifest
which_trait_type <- function(x){
  if (x %in% c("biomarkers", "continuous")) {
    return("quantitative")
  } else if (x %in% c("prescriptions", "icd10", "phecode")) {
    return("binary")
  } else if (x %in% c("categorical")) {
    return("categorical")
  }
}


which_pops <- function(x){
  pop_names <- c("EUR" = "European ancestry",
                 "CSA" = "Central/South Asian ancestry",
                 "AFR" = "African ancestry",
                 "EAS" = "East Asian ancestry",
                 "MID" = "Middle Eastern ancestry",
                 "AMR" = "Admixed American ancestry")
  if (any(x == names(pop_names))) {
    unname(pop_names[which(x == names(pop_names))])
  } else if (x == "meta") {
    pops <- unlist(strsplit(phenotype_info$pops, split = ","))
    pops <- sapply(pops, function(x){unname(pop_names[x == names(pop_names)])})
    paste(pops, collapse = " | ")
  } else if (x == "metaHQ"){
    pops <- unlist(strsplit(phenotype_info$pops_pass_qc, split = ","))
    pops <- sapply(pops, function(x){unname(pop_names[x == names(pop_names)])})
    paste(pops, collapse = " | ")
  }
}


which_sampsize <- function(x){
  pop_names <- c("EUR" = "European ancestry",
                 "CSA" = "Central/South Asian ancestry",
                 "AFR" = "African ancestry",
                 "EAS" = "East Asian ancestry",
                 "MID" = "Middle Eastern ancestry",
                 "AMR" = "Admixed American ancestry")
  if (any(x == names(pop_names))) {
    
    n_controls <- unname(unlist(c(phenotype_info)[paste0("n_controls_", x)]))
    n_cases <- unname(unlist(c(phenotype_info)[paste0("n_cases_", x)]))
    
  } else if (x == "meta") {
    
    pops_in_meta <- str_split(unlist(c(phenotype_info)["pops"]), ",")
    n_controls <- sum(unlist(sapply(pops_in_meta,
                                    function(y){c(phenotype_info)[paste0("n_controls_", y)]})))
    n_cases <- sum(unlist(sapply(pops_in_meta,
                                 function(y){c(phenotype_info)[paste0("n_cases_", y)]})))
    if (n_cases != phenotype_info$n_cases_full_cohort_both_sexes) {
      warning("The number of cases calculated does not match the number of cases stated.")
    }
    
  } else if (x == "metaHQ"){
    
    pops_in_metaHQ <- str_split(unlist(c(phenotype_info)["pops_pass_qc"]), ",")
    n_controls <- sum(unlist(sapply(pops_in_metaHQ,
                                    function(y){c(phenotype_info)[paste0("n_controls_", y)]})))
    n_cases <- sum(unlist(sapply(pops_in_metaHQ,
                                 function(y){c(phenotype_info)[paste0("n_cases_", y)]})))
    if (n_cases != phenotype_info$n_cases_hq_cohort_both_sexes) {
      warning("The number of cases calculated does not match the number of cases stated.")
    }
    
  } else {
    stop("Error: you identified a population that does not exist.")
  }
  
  # make sure NULL or NA values are set to 0 to not interfere with addition
  n_controls <- ifelse(is.numeric(n_controls), n_controls, 0)
  n_cases <- ifelse(is.numeric(n_cases), n_cases, 0)
  
  if (phenotype_is_binary) {
    return(list("n_ctrl" = n_controls,
                "n_case" = n_cases,
                "n_samp" = n_controls + n_cases,
                "n_effective" = 4 / (1 / n_cases + 1 / n_controls)))
  } else {
    return(list("n_ctrl" = NA,
                "n_case" = NA,
                "n_samp" = n_cases,
                "n_effective" = n_cases))
  }
}


which_transf <- function(x){
  common <- c("log" = "",
              "inverse-rank normal transformation" = "irnt",
              "winsorization" = "",
              "standardization" = "")
  if (is.na(x)) {
    "none"
  } else if (sum(x == common) == 1) {
    names(which(common == x))
  } else {
    x
  }
}


# download the Pan-UK Biobank phenotype manifest
url <- "https://pan-ukb-us-east-1.s3.amazonaws.com/"
url_manifest <- paste0(url, "sumstats_release/phenotype_manifest.tsv.bgz")
manifest <- read_AWS(url_manifest)
rm(list = c("url_manifest"))


# download the variant manifest file
(url_VMF <- "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/full_variant_qc_metrics.txt.bgz")
data_VMF <- read_AWS(url_VMF)
rm(list = c("url_VMF"))
data_VMF <- data_VMF[, c("chrom", "pos", "ref", "alt", "rsid", "info")]
head(data_VMF)


# identify if the phenotype is binary
for (input in phenocode_list) {
  phenotype_info <- manifest[manifest$phenocode %in% input, ]
  phenotype_name <- phenotype_info$description
  phenotype_is_binary <- which_trait_type(phenotype_info$trait_type) == "binary"
  
  # Label step of the download
  print(paste0("Beginning download of <", phenotype_name, "> dataset from Pan UKBB."))


  # download the GSR data pertaining to the chosen phenotype
  (url_PHE <- paste0(url, gsub("s3://pan-ukb-us-east-1/", "", phenotype_info$aws_path)))
  data_PHE <- read_AWS(url_PHE)
  rm(list = c("url_PHE"))
  print(head(data_PHE))
  
  
  # check for unique matching ID
  print(colnames(data_PHE)[1:4])
  print(colnames(data_VMF)[1:4])
  
  
  # merge the per-phenotype data with the variant manifest file to get RSID
  data_PHE <- merge(x = data_PHE, y = data_VMF[, c("chrom", "pos", "ref", "alt", "rsid", "info")],
                    by.x = c("chr",   "pos", "ref", "alt"),
                    by.y = c("chrom", "pos", "ref", "alt"),
                    all.x = TRUE, sort = FALSE)
  print(head(data_PHE))
  
  
  ############################
  # variable transformations #
  ############################
  
  # remove natural log transformation of p-value
  pval_cols <- colnames(data_PHE)[grepl("pval_", colnames(data_PHE))]
  data_PHE[, c(pval_cols) := lapply(.SD, function(x){exp(x)}), .SDcols = pval_cols]
  pval_log10 <- gsub("pval", "pval_log10", pval_cols)
  data_PHE[, c(pval_log10) := lapply(.SD, function(x){-log10(x)}), .SDcols = pval_cols]
  rm(list = c("pval_cols", "pval_log10"))
  
  
  # list names in dataset by population
  meta_pop <- c("af_meta", "af_cases_meta", "af_controls_meta", "beta_meta", "se_meta", "pval_meta", "pval_heterogeneity", "pval_log10_meta")
  metaHQ_pop <- c("af_meta_hq", "af_cases_meta_hq", "af_controls_meta_hq", "beta_meta_hq", "se_meta_hq", "pval_meta_hq", "pval_heterogeneity_hq", "pval_log10_meta_hq")
  key_meta <- c("effect_allele_freq", "eaf_case", "eaf_ctrl", "beta", "se", "p_value", "heterogeneity_p_value", "p_value_log10")
  
  
  # identify which ancestry populations have data for that phenotype
  pop_list <- unlist(str_split(unlist(c(phenotype_info)["pops"]), ","))
  all_pop <- lapply(pop_list, function(x){paste(c("af_", "af_cases_", "af_controls_", "beta_", "se_", "pval_", "pval_log10_", "low_confidence_"), x, sep = "")})
  key_all <- c("effect_allele_freq", "eaf_case", "eaf_ctrl", "beta", "se", "p_value", "p_value_log10", NA)
  
  
  # create a list of meta, high-quality meta, and specific ancestry populations
  pop_col_subset <- append(list(meta_pop, metaHQ_pop), all_pop)
  names(pop_col_subset) <- c("meta", "metaHQ", pop_list)
  
  
  # if only a specific population is of interest, then only download data for that population
  if (!identical(population_list, "all_available")) {
    pop_col_subset <- pop_col_subset[intersect(unlist(population_list), names(pop_col_subset))]
  }
  rm(list = c("meta_pop", "metaHQ_pop", "pop_list", "all_pop"))
  
  
  # create dataset for each population of interest
  for (i in 1:length(pop_col_subset)) {
    # make a copy of the dataset
    data_temp <- copy(data_PHE)
  
    # identify which population we are working with
    pop <- names(pop_col_subset)[i]
    
    # set the variable key according to the population
    if (grepl("meta", pop)) {
      key <- key_meta
    } else {
      key <- key_all
    }
    
    # identify which variables pertain to the dataset
    pop_vars <- pop_col_subset[[i]]
    
    # make a key that returns either the matching value or NA
    L0NA <- function(x){
      y <- pop_vars[key == x]
      return(ifelse(length(y) > 0, y, NA))
    }
    
    # rename the population subset of phenotype data to match PRIMED GSR DD
    rename <- c(
      SNPID = NA,
      chromosome = "chr",
      position = "pos",
      rsID = "rsid",
      strand = "strand",
      effect_allele = "alt",
      other_allele = "ref",
      ref_allele = "ref2",
      alt_allele = "alt2",
      effect_allele_freq = L0NA("effect_allele_freq"),
      eaf_case = L0NA("eaf_case"),
      eaf_ctrl = L0NA("eaf_ctrl"),
      mac = NA,
      p_value = L0NA("p_value"),
      p_value_log10 = L0NA("p_value_log10"),
      beta = L0NA("beta"),
      se = L0NA("se"),
      beta_ci_lower = NA,
      beta_ci_upper = NA,
      odds_ratio = NA,
      OR_ci_lower = NA,
      OR_ci_upper = NA,
      direction_of_effect = NA,
      n_samp = NA,
      n_case = NA,
      n_ctrl = NA,
      eff_sample_size = NA,
      is_imputed = NA,
      imputation_quality_score = "info",
      heterogeneity_p_value = L0NA("heterogeneity_p_value"),
      heterogeneity_I2 = NA
    )
    rm(list = c("pop_vars", "L0NA"))
    
    # create duplicate columns for ref and alt since it is used twice
    data_temp[, ':='(ref2 = ref,
                     alt2 = alt)]
  
    # rename the dataset to match PRIMED notation
    setnames(data_temp,
             old = unname(rename[!is.na(rename)]),
             new =  names(rename[!is.na(rename)]),
             skip_absent = TRUE)
    
    # subset to only columns in the population
    data_temp <- data_temp[, intersect(c(c("chromosome", "position", "effect_allele", "other_allele", "rsID", "SNPID"), key),
                                       colnames(data_temp)), with = FALSE]
    rm(list = "key")
    
    # remove all rows where the p-value is NA
    print(data_temp)
    data_temp <- subset(data_temp, !is.na(data_temp$p_value))
    print(data_temp)
    
    # construct 95% confidence intervals
    data_temp[, ':='(beta_ci_lower = beta + qnorm(0.025) * se,
                     beta_ci_upper = beta + qnorm(0.975) * se)]
    
    # if the phenotype is binary, compute the odds ratio and corresponding confidence intervals
    if (phenotype_is_binary) {
      data_temp[, ':='(odds_ratio  = exp(beta),
                       OR_ci_lower = exp(beta_ci_lower),
                       OR_ci_upper = exp(beta_ci_upper))]
    }
    
    # create new column for DNA strand
    data_temp[, ':='(strand = "forward")]
    
    # create new column for SNPID
    data_temp[, ':='(SNPID = paste(chromosome, position, other_allele, effect_allele, sep = ":"))]
    
    # direction of the beta estimate (code will not pass validation if exactly equal to 0)
    data_temp[, ':='(direction_of_effect = ifelse(beta < 0, "-",
                                                 ifelse(beta > 0, "+", "0")))]
    if (any(data_temp$direction_of_effect %in% "0")) {
      warning('Beta estimate is exactly equal to 0, so sign is neither "-" or "+"')
    }
    
    
    # put in the order of the data dictionary
    ddorder <- c(names(rename)[names(rename) %in% colnames(data_temp)], # order of columns in data dictionary
                 setdiff(colnames(data_temp), names(rename))) # columns not present in data dictionary
    setcolorder(data_temp, ddorder)
    rm(list = c("ddorder", "rename"))
    
    
    # check if any required columns are missing
    required <- c("chromosome", "position", "strand", "effect_allele", "other_allele",
                  "effect_allele_freq", "p_value", "beta", "se", "is_imputed")
    if (any(!(required %in% colnames(data_temp)))) {
      warning(paste0("The following required mappings are missing in the provided dataset:\n       ",
                  paste0(required[!(required %in% colnames(data_temp))], collapse = "\n       ")))
    }
    rm(list = c("required"))
    
    
    ###############################################################
    ### run the analysis and file workflows for this population ###
    ###############################################################
    
    fields <- list(
      "gsr_source" = "Pan UKBB",
      "gsr_source_url" = "https://pan.ukbb.broadinstitute.org/",
      "gwas_catalog_study_id" = NA,
      "dbgap_analysis_accession" = NA,
      "pubmed_id" = NA,
      "first_author" = NA,
      "link" = NA,
      "release_date" = NA,
      "consent_code" = "NRES",
      "upload_date" = format(Sys.Date(), "%m/%d/%Y"),
      "contributor_contact" = "primedconsortium@uw.edu",
      "trait" = phenotype_name,
      "trait_type" = which_trait_type(phenotype_info$trait_type),
      "trait_unit" = ifelse(phenotype_is_binary, "binary", "unknown"),
      "trait_transformation" = which_transf(phenotype_info$modifier),
      "trait_definition" = ifelse(!is.na(phenotype_info$description_more),
                                  phenotype_info$description_more,
                                  phenotype_info$description), ##### Longer definition not always given. #####
      "covariates" = "Age | Sex | Age * Sex | Age^2 | Age^2 * Sex | First 10 principal components",
      "concept_id" = NA, # We may search at this link, but need a formal rule for which to use: https://athena.ohdsi.org/search-terms/start
      "mapped_trait" = NA,
      "reference_assembly" = "GRCh37", # found on https://pan.ukbb.broadinstitute.org/docs/per-phenotype-files/index.html
      "dbsnp_build_version" = NA,
      "n_variants" = nrow(data_temp), # use total number of variants with p-values
      "min_MAF_filter" = NA,
      "min_MAC_filter" = 20,
      "genotyping_technology" = "genome-wide array",
      "genotyping_platform" = "UK Biobank Axiom Array",
      "is_imputed" = "TRUE",
      "imputation_reference_panel" = "Other", # details provided at https://www.nature.com/articles/s41586-018-0579-z#Sec5
      "imputation_reference_panel_detail" = "Primary imputation using HRC combined with imputation using UK10K + 1000G; see https://www.nature.com/articles/s41586-018-0579-z",
      "imputation_quality_filter" = 0.8, # minimum info score accepted by the Broad Institute
      "n_samp" = which_sampsize(pop)$n_samp,
      "n_case" = which_sampsize(pop)$n_case,
      "n_ctrl" = which_sampsize(pop)$n_ctrl,
      "n_effective" = which_sampsize(pop)$n_effective,
      "proportion_female" = NA,
      "age_mean" = NA,
      "age_median" = NA,
      "age_min" = NA,
      "age_max" = NA,
      "case_age_mean" = NA,
      "case_age_median" = NA,
      "case_age_min" = NA,
      "case_age_max" = NA,
      "ctrl_age_mean" = NA,
      "ctrl_age_median" = NA,
      "ctrl_age_min" = NA,
      "ctrl_age_max" = NA,
      "cohorts" = "UKBB",
      "is_meta_analysis" = ifelse(grepl("meta", pop), "TRUE", "FALSE"),
      "population_descriptors" = which_pops(pop),
      "population_proportions" = NA,
      "countries_of_recruitment" = "England | Scotland | Wales", # I believe that Northern Ireland is not included in UK Biobank
      "countries_of_birth" = NA,
      "analysis_method" = ifelse(phenotype_is_binary, "logistic mixed model", "linear mixed model"), # https://pan.ukbb.broadinstitute.org/docs/technical-overview
      "analysis_software" = "SAIGE implemented in Hail Batch")
    
    analysis <- tibble(field = names(fields),
                       value = unlist(fields))
    
    # remove missing, non-required fields
    required <- c("gsr_source", "consent_code", "upload_date", "contributor_contact", "trait", "trait_type",
                  "trait_unit", "trait_transformation", "trait_definition", "covariates", "concept_id",
                  "reference_assembly", "n_variants", "genotyping_technology", "genotyping_platform",
                  "is_imputed", "n_samp", "n_effective", "cohorts",  "is_meta_analysis", "population_descriptors",
                  "countries_of_recruitment", "analysis_method")
    if (analysis$value[analysis$field == "is_imputed"] %in% "TRUE") {
      required <- c(required, "imputation_reference_panel", "imputation_reference_panel_detail", "imputation_quality_filter")
    }
    if (analysis$value[analysis$field == "trait_type"] %in% "binary") {
      required <- c(required, "n_case", "n_ctrl")
    }
    analysis <- analysis[!is.na(analysis$value) | c(analysis$field %in% required), ]
    
    
    # save the analysis table
    outfile2 <- paste0(gsub(" ", "", phenotype_name), "_", pop, "_analysis.tsv")
    fwrite(analysis, outfile2, sep = "\t")
    
    
    # save the dataset and file tables split by chromosome
    setkey(data_temp, chromosome)
    
    print(unique(data_temp$chromosome))
    print(head(data_temp))
    
    for (chr in unique(data_temp$chromosome)) {
      # save the wrangled data
      outfile1 <- paste0(gsub(" ", "", phenotype_name), "_", pop, "_", chr, "_data.tsv.gz")
      fwrite(data_temp[as.character(chr)], outfile1, sep = "\t") # save to the local directory

      
      # save the file table
      file_table <- tibble(md5sum = md5sum(files = outfile1),
                           file_path = file.path(outfile1),
                           file_type = "data",
                           n_variants = nrow(data_temp[as.character(chr)]), # adjusted to exclude missing p-values
                           chromosome = chr)
      outfile3 <- paste0(gsub(" ", "", phenotype_name), "_", pop, "_", chr, "_file.tsv")
      fwrite(file_table, outfile3, sep = "\t")
    }
    
    # clear out the things that change with the dataset
    rm(list = c("outfile1", "outfile2", "outfile3",
                "pop", "data_temp", "fields", "analysis", "file_table"))
  }
  
  rm(list = c("key_all", "key_meta", "pop_col_subset", "i", "phenotype_name", "phenotype_is_binary", "bucket"))
}
